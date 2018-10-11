
#'   
#' ***  
#' ### Description:  
#' Determine if glycogen synthase 1 gene (GYS1) is differentially expressed 
#' between NN, PN and PP genotypes.  
#'   
#' ***  

#+  echo=FALSE
rm(list=ls())
library(prettydoc)


#+  echo=FALSE
dir="/mnt/research/NMDL/qPCR/GYS1"
qpcr <- read.table(paste(dir, "GYS1_qpcr.txt", sep="/"), 
    header=TRUE)

### RT-qPCR Data
qpcr

#+  echo=FALSE
geoMean <- function(x){
        prod(x[!is.na(x)])^(1/length(x[!is.na(x)]))
}

#+  echo=FALSE
meanCT <- function(gene, nm, dataT, byRow=FALSE){
    anim <- as.character(unique(dataT$Sample))

    if(byRow == TRUE){
    mCT <- unlist(lapply(as.character(anim), function(x) 
        mean(unlist(dataT[dataT$Sample %in% x, gene]), na.rm=TRUE)))

    } else {    
    
    mCT <- unlist(lapply(as.character(anim), function(x) 
        mean(dataT[dataT$Sample %in% x & dataT$Target %in% nm, gene], na.rm=TRUE)))
    }
    
    names(mCT) <- anim
    return(mCT)
}


#+  echo=FALSE
deltaCT <- function(gene, nm, control, dataT, byRow=FALSE){
    # Animal IDs
    anim <- as.character(unique(dataT$Sample))

    if(byRow == TRUE){
        # Geometric mean of control genes
        ref <- unlist(lapply(anim, function(x) 
            geoMean(unlist(dataT[dataT$Sample %in% x, control]))))
        names(ref) <- anim
        # Average of target gene
        MG <- meanCT(gene=gene, nm=nm, dataT=dataT, byRow=byRow)

        } else {

        # Geometric mean of control genes
        ref <- unlist(lapply(anim, function(x) 
            geoMean(dataT[dataT$Sample %in% x & dataT$Target %in% control, gene])))
        names(ref) <- anim 
        # Average of target gene
        MG <- meanCT(gene=gene, nm=nm, dataT=dataT, byRow=byRow)
    }

    # DeltaCT
    dCT <- MG - ref
    names(dCT) <- anim
    # Return results
    rst <- data.frame(Target=MG, Control=ref, dCT=dCT)
    return(rst)
}


#' ### Calculate deltaCT for target gene  
#' Results negative controls  
#+  echo=FALSE 
idx <- grep("NTC", qpcr$Sample)
qpcr[idx,]

#+  echo=FALSE
qpcr <- qpcr[-idx,]

#' Genotypes  
#+  echo=FALSE  
geno <- unique(qpcr[,c("Sample","Genotype")])
geno$Genotype <- factor(geno$Genotype)
qpcr <- qpcr[,-6]
table(geno$Genotype)

#' Control genes (three)  
#+  echo=FALSE  
control <- colnames(qpcr)[-c(1,5)]
control

#+  echo=FALSE
qpcr <- data.frame(apply(qpcr, 2, as.numeric))


#' Average DeltaCT GYS1  
#+  echo=FALSE
dCT <- deltaCT(gene="GYS1", nm="GYS1", control=control, dataT=qpcr, byRow=TRUE)
dCT <- data.frame(dCT, Genotype=geno$Genotype)
dCT


#' ### Test Significance  
#' Average dCT per Genotype  
#+  echo=FALSE  
avg.geno <- do.call(rbind, lapply(as.character(unique(dCT$Genotype)), function(x) 
 data.frame(AVG=mean(dCT[dCT$Genotype == x, "dCT"]), StDev=sd(dCT[dCT$Genotype == x, "dCT"]))))
rownames(avg.geno) <- as.character(unique(dCT$Genotype))
avg.geno

#' Simple linear regression  
#+  echo=FALSE  
anova(lm(dCT ~ Genotype, data=dCT))

#' Calculata delta-delta CT and fold change comparing normal (NN) with affected (PN and PP)  
#+  echo=FALSE  
fold <- do.call(rbind, lapply(c("PN","PP"), function(x) data.frame(
    DeltaDeltaCT=avg.geno[rownames(avg.geno) == x, "AVG"] - 
        avg.geno[rownames(avg.geno) == "NN", "AVG"],
    FoldChange=2^-(avg.geno[rownames(avg.geno) == x, "AVG"] - 
        avg.geno[rownames(avg.geno) == "NN", "AVG"]))))
rownames(fold) <- c("PN","PP")
fold


 
#+ echo = FALSE, eval = FALSE  
htmlRunR
Summary.R  nodes=1,cpus-per-task=1,time=00:30:00,mem=5G +GYS1 RT-qPCR

