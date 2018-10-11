#' ### Description:  
#' Determine if glycogen synthase 1 gene (GYS1) is differentially expressed 
#' between NN, PN and PP genotypes.  
#'   
#' ***  
#' **Code:**  
#' Directory:  
#'   
#' > &nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/qPCR/GYS1  
#'   
#' R Script:  
#'   
#' > &nbsp;&nbsp;&nbsp;&nbsp;GYS1.R  
#'  
#' **Input files:**  
#'   
#' > &nbsp;&nbsp;&nbsp;&nbsp;GYS1_qpcr.txt  
#'    
#' **Output files:**  
#'   
#' > &nbsp;&nbsp;&nbsp;&nbsp;GYS1.html    
#'   
#' ***  

#' ### Code  
#' Clear Environment  
rm(list=ls())

#' **Load required packages**
library(prettydoc)

#' **Session Information**  
sessionInfo()

#' ### Read in raw CT for GYS1  
system("mv ../GYS1_qpcr.txt .")
dir="/mnt/research/NMDL/qPCR/GYS1"
qpcr <- read.table(paste(dir, "GYS1_qpcr.txt", sep="/"), 
    header=TRUE)



#' ### Required Functions  
#' **Function to calculate the geometric mean of an array**  
#'   
#' > INPUT  
#' > &nbsp;&nbsp;&nbsp;&nbsp;`x` - vector of numeric values (NAs are removed)  
#'   
#' > OUTPUT  
#' > &nbsp;&nbsp;&nbsp;&nbsp;vector with geometric mean per sample  

geoMean <- function(x){
        prod(x[!is.na(x)])^(1/length(x[!is.na(x)]))
}

#' **Function to calculate average of raw intensity for each sample**  
#'   
#' > INPUT  
#' > &nbsp;&nbsp;&nbsp;&nbsp;`gene` - column name of target gene CT values, character scaler if repeated  
#' > &nbsp;&nbsp;&nbsp;&nbsp;measures are in a single column or character vector if repeated measures  
#' > &nbsp;&nbsp;&nbsp;&nbsp;are in multiple columns.  
#' > &nbsp;&nbsp;&nbsp;&nbsp;`nm` - character of gene name. Only uesed when byRow is TRUE.  
#' > &nbsp;&nbsp;&nbsp;&nbsp;The column with gene names must be called **Target** and the "nm" must match  
#' > &nbsp;&nbsp;&nbsp;&nbsp;a gene within the **Target** column.  
#' > &nbsp;&nbsp;&nbsp;&nbsp;`dataT` - data matrix with column of animal ID (must be named **Sample**)  
#' > &nbsp;&nbsp;&nbsp;&nbsp;and CT values for gene.   
#' > &nbsp;&nbsp;&nbsp;&nbsp;`byRow` - logical indicating if the repeated measures are in a single  
#' > &nbsp;&nbsp;&nbsp;&nbsp;column (multiple rows per sample), or in different columns (single row per  
#' > &nbsp;&nbsp;&nbsp;&nbsp;sample). Default=FALSE  
#'   
#' > OUTPUT  
#' > &nbsp;&nbsp;&nbsp;&nbsp;vector with mean CT value for repeated measures of a gene per sample  

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


#' **Function to calculate deltaCT values**   
#' Target gene will been the arithmetric mean of technical replicates and 
#' reference (control) genes will be the geometric mean.  
#'   
#' > INPUT  
#' > &nbsp;&nbsp;&nbsp;&nbsp;`gene` - column name of target gene CT values, character  
#' > &nbsp;&nbsp;&nbsp;&nbsp;scaler if repeated measures are in a single column or character  
#' > &nbsp;&nbsp;&nbsp;&nbsp;vector if repeated measures are in multiple columns.  
#' > &nbsp;&nbsp;&nbsp;&nbsp;`nm` - character of gene name. Only uesed when byRow is FALSE.  
#' > &nbsp;&nbsp;&nbsp;&nbsp;The column with gene names must be called "Target" and the "nm"  
#' > &nbsp;&nbsp;&nbsp;&nbsp;must match a gene within the "Target" column.  
#' > &nbsp;&nbsp;&nbsp;&nbsp;`control` - column names of control gene CT values, character  
#' > &nbsp;&nbsp;&nbsp;&nbsp;scaler if repeated measures are in a single column or character vector  
#' > &nbsp;&nbsp;&nbsp;&nbsp;if repeated measures are in multiple columns.  
#' > &nbsp;&nbsp;&nbsp;&nbsp;`dataT` - data matrix with column of animal ID (must be named "Sample")  
#' > &nbsp;&nbsp;&nbsp;&nbsp;and CT values for gene.  
#' > &nbsp;&nbsp;&nbsp;&nbsp;`byRow` - logical indicating if the repeated measures are in a single  
#' > &nbsp;&nbsp;&nbsp;&nbsp;column (multiple rows per sample), or in different columns (single  
#' > &nbsp;&nbsp;&nbsp;&nbsp;row per sample). Default=FALSE  
#'   
#' > OUTPUT  
#' > &nbsp;&nbsp;&nbsp;&nbsp;Data frame with:  
#' > &nbsp;&nbsp;&nbsp;&nbsp;`Target` - mean CT per sample for target gene per sample  
#' > &nbsp;&nbsp;&nbsp;&nbsp;`Control` - geometric mean of all control genes per sample  
#' > &nbsp;&nbsp;&nbsp;&nbsp;`dCT` - delta CT value for target gene normalized by control gene  
#' > &nbsp;&nbsp;&nbsp;&nbsp;expression per sample   

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
idx <- grep("NTC", qpcr$Sample)
qpcr[idx,]

#' Remove negative controls  
qpcr <- qpcr[-idx,]

#' Genotypes  
geno <- unique(qpcr[,c("Sample","Genotype")])
geno$Genotype <- factor(geno$Genotype)
qpcr <- qpcr[,-6]
table(geno$Genotype)

#' Control genes (three)  
control <- colnames(qpcr)[-c(1,5)]
control

#' Check qpcr data  
#' Read in as factor  
str(qpcr)

#' Change to numeric  
qpcr <- data.frame(apply(qpcr, 2, as.numeric))
str(qpcr)

#' DeltaCT GYS1  
dCT <- deltaCT(gene="GYS1", nm="GYS1", control=control, dataT=qpcr, byRow=TRUE)
dCT <- data.frame(dCT, Genotype=geno$Genotype)
dCT


#' ### Test Significance  
#' Average dCT per Genotype  
avg.geno <- do.call(rbind, lapply(as.character(unique(dCT$Genotype)), function(x) 
 data.frame(AVG=mean(dCT[dCT$Genotype == x, "dCT"]), StDev=sd(dCT[dCT$Genotype == x, "dCT"]))))
rownames(avg.geno) <- as.character(unique(dCT$Genotype))
avg.geno

#' Simple linear regression  
anova(lm(dCT ~ Genotype, data=dCT))

#' Calculata delta-delta CT and fold change comparing normal (NN) with affected (PN and PP)  
fold <- do.call(rbind, lapply(c("PN","PP"), function(x) data.frame(
    DeltaDeltaCT=avg.geno[rownames(avg.geno) == x, "AVG"] - 
        avg.geno[rownames(avg.geno) == "NN", "AVG"],
    FoldChange=2^-(avg.geno[rownames(avg.geno) == x, "AVG"] - 
        avg.geno[rownames(avg.geno) == "NN", "AVG"]))))
rownames(fold) <- c("PN","PP")
fold



#' ### Run R Script  
#+ eval = FALSE  
htmlRunR
GYS1.R  nodes=1,cpus-per-task=1,time=00:30:00,mem=5G +GYS1 RT-qPCR

