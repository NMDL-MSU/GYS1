---
title: GYS1 RT-qPCR
author: Deborah Velez-Irizarry
date: Wed Oct 10 17:09:23 EDT 2018
output:
  prettydoc::html_pretty:
    theme: tactile
    highlight: github
    toc: true
---
### Description:  
Determine if glycogen synthase 1 gene (GYS1) is differentially expressed 
between NN, PN and PP genotypes.  
  
***  
**Code:**  
Directory:  
  
> &nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/qPCR/GYS1  
  
R Script:  
  
> &nbsp;&nbsp;&nbsp;&nbsp;GYS1.R  
 
**Input files:**  
  
> &nbsp;&nbsp;&nbsp;&nbsp;GYS1_qpcr.txt  
   
**Output files:**  
  
> &nbsp;&nbsp;&nbsp;&nbsp;GYS1.html    
  
***  
### Code  
Clear Environment  


```r
rm(list=ls())
```

**Load required packages**


```r
library(prettydoc)
```

**Session Information**  


```r
sessionInfo()
```

```
## R version 3.5.1 (2018-07-02)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: CentOS Linux 7 (Core)
## 
## Matrix products: default
## BLAS/LAPACK: /opt/software/OpenBLAS/0.2.20-GCC-6.4.0-2.28/lib/libopenblas_haswellp-r0.2.20.so
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] prettydoc_0.2.1 knitr_1.20     
## 
## loaded via a namespace (and not attached):
##  [1] compiler_3.5.1  backports_1.1.2 magrittr_1.5    rprojroot_1.3-2
##  [5] htmltools_0.3.6 tools_3.5.1     Rcpp_0.12.17    rmarkdown_1.10 
##  [9] stringi_1.2.3   digest_0.6.15   stringr_1.3.1   evaluate_0.10.1
```

### Read in raw CT for GYS1  


```r
system("mv ../GYS1_qpcr.txt .")
dir="/mnt/research/NMDL/qPCR/GYS1"
qpcr <- read.table(paste(dir, "GYS1_qpcr.txt", sep="/"), 
    header=TRUE)
```

### Required Functions  
**Function to calculate the geometric mean of an array**  
  
> INPUT  
> &nbsp;&nbsp;&nbsp;&nbsp;`x` - vector of numeric values (NAs are removed)  
  
> OUTPUT  
> &nbsp;&nbsp;&nbsp;&nbsp;vector with geometric mean per sample  


```r
geoMean <- function(x){
        prod(x[!is.na(x)])^(1/length(x[!is.na(x)]))
}
```

**Function to calculate average of raw intensity for each sample**  
  
> INPUT  
> &nbsp;&nbsp;&nbsp;&nbsp;`gene` - column name of target gene CT values, character scaler if repeated  
> &nbsp;&nbsp;&nbsp;&nbsp;measures are in a single column or character vector if repeated measures  
> &nbsp;&nbsp;&nbsp;&nbsp;are in multiple columns.  
> &nbsp;&nbsp;&nbsp;&nbsp;`nm` - character of gene name. Only uesed when byRow is TRUE.  
> &nbsp;&nbsp;&nbsp;&nbsp;The column with gene names must be called **Target** and the "nm" must match  
> &nbsp;&nbsp;&nbsp;&nbsp;a gene within the **Target** column.  
> &nbsp;&nbsp;&nbsp;&nbsp;`dataT` - data matrix with column of animal ID (must be named **Sample**)  
> &nbsp;&nbsp;&nbsp;&nbsp;and CT values for gene.   
> &nbsp;&nbsp;&nbsp;&nbsp;`byRow` - logical indicating if the repeated measures are in a single  
> &nbsp;&nbsp;&nbsp;&nbsp;column (multiple rows per sample), or in different columns (single row per  
> &nbsp;&nbsp;&nbsp;&nbsp;sample). Default=FALSE  
  
> OUTPUT  
> &nbsp;&nbsp;&nbsp;&nbsp;vector with mean CT value for repeated measures of a gene per sample  


```r
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
```

**Function to calculate deltaCT values**   
Target gene will been the arithmetric mean of technical replicates and 
reference (control) genes will be the geometric mean.  
  
> INPUT  
> &nbsp;&nbsp;&nbsp;&nbsp;`gene` - column name of target gene CT values, character  
> &nbsp;&nbsp;&nbsp;&nbsp;scaler if repeated measures are in a single column or character  
> &nbsp;&nbsp;&nbsp;&nbsp;vector if repeated measures are in multiple columns.  
> &nbsp;&nbsp;&nbsp;&nbsp;`nm` - character of gene name. Only uesed when byRow is FALSE.  
> &nbsp;&nbsp;&nbsp;&nbsp;The column with gene names must be called "Target" and the "nm"  
> &nbsp;&nbsp;&nbsp;&nbsp;must match a gene within the "Target" column.  
> &nbsp;&nbsp;&nbsp;&nbsp;`control` - column names of control gene CT values, character  
> &nbsp;&nbsp;&nbsp;&nbsp;scaler if repeated measures are in a single column or character vector  
> &nbsp;&nbsp;&nbsp;&nbsp;if repeated measures are in multiple columns.  
> &nbsp;&nbsp;&nbsp;&nbsp;`dataT` - data matrix with column of animal ID (must be named "Sample")  
> &nbsp;&nbsp;&nbsp;&nbsp;and CT values for gene.  
> &nbsp;&nbsp;&nbsp;&nbsp;`byRow` - logical indicating if the repeated measures are in a single  
> &nbsp;&nbsp;&nbsp;&nbsp;column (multiple rows per sample), or in different columns (single  
> &nbsp;&nbsp;&nbsp;&nbsp;row per sample). Default=FALSE  
  
> OUTPUT  
> &nbsp;&nbsp;&nbsp;&nbsp;Data frame with:  
> &nbsp;&nbsp;&nbsp;&nbsp;`Target` - mean CT per sample for target gene per sample  
> &nbsp;&nbsp;&nbsp;&nbsp;`Control` - geometric mean of all control genes per sample  
> &nbsp;&nbsp;&nbsp;&nbsp;`dCT` - delta CT value for target gene normalized by control gene  
> &nbsp;&nbsp;&nbsp;&nbsp;expression per sample   


```r
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
```

### Calculate deltaCT for target gene  
Results negative controls  


```r
idx <- grep("NTC", qpcr$Sample)
qpcr[idx,]
```

```
##    Sample         ACTB          B2M        GAPDH         GYS1     Genotype
## 12    NTC Undetermined Undetermined Undetermined Undetermined Undetermined
## 24    NTC Undetermined           37 Undetermined Undetermined Undetermined
```

Remove negative controls  


```r
qpcr <- qpcr[-idx,]
```

Genotypes  


```r
geno <- unique(qpcr[,c("Sample","Genotype")])
geno$Genotype <- factor(geno$Genotype)
qpcr <- qpcr[,-6]
table(geno$Genotype)
```

```
## 
## NN PN PP 
##  5  3  3
```

Control genes (three)  


```r
control <- colnames(qpcr)[-c(1,5)]
control
```

```
## [1] "ACTB"  "B2M"   "GAPDH"
```

Check qpcr data  
Read in as factor  


```r
str(qpcr)
```

```
## 'data.frame':	22 obs. of  5 variables:
##  $ Sample: Factor w/ 12 levels "53","6000","6001",..: 2 3 7 8 10 1 4 5 6 9 ...
##  $ ACTB  : Factor w/ 23 levels "21.43","21.802",..: 22 5 17 14 16 19 12 7 15 2 ...
##  $ B2M   : Factor w/ 24 levels "17.134","17.155",..: 19 2 13 5 22 6 8 12 15 9 ...
##  $ GAPDH : Factor w/ 22 levels "14.691","14.746",..: 15 3 17 2 8 19 13 5 11 10 ...
##  $ GYS1  : Factor w/ 23 levels "19.518","19.542",..: 10 22 8 4 11 20 2 16 14 3 ...
```

Change to numeric  


```r
qpcr <- data.frame(apply(qpcr, 2, as.numeric))
str(qpcr)
```

```
## 'data.frame':	22 obs. of  5 variables:
##  $ Sample: num  6000 6001 6008 6055 68 ...
##  $ ACTB  : num  24.3 22 22.9 22.6 22.7 ...
##  $ B2M   : num  18.7 17.2 18.3 18 19.1 ...
##  $ GAPDH : num  15.6 14.8 15.7 14.7 15.4 ...
##  $ GYS1  : num  20 20.6 20 19.6 20 ...
```

DeltaCT GYS1  


```r
dCT <- deltaCT(gene="GYS1", nm="GYS1", control=control, dataT=qpcr, byRow=TRUE)
dCT <- data.frame(dCT, Genotype=geno$Genotype)
dCT
```

```
##       Target  Control       dCT Genotype
## 6000 20.0035 19.04692 0.9565778       NN
## 6001 20.6135 17.75067 2.8628289       NN
## 6008 19.9455 18.71604 1.2294625       PN
## 6055 19.6585 18.15116 1.5073406       PN
## 68   19.9930 18.80134 1.1916632       PP
## 53   20.3950 18.64776 1.7472448       PP
## 6003 19.5300 18.46735 1.0626529       NN
## 6004 20.1895 18.34509 1.8444139       NN
## 6007 20.1040 18.57896 1.5250356       NN
## 6062 19.6305 18.24491 1.3855908       PN
## 962  20.2980 18.90125 1.3967532       PP
```

### Test Significance  
Average dCT per Genotype  


```r
avg.geno <- do.call(rbind, lapply(as.character(unique(dCT$Genotype)), function(x) 
 data.frame(AVG=mean(dCT[dCT$Genotype == x, "dCT"]), StDev=sd(dCT[dCT$Genotype == x, "dCT"]))))
rownames(avg.geno) <- as.character(unique(dCT$Genotype))
avg.geno
```

```
##         AVG     StDev
## NN 1.650302 0.7665144
## PN 1.374131 0.1392931
## PP 1.445220 0.2809440
```

Simple linear regression  


```r
anova(lm(dCT ~ Genotype, data=dCT))
```

```
## Analysis of Variance Table
## 
## Response: dCT
##           Df  Sum Sq Mean Sq F value Pr(>F)
## Genotype   2 0.16549 0.08275  0.2599 0.7774
## Residuals  8 2.54684 0.31836
```

Calculata delta-delta CT and fold change comparing normal (NN) with affected (PN and PP)  


```r
fold <- do.call(rbind, lapply(c("PN","PP"), function(x) data.frame(
    DeltaDeltaCT=avg.geno[rownames(avg.geno) == x, "AVG"] - 
        avg.geno[rownames(avg.geno) == "NN", "AVG"],
    FoldChange=2^-(avg.geno[rownames(avg.geno) == x, "AVG"] - 
        avg.geno[rownames(avg.geno) == "NN", "AVG"]))))
rownames(fold) <- c("PN","PP")
fold
```

```
##    DeltaDeltaCT FoldChange
## PN   -0.2761705   1.210976
## PP   -0.2050814   1.152751
```

### Run R Script  


```r
htmlRunR
GYS1.R  nodes=1,cpus-per-task=1,time=00:30:00,mem=5G +GYS1 RT-qPCR
```

