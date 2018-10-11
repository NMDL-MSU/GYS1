---
title: GYS1 RT-qPCR
author: Deborah Velez-Irizarry
date: Wed Oct 10 17:09:50 EDT 2018
output:
  prettydoc::html_pretty:
    theme: tactile
    highlight: github
    toc: true
---
  
***  
### Description:  
Determine if glycogen synthase 1 gene (GYS1) is differentially expressed 
between NN, PN and PP genotypes.  
  
***  



```
##    Sample         ACTB          B2M        GAPDH         GYS1     Genotype
## 1    6000        24.32       18.709       15.553       19.993           NN
## 2    6001       22.005       17.155       14.811       20.643           NN
## 3    6008       22.873       18.321       15.668       19.974           PN
## 4    6055       22.593       18.004       14.746       19.648           PN
## 5      68       22.716       19.143       15.376       19.996           PP
## 6      53       22.904       18.007       15.744       20.417           PP
## 7    6003       22.504       18.178       15.517       19.542           NN
## 8    6004       22.199       18.314        15.21       20.196           NN
## 9    6007       22.669       18.559       15.442       20.116           NN
## 10   6062       21.802       18.219       15.419       19.609           PN
## 11    962        23.67        18.68       15.653       20.284           PP
## 12    NTC Undetermined Undetermined Undetermined Undetermined Undetermined
## 13   6000       23.401       18.786       15.348       20.014           NN
## 14   6001       21.995       17.134       14.846       20.584           NN
## 15   6008       22.237       18.326       16.064       19.917           PN
## 16   6055       22.592       17.964       14.691       19.669           PN
## 17     68       22.455       19.106       15.398        19.99           PP
## 18     53       22.874       17.975        15.75       20.373           PP
## 19   6003        22.23       18.101        15.53       19.518           NN
## 20   6004       22.133       18.301       15.218       20.183           NN
## 21   6007       21.977       18.566       15.515       20.092           NN
## 22   6062        21.43       18.251       15.398       19.652           PN
## 23    962       22.482       18.657       15.707       20.312           PP
## 24    NTC Undetermined           37 Undetermined Undetermined Undetermined
```




### Calculate deltaCT for target gene  
Results negative controls  


```
##    Sample         ACTB          B2M        GAPDH         GYS1     Genotype
## 12    NTC Undetermined Undetermined Undetermined Undetermined Undetermined
## 24    NTC Undetermined           37 Undetermined Undetermined Undetermined
```


Genotypes  


```
## 
## NN PN PP 
##  5  3  3
```

Control genes (three)  


```
## [1] "ACTB"  "B2M"   "GAPDH"
```


Average DeltaCT GYS1  


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


```
##         AVG     StDev
## NN 1.650302 0.7665144
## PN 1.374131 0.1392931
## PP 1.445220 0.2809440
```

Simple linear regression  


```
## Analysis of Variance Table
## 
## Response: dCT
##           Df  Sum Sq Mean Sq F value Pr(>F)
## Genotype   2 0.16549 0.08275  0.2599 0.7774
## Residuals  8 2.54684 0.31836
```

Calculata delta-delta CT and fold change comparing normal (NN) with affected (PN and PP)  


```
##    DeltaDeltaCT FoldChange
## PN   -0.2761705   1.210976
## PP   -0.2050814   1.152751
```


