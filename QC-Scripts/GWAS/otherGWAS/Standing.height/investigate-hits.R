source('/well/donnelly/ukbiobank_project_8874/clare/commonScripts/myManhattan.R')
source('/well/donnelly/ukbiobank_project_8874/clare/commonScripts/usefulFunctions.R')
h = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")
plink="/well/ukbiobank/qcoutput/Src/plink-v1.90b3k-x86_64/plink"
for(s in h) source(s)

library(dplyr)
library(qqman)
library(stringr)
library(binom)
library(sp)

vers ="v3"
phenoe = "Standing.height"

dataFile = paste0("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/",pheno,"/BOLTLMM.",vers,"/",pheno,"-BOLT-LMM-",vers,".out")

# read phenodata
phenoData = paste0(baseSampleQCDir,"/data/GWAS/otherGWAS/PhenotypesForBOLT-",vers,".txt")
pheno = read.table(phenoData,header=TRUE,stringsAsFactors=FALSE)


plotOutDir = "plots"
extraTitle=""

snpList = read.SNPQC.files()
QCexclusions = unique(c(snpList$arraySNPs,snpList$concordanceSNPs))

########
# READ IN GWAS RESULTS
