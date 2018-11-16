                                        # This is for the paper!
h = c("/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.snpqc-tests.R","/well/ukbiobank/expt/V2_QCed.SNP-QC/src/V2_QCed.bin2clusterplots.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/readPSperformance.R","/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/auxFunctions.R")
for(s in h) source(s)
source("/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/R/scripts/ukbbcolors.R")

library(igraph)
library(xtable)
library(hexbin)
library(latticeExtra)
library(viridis)


setwd('/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/forPaper')
###########
releaseSampleQCPrefix = 'b1__b11-b001__b095-sampleTable_v4'
###########


load(paste0('../../data/ForRelease/',releaseSampleQCPrefix,'_allColumns.RData'),verbose=TRUE)

cob =  get.place.of.birth(outTable)
eth2 = ethnicity2pop(outTable$Ethnic.background)
pop = ethnicity2pop2(outTable$Ethnic.background)

