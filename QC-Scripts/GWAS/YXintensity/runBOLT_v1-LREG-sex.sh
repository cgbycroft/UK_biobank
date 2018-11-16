#!/bin/bash
#$ -cwd
#$ -N bolt_lreg-sex-quant-Age-Ychrom_v1
#$ -o Logs
#$ -j y
#$ -pe shmem 7
#$ -V

# Downloaded from https://data.broadinstitute.org/alkesgroup/BOLT-LMM/downloads/ on May 12th 2016
bolt=/well/donnelly/ukbiobank_project_8874/clare/src/BOLT-LMM_v2.2/bolt 
basedir=/well/ukbiobank/qcoutput.V2_QCed.sample-QC

vers=v1
pheno=Y.X.intensity.ratio
### this sub-run is:
# BOLT version: LREG 

mkdir $basedir/data/GWAS/YXintensity/BOLTLMM.$vers

genoData=$basedir/data/Combined/byChrom/b1__b11-b001__b095-sexchrom-oxfordqc-recoded  # Sex chromosomes recoded to 1,2,3 and 4

phenoFile=$basedir/data/GWAS/YXintensity/YchromPhenotypesForBOLT-$vers.txt

# phenoFile and snpSet was created using subsetPhenotypes.R
snpSet=$basedir/data/b1__b11-b001__b095-autosome-sampleqc-fastpca-init-highquality.prune.in  # just use SNPs from original PCA :) NOT ACTUALLY USED IN LREG!!

outFile=$basedir/data/GWAS/YXintensity/BOLTLMM.$vers/Ychrom-BOLT-LREG-quant-Age-$vers-sexchroms.out

# this run just outputs linear regression stats for all snps. i.e no lmm. Use PCs as covariates here.

# MEMORY NEEDED: 74,100,452,428 bytes = 5 C-NODES  <==== Sample QC SNPs
# MEMORY NEEDED: 96,739,226,940  = 7 C-NODES   <==== Oxford QC SNPs

$bolt \
--fam=$genoData-1.fam \
--bed=$genoData-{1:4}.bed \
--bim=$genoData-{1:4}.bim \
--phenoFile=$phenoFile \
--phenoCol=$pheno \
--covarFile=$phenoFile \
--covarCol=Array \
--qCovarCol=Age.when.attended.assessment.centre \
--qCovarCol=PC{1:8} \
--numThreads=7 \
--statsFile=$outFile \
--covarMaxLevels=111 \
--verboseStats

# qsub -q short.qa -P donnelly.prja runBOLT_v1-LREG-sex.sh

