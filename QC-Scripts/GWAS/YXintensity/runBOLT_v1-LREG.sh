#!/bin/bash
#$ -cwd
#$ -N bolt_lreg-quant-Age-Ychrom_v1
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

genoData=$basedir/data/Combined/b1__b11-b001__b095-autosome-oxfordqc  #### JUST USING THE SAMPLE QC DATA (600,000 SNPs) Want to re-do this with the full QC'd 800,000
phenoFile=$basedir/data/GWAS/YXintensity/YchromPhenotypesForBOLT-$vers.txt

# phenoFile and snpSet was created using subsetPhenotypes.R
snpSet=$basedir/data/b1__b11-b001__b095-autosome-sampleqc-fastpca-init-highquality.prune.in  # just use SNPs from original PCA :) NOT ACTUALLY USED IN LREG!!

outFile=$basedir/data/GWAS/YXintensity/BOLTLMM.$vers/Ychrom-BOLT-LREG-quant-Age-$vers.out

# this run just outputs linear regression stats for all snps. i.e no lmm. Use PCs as covariates here.

# MEMORY NEEDED: 74,100,452,428 bytes = 5 C-NODES  <==== Sample QC SNPs
# MEMORY NEEDED: 96,739,226,940  = 7 C-NODES   <==== Oxford QC SNPs

$bolt \
--bfile=$genoData \
--phenoFile=$phenoFile \
--phenoCol=$pheno \
--covarFile=$phenoFile \
--covarCol=Array \
--qCovarCol=Age.when.attended.assessment.centre \
--qCovarCol=PC{1:8} \
--numThreads=7 \
--statsFile=$outFile \
--covarMaxLevels=111 \
--geneticMapFile=/well/donnelly/ukbiobank_project_8874/clare/src/BOLT-LMM_v2.2/tables/genetic_map_hg19.txt.gz \
--verboseStats

# qsub -q short.qc -P donnelly.prjc runBOLT_v1-LREG.sh

