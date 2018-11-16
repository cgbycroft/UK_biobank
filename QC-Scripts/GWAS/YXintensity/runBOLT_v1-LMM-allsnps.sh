#!/bin/bash
#$ -cwd
#$ -N bolt_lmm-all-snps-quant-Age-Ychrom_v1
#$ -o Logs
#$ -j y
#$ -pe shmem 8
#$ -V

# Downloaded from https://data.broadinstitute.org/alkesgroup/BOLT-LMM/downloads/ on May 12th 2016
bolt=/well/donnelly/ukbiobank_project_8874/clare/src/BOLT-LMM_v2.2/bolt 
basedir=/well/ukbiobank/qcoutput.V2_QCed.sample-QC

vers=v1
pheno=Y.X.intensity.ratio
### this sub-run is:
# BOLT version: LMM default - let BOLT decide on effect sizes.

mkdir $basedir/data/GWAS/YXintensity/BOLTLMM.$vers

genoData=$basedir/data/Combined/b1__b11-b001__b095-autosome-oxfordqc

phenoFile=$basedir/data/GWAS/YXintensity/YchromPhenotypesForBOLT-$vers.txt

# phenoFile was created using subsetPhenotypes.R
# snpSet was created using selectSNPs.R for Sample QC
snpSet=$basedir/QC-Scripts/Sample-QC/SelectSNPs/snpIncludeList-autosome.txt

outFile=$basedir/data/GWAS/YXintensity/BOLTLMM.$vers/Ychrom-BOLT-LMM-all-snps-quant-Age-$vers.out

# MEMORY NEEDED: 74,100,452,428 bytes = 5 C-NODES  <==== Sample QC SNPs
# MEMORY NEEDED: 96,739,226,940  = 7 C-NODES   <==== Oxford QC SNPs


$bolt \
--bfile $genoData \
--lmmInfOnly \
--phenoFile=$phenoFile \
--phenoCol=$pheno \
--covarFile=$phenoFile \
--covarCol=Array \
--qCovarCol=Age.when.attended.assessment.centre \
--numThreads=8 \
--statsFile=$outFile \
--covarMaxLevels=111 \
--geneticMapFile=/well/donnelly/ukbiobank_project_8874/clare/src/BOLT-LMM_v2.2/tables/genetic_map_hg19.txt.gz \
--verboseStats


# qsub -q coolibah.q -P donnelly.prjc runBOLT_v1-LMM-allsnps.sh

