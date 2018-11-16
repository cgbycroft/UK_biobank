#!/bin/bash
#$ -cwd
#$ -N bolt_lmm-quant-Age-Ychrom_v1
#$ -o Logs
#$ -j y
#$ -pe shmem 10
#$ -V

# Downloaded from https://data.broadinstitute.org/alkesgroup/BOLT-LMM/downloads/ on May 12th 2016
bolt=/well/donnelly/ukbiobank_project_8874/clare/src/BOLT-LMM_v2.2/bolt 
basedir=/well/ukbiobank/qcoutput.V2_QCed.sample-QC

vers=v1
pheno=Y.X.intensity.ratio
### this sub-run is:
# BOLT version: LMM default - let BOLT decide on effect sizes.

mkdir $basedir/data/GWAS/YXintensity/BOLTLMM.$vers

genoDataAuto=$basedir/data/Combined/b1__b11-b001__b095-autosome-oxfordqc
genoDataSex=$basedir/data/Combined/b1__b11-b001__b095-sexchrom-oxfordqc

phenoFile=$basedir/data/GWAS/YXintensity/YchromPhenotypesForBOLT-$vers.txt

# phenoFile and snpSet was created using subsetPhenotypes.R
snpSet=$basedir/data/PCA/b1__b11-b001__b095-autosome-sampleqc-fastpca-init-highquality.prune.in  # just use SNPs from original PCA :) NOT ACTUALLY USED IN LREG!!

outFile=$basedir/data/GWAS/YXintensity/BOLTLMM.$vers/Ychrom-BOLT-LMM-quant-Age-$vers.out

# MEMORY NEEDED: 74,100,452,428 bytes = 5 C-NODES  <==== Sample QC SNPs
# MEMORY NEEDED: 96,739,226,940  = 7 C-NODES   <==== Oxford QC SNPs

$bolt \
--fam $genoDataAuto.fam \
--bed $genoDataAuto.bed \
--bim $genoDataAuto.bim \
--bed $genoDataSex.bed \
--bim $genoDataSex.bim \
--lmmInfOnly \
--modelSnps=$snpSet \
--phenoFile=$phenoFile \
--phenoCol=$pheno \
--covarFile=$phenoFile \
--covarCol=Array \
--qCovarCol=Age.when.attended.assessment.centre \
--numThreads=10 \
--statsFile=$outFile \
--covarMaxLevels=111 \
--geneticMapFile=/well/donnelly/ukbiobank_project_8874/clare/src/BOLT-LMM_v2.2/tables/genetic_map_hg19.txt.gz \
--verboseStats


# qsub -q short.qc -P donnelly.prjc runBOLT_v1-LMM.sh

