#!/bin/bash
#$ -cwd
#$ -N bolt_lmm_v1-all-snps
#$ -o Logs
#$ -j y
#$ -pe shmem 13
#$ -V

# NOTE: if running on jarrah, need to ask for at least 13 cores, at 8GB per core.

# Downloaded from https://data.broadinstitute.org/alkesgroup/BOLT-LMM/downloads/ on May 12th 2016
bolt=/well/donnelly/ukbiobank_project_8874/clare/src/BOLT-LMM_v2.2/bolt 
basedir=/well/ukbiobank/qcoutput.V2_QCed.sample-QC

vers=v1
pheno=Place.of.birth.in.UK...north.co.ordinate

### this sub-run is:
# BOLT version: LMM default - let BOLT decide on effect sizes.

mkdir $basedir/data/GWAS/otherGWAS/$pheno
mkdir $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers

genoDataAuto=$basedir/data/Combined/b1__b11-b001__b095-autosome-oxfordqc
#genoDataSex=$basedir/data/Combined/b1__b11-b001__b095-sexchrom-oxfordqc

phenoFile=$basedir/data/GWAS/otherGWAS/PhenotypesForBOLT-$vers.txt

# phenoFile was created using subsetPhenotypes.R
# snpSet was created using selectSNPs.R for Sample QC

snpSet=$basedir/QC-Scripts/Sample-QC/SelectSNPs/snpIncludeList-autosome.txt

outFile=$basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-all-snps-$vers.out

# MEMORY NEEDED: 74,100,452,428 bytes = 5 C-NODES  <==== Sample QC SNPs
# MEMORY NEEDED: 96,739,226,940  = 7 C-NODES   <==== Oxford QC SNPs
# ALL CHROMOSOMES AT ONCE!!


$bolt \
--bfile $genoDataAuto \
--lmmInfOnly \
--modelSnps=$snpSet \
--phenoFile=$phenoFile \
--phenoCol=$pheno \
--covarFile=$phenoFile \
--covarCol=Array \
--qCovarCol=Age.when.attended.assessment.centre \
--covarCol=Inferred.Gender \
--numThreads=13 \
--statsFile=$outFile \
--covarMaxLevels=111 \
--geneticMapFile=/well/donnelly/ukbiobank_project_8874/clare/src/BOLT-LMM_v2.2/tables/genetic_map_hg19.txt.gz \
--verboseStats


# qsub -q short.qc -P donnelly.prjc runBOLT_v1-LMM.sh

