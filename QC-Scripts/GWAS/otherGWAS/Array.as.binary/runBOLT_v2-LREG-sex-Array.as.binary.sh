#!/bin/bash
#$ -cwd
#$ -N bolt_lreg-sex_v2
#$ -o Logs
#$ -j y
#$ -pe shmem 5
#$ -V

# Downloaded from https://data.broadinstitute.org/alkesgroup/BOLT-LMM/downloads/ on May 12th 2016
bolt=/well/donnelly/ukbiobank_project_8874/clare/src/BOLT-LMM_v2.2/bolt 
basedir=/well/ukbiobank/qcoutput.V2_QCed.sample-QC

vers=v2
pheno=Array.as.binary

### this sub-run is:
# BOLT version: LREG - sex chromosomes 

mkdir $basedir/data/GWAS/otherGWAS/$pheno
mkdir $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers

genoData=$basedir/data/Combined/byChrom/b1__b11-b001__b095-sexchrom-oxfordqc-recoded  # Sex chromosomes recoded to 1,2,3 and 4

phenoFile=$basedir/data/GWAS/otherGWAS/PhenotypesForBOLT-$vers.txt

# phenoFile and snpSet was created using subsetPhenotypes.R
snpSet=$basedir/data/PCA/b1__b11-b001__b095-autosome-sampleqc-fastpca-init-highquality.prune.in  # just use SNPs from original PCA :) NOT ACTUALLY USED IN LREG!!  Note also that this does not include any SNPs on the sex chromosomes.

outFile=$basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LREG-sex-$vers.out

# MEMORY NEEDED: 74,100,452,428 bytes = 5 C-NODES  <==== Sample QC SNPs
# MEMORY NEEDED: 96,739,226,940  = 7 C-NODES   <==== Oxford QC SNPs
# ALL CHROMOSOMES AT ONCE!!

# SEX CHROMS: 2,326,447,666 bytes needed = 1 node!

# MUST MAKE SURE THAT THE .fam files for the autosomes and sex chroms are the same (they are as of May 16th)!!

$bolt \
--fam=$genoData-1.fam \
--bed=$genoData-{1:4}.bed \
--bim=$genoData-{1:4}.bim \
--phenoFile=$phenoFile \
--phenoCol=$pheno \
--covarFile=$phenoFile \
--qCovarCol=Age.when.attended.assessment.centre \
--qCovarCol=PC{1:20} \
--covarCol=Inferred.Gender \
--numThreads=5 \
--statsFile=$outFile \
--covarMaxLevels=111 \
--verboseStats


# qsub -q short.qa -P donnelly.prja runBOLT_v2-LREG-sex.sh

