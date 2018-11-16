#!/bin/bash
#$ -cwd
#$ -N bolt_lmm_v7
#$ -o Logs
#$ -j y
#$ -pe shmem 7
#$ -V

# Downloaded from https://data.broadinstitute.org/alkesgroup/BOLT-LMM/downloads/ on May 12th 2016
bolt=/well/donnelly/ukbiobank_project_8874/clare/src/BOLT-LMM_v2.2/bolt 
basedir=/well/ukbiobank/qcoutput.V2_QCed.sample-QC

vers=v7
pheno=Intra.ocular.pressure..Goldmann.correlated..mean

### this sub-run is:
# BOLT version: LMM default - let BOLT decide on effect sizes.

mkdir $basedir/data/GWAS/otherGWAS/$pheno
mkdir $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers

genoDataAuto=/well/ukbiobank/expt/V2_QCed.export/data/genotype_gwas/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.genotype_gwas

genoDataSex=/well/ukbiobank/expt/V2_QCed.export/data/genotype_gwas/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.genotype_gwas

phenoFile=$basedir/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt  # using White british subset. Note: some individuals will be excluded in the genotype data

# phenoFile and snpSet was created using subsetPhenotypes.R
snpSet=$basedir/data/PCA/b1__b11-b001__b095-autosome-sampleqc-fastpca-White_British-highquality.prune.in  # just use SNPs from original PCA :) NOT ACTUALLY USED IN LREG!!  Note also that this does not include any SNPs on the sex chromosomes.

outFile=$basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out

# MEMORY NEEDED: 74,100,452,428 bytes = 5 C-NODES  <==== Sample QC SNPs
# MEMORY NEEDED: 96,739,226,940  = 7 C-NODES   <==== Oxford QC SNPs
# ALL CHROMOSOMES AT ONCE!!

# MUST MAKE SURE THAT THE .fam files for the autosomes and sex chroms are the same (they are as of May 16th)!!

# NOTE: This version uses 20 PCs as covariates - for LMM version as well!

$bolt \
--bfile=$genoDataAuto \
--lmmInfOnly \
--modelSnps=$snpSet \
--phenoFile=$phenoFile \
--phenoCol=$pheno \
--covarFile=$phenoFile \
--covarCol=Array \
--qCovarCol=Age.when.attended.assessment.centre \
--covarCol=Inferred.Gender \
--qCovarCol=PC{1:20} \
--numThreads=7 \
--statsFile=$outFile \
--covarMaxLevels=111 \
--geneticMapFile=/well/donnelly/ukbiobank_project_8874/clare/src/BOLT-LMM_v2.2/tables/genetic_map_hg19.txt.gz \
--verboseStats \
--maxMissingPerSnp=1 \
--maxMissingPerIndiv=1

# qsub -q short.qc -P donnelly.prjc runBOLT_v5-LMM.sh

