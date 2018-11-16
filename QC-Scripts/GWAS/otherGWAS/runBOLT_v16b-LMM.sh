#!/bin/bash
#$ -cwd
#$ -N bolt_lmm_versionNumber
#$ -o Logs
#$ -j y
#$ -pe shmem 7
#$ -V

# Downloaded from https://data.broadinstitute.org/alkesgroup/BOLT-LMM/downloads/ on May 12th 2016
bolt=/well/donnelly/ukbiobank_project_8874/clare/src/BOLT-LMM_v2.2/bolt 
basedir=/well/ukbiobank/qcoutput.V2_QCed.sample-QC

vers=versionNumber
pheno=phenotypeForThisFile

### this sub-run is:
# BOLT version: LMM default - let BOLT decide on effect sizes.

mkdir $basedir/data/GWAS/otherGWAS/$pheno
mkdir $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers

#genoDataAuto=/well/ukbiobank/expt/V2_QCed.export/data/genotype_gwas/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.genotype_gwas
genoDataAuto=/well/ukbiobank/expt/V2_QCed.export/data/calls_export/v3/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.chr{1:22}  # <=== Release data as of March 6 2017. Separated by chromosome.

phenoFile=$basedir/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt  # using White british subset. Note: some individuals will be excluded in the genotype data

# phenoFile and snpSet was created using subsetPhenotypes.R
#snpSet=$basedir/data/PCA/b1__b11-b001__b095-autosome-sampleqc-fastpca-White_British-highquality.prune.in  # just use SNPs from original PCA :) NOT ACTUALLY USED IN LREG!!  Note also that this does not include any SNPs on the sex chromosomes.
snpSet=$basedir/data/GWAS/otherGWAS/b1__b11-b001__b095-autosome-sampleqc-fastpca-White_British-highquality-pruned.snpsetRSIDs.txt   # <=== snps in rsid format to match the release genotype data.


# file for LD score regression (matched on bp)
LDscoresFile=/well/donnelly/ukbiobank_project_8874/clare/src/BOLT-LMM_v2.2/tables/LDSCORE.1000G_EUR.tab.gz

###############
# THIS IS THE ONLY DIFFERENCE FROM v16.  Sampling 10,000 random individuals.  The same set as for v19b

# List of sample exclusions so that the input matches the imputed data run (v15)
samplesToRemove=$basedir/QC-Scripts/GWAS/otherGWAS/samplesToExclude.$vers.txt
###############


outFile=$basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out

# MEMORY NEEDED: 74,100,452,428 bytes = 5 C-NODES  <==== Sample QC SNPs
# MEMORY NEEDED: 96,739,226,940  = 7 C-NODES   <==== Oxford QC SNPs
# ALL CHROMOSOMES AT ONCE!!

# MUST MAKE SURE THAT THE .fam files for the autosomes and sex chroms are the same (they are as of May 16th)!!

# NOTE: This version uses 20 PCs as covariates - for LMM version as well!
# To adjust default snp missing rate threshold (0.1) use --maxMissingPerSnp. e.g --maxMissingPerSnp=1

# NOTE: Use of one *fam file assumes that all fam files are int he same order. To show this run:
# for i in `seq 1 21`; do diff /well/ukbiobank/expt/V2_QCed.export/data/calls_export/v3/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.chr22.fam /well/ukbiobank/expt/V2_QCed.export/data/calls_export/v3/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.chr$i.fam; done;


$bolt \
--fam=/well/ukbiobank/expt/V2_QCed.export/data/calls_export/v3/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.chr22.fam \
--bed=$genoDataAuto.bed \
--bim=$genoDataAuto.bim \
--lmmInfOnly \
--remove $samplesToRemove \
--modelSnps=$snpSet \
--phenoFile=$phenoFile \
--phenoCol=$pheno \
--covarFile=$phenoFile \
--covarCol=Array \
--qCovarCol=Age.when.attended.assessment.centre \
--covarCol=Inferred.Gender \
--qCovarCol=PC{1:20} \
--maxMissingPerSnp 1 \
--numThreads=7 \
--statsFile=$outFile \
--covarMaxLevels=111 \
--geneticMapFile=/well/donnelly/ukbiobank_project_8874/clare/src/BOLT-LMM_v2.2/tables/genetic_map_hg19.txt.gz \
--LDscoresFile=$LDscoresFile \
--LDscoresMatchBp \
--verboseStats


# qsub -q short.qc -P donnelly.prjc runBOLT_versionNumber-LMM.sh

