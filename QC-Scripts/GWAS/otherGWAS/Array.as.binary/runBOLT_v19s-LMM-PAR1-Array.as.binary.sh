#!/bin/bash
#$ -cwd
#$ -N bolt_lmm_v19s_PAR1
#$ -o Logs
#$ -j y
#$ -pe shmem 2
#$ -V
#$ -t 1-6

# Downloaded from https://data.broadinstitute.org/alkesgroup/BOLT-LMM/downloads/ on May 12th 2016
bolt=/well/donnelly/ukbiobank_project_8874/clare/src/BOLT-LMM_v2.2/bolt 
basedir=/well/ukbiobank/qcoutput.V2_QCed.sample-QC

# RUNS CHROM X IN CHUNKS!
vers=v19s
pheno=Array.as.binary

### this sub-run is:
# BOLT version: LMM default - let BOLT decide on effect sizes.

mkdir $basedir/data/GWAS/otherGWAS/$pheno
mkdir $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers

#genoDataAuto=/well/ukbiobank/expt/V2_QCed.export/data/genotype_gwas/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.genotype_gwas

#genoDataSex=/well/ukbiobank/expt/V2_QCed.export/data/genotype_gwas/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.genotype_gwas

i=$SGE_TASK_ID  # These are the bgen chunks
#chr=2 # For when running directly on a server

chr=PAR1

sampleFile=/well/ukbiobank/imputation/final/full/bgen/hrc+uk10k_sorted_8bit_rsids/v2/chrPAR.sample
genoDataPlink=$basedir/data/GWAS/otherGWAS/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.ImputedBgenSampleSet-for-$vers.PAR.modelSnps  # Subset the samples in genotype data to match the imputed data set samples. I.e genotypes for release but subset for model SNPs and set of samples in .bgen data. This was done in prelims.sh (requirement of BOLT!)


# Chunked X chrom files created in subset-sex-chroms.sh
genoBgen=$basedir/data/imputed/chr$chr.v2.$i.recoded.bgen


phenoFile=$basedir/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt  # using White british subset. Note: some individuals will be excluded in the imputed data. CHECK bolt output for exact numbers.

# phenoFile was created using subsetPhenotypes.R

# file for LD score regression
LDscoresFile=/well/donnelly/ukbiobank_project_8874/clare/src/BOLT-LMM_v2.2/tables/LDSCORE.1000G_EUR.tab.gz

snpSet=$basedir/QC-Scripts/GWAS/otherGWAS/dummy.snps.txt # just use SNPs from original PCA :) NOT ACTUALLY USED IN LREG!!  Might want to try different set... e.g with lower MAF threshold.

outFile=$basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-${vers}-chr$chr.$i.out

# MEMORY NEEDED: 74,100,452,428 bytes = 5 C-NODES  <==== Sample QC SNPs
# MEMORY NEEDED: 96,739,226,940  = 7 C-NODES   <==== Oxford QC SNPs
# ALL CHROMOSOMES AT ONCE!!

# MUST MAKE SURE THAT THE .fam files for the autosomes and sex chroms are the same (they are as of May 16th)!!

# NOTE: This version computes both LREG and LMM versions. effect sizes only available for betas.
# ALSO NOTE: genotypes and imputed data should all be aligned to the reference allele!?
# Forced bolt's snp maxmissing filter to 1. Should not actually matter with .bgen imputed files...

# NOTE: LMM version doesn't make sense as the chromosome number has been changed!


$bolt \
--bfile=$genoDataPlink \
--bgenFile=$genoBgen \
--sampleFile=$sampleFile \
--phenoFile=$phenoFile \
--phenoCol=$pheno \
--covarFile=$phenoFile \
--qCovarCol=Age.when.attended.assessment.centre \
--covarCol=Inferred.Gender \
--qCovarCol=PC{1:20} \
--maxMissingPerSnp 1 \
--numThreads=2 \
--statsFileBgenSnps=$outFile \
--statsFile=$outFile.geno \
--covarMaxLevels=111 \
--geneticMapFile=/well/donnelly/ukbiobank_project_8874/clare/src/BOLT-LMM_v2.2/tables/genetic_map_hg19.txt.gz \
--LDscoresFile=$LDscoresFile \
--LDscoresMatchBp \
--verboseStats


# qsub -q short.qc -P donnelly.prjc runBOLT_v19s-LMM.sh

