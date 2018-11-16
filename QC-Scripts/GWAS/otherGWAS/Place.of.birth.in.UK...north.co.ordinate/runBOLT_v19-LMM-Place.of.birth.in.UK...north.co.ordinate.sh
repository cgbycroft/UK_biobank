#!/bin/bash
#$ -cwd
#$ -N bolt_lmm_v19
#$ -o Logs
#$ -j y
#$ -pe shmem 5

#$ -V

# Downloaded from https://data.broadinstitute.org/alkesgroup/BOLT-LMM/downloads/ on May 12th 2016
bolt=/well/donnelly/ukbiobank_project_8874/clare/src/BOLT-LMM_v2.2/bolt 
basedir=/well/ukbiobank/qcoutput.V2_QCed.sample-QC

vers=v19
pheno=Place.of.birth.in.UK...north.co.ordinate

### this sub-run is:
# BOLT version: LMM default - let BOLT decide on effect sizes.

mkdir $basedir/data/GWAS/otherGWAS/$pheno
mkdir $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers

#genoDataAuto=/well/ukbiobank/expt/V2_QCed.export/data/genotype_gwas/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.genotype_gwas

#genoDataSex=/well/ukbiobank/expt/V2_QCed.export/data/genotype_gwas/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.genotype_gwas

chr=$SGE_TASK_ID
#chr=2 # For when running directly on a server


### chromosomes 16-22 were run using files in the below directory (November 2016). 
#sampleFile=/well/ukbiobank/imputation/final/full/bgen/chr$chr.hrc+uk10k.I4.v1.1.sample
#genoBgen=/well/ukbiobank/imputation/final/full/bgen/chr$chr.hrc+uk10k.I4.v1.1.bgen
## Then chromosomes 1-15 were run using the new directory set up by Jonathan (December 2016). Data should be the same, just shifted into the directory hrc+uk10k_sorted_v1.1
#sampleFile=/well/ukbiobank/imputation/final/full/bgen/hrc+uk10k_sorted_v1.1/chr$chr.hrc+uk10k.I4.v1.1.sample
#genoBgen=/well/ukbiobank/imputation/final/full/bgen/hrc+uk10k_sorted_v1.1/chr$chr.hrc+uk10k.I4.v1.1.bgen
## ===> Files created by Colin from bgen files v1.2
#sampleFile=/well/ukbiobank/imputation/final/full/bgen/hrc+uk10k_sorted_8bit_rsids/imputed.sample
#genoBgen=/well/ukbiobank/expt/V2_QCed.imputation.sanity_check/data/chr$chr.hrc+uk10k_sorted_8bit_rsids.v1.1.bgen
sampleFile=/well/ukbiobank/imputation/final/full/bgen/hrc+uk10k_sorted_8bit_rsids/imputed.sample
#genoBgen=/well/ukbiobank/expt/full_release_issues/uk10k_imputation_annotation/hrc+uk10k_sorted_8bit_rsids_chr$chr.v2.v1_1.bgen
genoBgen=/well/ukbiobank/expt/full_release_issues/uk10k_imputation_annotation/hrc+uk10k_sorted_8bit_rsids_chr$chr.v4.v1_1.bgen

genoDataPlink=$basedir/data/GWAS/otherGWAS/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.ImputedBgenSampleSet-for-v15.modelSnps  # Subset the samples in genotype data to match the imputed data set samples. I.e genotypes for release but subset for model SNPs and set of samples in .bgen data. This was done in prelims.sh (requirement of BOLT!)

phenoFile=$basedir/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt  # using White british subset. Note: some individuals will be excluded in the imputed data. CHECK bolt output for exact numbers.

# phenoFile was created using subsetPhenotypes.R

# file for LD score regression
LDscoresFile=/well/donnelly/ukbiobank_project_8874/clare/src/BOLT-LMM_v2.2/tables/LDSCORE.1000G_EUR.tab.gz

snpSet=$basedir/data/GWAS/otherGWAS/b1__b11-b001__b095-autosome-sampleqc-fastpca-White_British-highquality-pruned.snpsetRSIDs.txt  # just use SNPs from original PCA :) NOT ACTUALLY USED IN LREG!!  Might want to try different set... e.g with lower MAF threshold.

outFile=$basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-${vers}-chr$chr.out

# MEMORY NEEDED: 74,100,452,428 bytes = 5 C-NODES  <==== Sample QC SNPs
# MEMORY NEEDED: 96,739,226,940  = 7 C-NODES   <==== Oxford QC SNPs
# ALL CHROMOSOMES AT ONCE!!

# MUST MAKE SURE THAT THE .fam files for the autosomes and sex chroms are the same (they are as of May 16th)!!

# NOTE: This version computes both LREG and LMM versions. effect sizes only available for betas.
# ALSO NOTE: genotypes and imputed data should all be aligned to the reference allele!?
# Forced bolt's snp maxmissing filter to 1. Should not actually matter with .bgen imputed files...


$bolt \
--bfile=$genoDataPlink \
--bgenFile=$genoBgen \
--sampleFile=$sampleFile \
--lmmInfOnly \
--modelSnps=$snpSet \
--phenoFile=$phenoFile \
--phenoCol=$pheno \
--covarFile=$phenoFile \
--covarCol=Array \
--qCovarCol=Age.when.attended.assessment.centre \
--covarCol=Inferred.Gender \
--qCovarCol=PC{1:20} \
--maxMissingPerSnp 1 \
--numThreads=5 \
--statsFileBgenSnps=$outFile \
--statsFile=$outFile.geno \
--covarMaxLevels=111 \
--geneticMapFile=/well/donnelly/ukbiobank_project_8874/clare/src/BOLT-LMM_v2.2/tables/genetic_map_hg19.txt.gz \
--LDscoresFile=$LDscoresFile \
--LDscoresMatchBp \
--verboseStats


# qsub -q short.qc -P donnelly.prjc runBOLT_v19-LMM.sh

