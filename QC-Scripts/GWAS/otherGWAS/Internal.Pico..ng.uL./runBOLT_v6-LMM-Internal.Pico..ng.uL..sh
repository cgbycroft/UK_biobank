#!/bin/bash
#$ -cwd
#$ -N bolt_lmm_v6
#$ -o Logs
#$ -j y
#$ -pe shmem 4
#$ -V

# Downloaded from https://data.broadinstitute.org/alkesgroup/BOLT-LMM/downloads/ on May 12th 2016
bolt=/well/donnelly/ukbiobank_project_8874/clare/src/BOLT-LMM_v2.2/bolt 
basedir=/well/ukbiobank/qcoutput.V2_QCed.sample-QC

vers=v6
pheno=Internal.Pico..ng.uL.

### this sub-run is:
# BOLT version: LMM default - let BOLT decide on effect sizes.

mkdir $basedir/data/GWAS/otherGWAS/$pheno
mkdir $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers

#genoDataAuto=/well/ukbiobank/expt/V2_QCed.export/data/genotype_gwas/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.genotype_gwas

#genoDataSex=/well/ukbiobank/expt/V2_QCed.export/data/genotype_gwas/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.genotype_gwas

chr=$SGE_TASK_ID

sampleFile=/well/ukbiobank/phasing/final/phased_chunks/chr$chr.test2.sample
genoBgen=/well/ukbiobank/imputation/final/chr$chr/out/chr$chr.test2.hrc.I4.bgen 
genoDataPlink=$basedir/data/GWAS/otherGWAS/ImputedBgenSampleSet.$vers.dummy

phenoFile=$basedir/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt  # using White british subset. Note: some individuals will be excluded in the imputed data

# phenoFile was created using subsetPhenotypes.R

outFile=$basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-${vers}-chr$chr.out

# MEMORY NEEDED: 74,100,452,428 bytes = 5 C-NODES  <==== Sample QC SNPs
# MEMORY NEEDED: 96,739,226,940  = 7 C-NODES   <==== Oxford QC SNPs
# ALL CHROMOSOMES AT ONCE!!

# MUST MAKE SURE THAT THE .fam files for the autosomes and sex chroms are the same (they are as of May 16th)!!

# NOTE: This version only computes LREG because it's imputed data.


$bolt \
--bfile=$genoDataPlink \
--bgenFile=$genoBgen \
--sampleFile=$sampleFile \
--phenoFile=$phenoFile \
--phenoCol=$pheno \
--covarFile=$phenoFile \
--covarCol=Array \
--qCovarCol=Age.when.attended.assessment.centre \
--covarCol=Inferred.Gender \
--qCovarCol=PC{1:20} \
--numThreads=4 \
--statsFileBgenSnps=$outFile \
--statsFile=$outFile.dummy \
--covarMaxLevels=111 \
--verboseStats


# qsub -q short.qc -P donnelly.prjc runBOLT_v6-LMM.sh

