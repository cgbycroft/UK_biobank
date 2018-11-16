#!/bin/bash
#$ -cwd
#$ -N bolt_lmm_versionNumber
#$ -o Logs
#$ -j y
#$ -pe shmem 4
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

#genoDataSex=/well/ukbiobank/expt/V2_QCed.export/data/genotype_gwas/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.genotype_gwas

chr=$SGE_TASK_ID

# Files from Jonathan
sampleFile=/well/ukbiobank/phasing/final/phased_chunks/chr$chr.test2.sample
genoBgen=/well/ukbiobank/imputation/final/chr$chr/out/chr$chr.test2.hrc.I4.bgen 

# Dummy file so bolt runs.
genoDataPlink=$basedir/data/GWAS/otherGWAS/ImputedBgenSampleSet.$vers.dummy

# phenoFile was created using subsetPhenotypes.R. Use v3 for WhiteBritish.
phenoFile=$basedir/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt  # using White british subset. Note: some individuals will be excluded in the imputed data

# Keep this naming convention as it is recognised by R scripts used to plot results.
outFile=$basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-${vers}-chr$chr.out


# NOTE: This version only computes LREG because it's imputed data.

# To adjust default snp missing rate threshold (0.1) use --maxMissingPerSnp. e.g --maxMissingPerSnp=0

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


# qsub -q short.qc -P donnelly.prjc runBOLT_versionNumber-LMM.sh

