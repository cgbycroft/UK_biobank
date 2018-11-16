#!/bin/bash
#$ -cwd
#$ -N bolt_lmm_v14s
#$ -o Logs
#$ -j y
#$ -pe shmem 3
#$ -V

# Downloaded from https://data.broadinstitute.org/alkesgroup/BOLT-LMM/downloads/ on May 12th 2016
bolt=/well/donnelly/ukbiobank_project_8874/clare/src/BOLT-LMM_v2.2/bolt
basedir=/well/ukbiobank/qcoutput.V2_QCed.sample-QC

vers=v14s
pheno=Place.of.birth.in.UK...east.co.ordinate

### this sub-run is:
# BOLT version: LMM default - let BOLT decide on effect sizes.

mkdir $basedir/data/GWAS/otherGWAS/$pheno
mkdir $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers

genoDataAll=$basedir/data/Combined/byChrom/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095-sexchroms

chr=$SGE_TASK_ID

thisChromGenoData=$genoDataAll-$chr

# chr: 23 = X; 24 = Y; 25 = XPAR; 26 = MT.
# sex chromosomes are subset in subset-sex-chroms.sh

phenoFile=$basedir/data/GWAS/otherGWAS/PhenotypesForBOLT-$vers-$chr.txt  # using White british subset. Note: some individuals will be excluded in the genotype data

# phenoFile and snpSet was created using subsetPhenotypes.R
snpSet=$basedir/data/PCA/b1__b11-b001__b095-autosome-sampleqc-fastpca-White_British-highquality.prune.in  # just use SNPs from original PCA :) NOT ACTUALLY USED IN LREG!!  Note also that this does not include any SNPs on the sex chromosomes.


outFile=$basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers-$chr.out

# MEMORY NEEDED: 74,100,452,428 bytes = 5 C-NODES  <==== Sample QC SNPs
# MEMORY NEEDED: 96,739,226,940  = 7 C-NODES   <==== Oxford QC SNPs
# ALL CHROMOSOMES AT ONCE!!

# MUST MAKE SURE THAT THE .fam files for the autosomes and sex chroms are the same (they are as of May 16th)!!

# NOTE: This version uses 20 PCs as covariates - for LMM version as well!
# To adjust default snp missing rate threshold (0.1) use --maxMissingPerSnp. e.g --maxMissingPerSnp=0


# X chromosome model:
# assume full inactivation (females: <0,0.5,1>, males: <0,1>). This is encoded in the plink data.

# PAR-X chromosome model: as with autosomes
# Y model: males: <0,1>, and females excluded.
# MT model: all samples: <0,1>

# NOTE: Compute only LREG in this version (so can just input the single chromosome)

# ALSO:  This version overrides BOLT's missing values defaults, and lets everything through

$bolt \
--fam=$thisChromGenoData.fam \
--bed=$thisChromGenoData.bed \
--bim=$thisChromGenoData-altcode.bim \
--phenoFile=$phenoFile \
--phenoCol=$pheno \
--covarFile=$phenoFile \
--covarCol=Array \
--qCovarCol=Age.when.attended.assessment.centre \
--covarCol=Inferred.Gender \
--qCovarCol=PC{1:20} \
--numThreads=3 \
--statsFile=$outFile \
--covarMaxLevels=111 \
--maxMissingPerSnp=1 \
--geneticMapFile=/well/donnelly/ukbiobank_project_8874/clare/src/BOLT-LMM_v2.2/tables/genetic_map_hg19.txt.gz \
--verboseStats


# qsub -t 23-26 -q short.qc -P donnelly.prjc runBOLT_v14s-LMM.sh
# change chromosome number in output back to sex chroms!

mv $outFile $outFile.orig
cat <(head -n 1 $outFile.orig) <( awk -v c=$chr '{OFS="\t"} {$2=c}1' $outFile.orig | tail -n +2) > $outFile

