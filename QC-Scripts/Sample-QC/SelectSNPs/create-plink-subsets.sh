#!/bin/bash
# The name of the job, can be anything, simply used when displaying the list of running jobs
#$ -N sampleQC-plink-subsets
# Giving the name of the output log file
#$ -o Logs
#$ -j y
# One needs to tell the queue system to use the current directory as the working directory
# Or else the script may fail as it will execute in your top level home directory
#$ -cwd
#$ -V
#$ -P donnelly.prjc -q short.qc


################################
# For UKBiobank sample QC pipeline step 1
################################
# Clare Bycroft, May 2016
#
# This script creates the base plink files required for Sample QC. Sex chromosomes and autosomes are created independently. There will be two sets: one with all batches combined, and another for each batch separately. All will contain the same set of QC'd SNPs (sex + autosome). This script does each batch separately.
################################

plink=/well/ukbiobank/qcoutput/Src/plink-v1.90b3k-x86_64/plink

############ %%%
inputdir=/well/ukbiobank/expt/V2_QCed.SNP-QC/data # this is the directory with the plink files created immediately after applying Affy's QC to the calls files
outputdir=/well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/ByBatch
batchlist=/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/batchList.txt
sampleExclusions=$outputdir/../../QC-Scripts/referenceList.txt
############

mkdir $outputdir

Batch=`head -n ${SGE_TASK_ID} $batchlist | tail -n 1`
#i=1
#Batch=`head -n ${i} $batchlist | tail -n 1`

echo $Batch

############ Autosomes

snplist=$outputdir/../../QC-Scripts/Sample-QC/SelectSNPs/snpIncludeList-autosome.txt # This will be the intersection of UKBiLEVE and UKBiobank snps + other filters. It is created using the file: selectSNPS.R


echo $plink --bfile $inputdir/$Batch/sorted_bed/$Batch-autosome --extract $snplist --remove-fam $sampleExclusions --keep-allele-order --make-bed --out $outputdir/$Batch-autosome-sampleqc

$plink --bfile $inputdir/$Batch/sorted_bed/$Batch-autosome --extract $snplist --remove-fam $sampleExclusions --keep-allele-order --make-bed --out $outputdir/$Batch-autosome-sampleqc

rm $outputdir/$Batch-autosome-sampleqc.log

# Now check that number of snps in plink files is as expected
echo "in plink files..."
wc $outputdir/$Batch-autosome-sampleqc.bim
echo "snp list..." 
wc $snplist



############ Sex chromosomes

snplist=$outputdir/../../QC-Scripts/Sample-QC/SelectSNPs/snpIncludeList-sexchrom.txt # This will be the intersection of UKBiLEVE and UKBiobank snps + other filters. It is created using the file: selectSNPS.R

echo $plink --bfile $inputdir/$Batch/sorted_bed/$Batch-sexchrom --extract $snplist --remove-fam $sampleExclusions --keep-allele-order --make-bed --out $outputdir/$Batch-sexchrom-sampleqc

$plink --bfile $inputdir/$Batch/sorted_bed/$Batch-sexchrom --extract $snplist --remove-fam $sampleExclusions --keep-allele-order --make-bed --out $outputdir/$Batch-sexchrom-sampleqc

rm $outputdir/$Batch-sexchrom-sampleqc.log


# Now check that number of snps in plink files is as expected
echo "in plink files..."
wc $outputdir/$Batch-sexchrom-sampleqc.bim
echo "snp list..."
wc $snplist

