#!/bin/bash
# The name of the job, can be anything, simply used when displaying the list of running jobs
#$ -N OxfordQC-plink-subsets
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
# This QC is applying Oxford QC by batch, so should look like the released dataset.
# Merging is done with merge-batches-Oxford-QC.sh
################################

plink=/well/ukbiobank/qcoutput/Src/plink-v1.90b3k-x86_64/plink
basedir=/well/ukbiobank/qcoutput.V2_QCed.sample-QC
############ %%%
inputdir=/well/ukbiobank/expt/V2_QCed.SNP-QC/data # this is the directory with the plink files created immediately after applying Affy's QC to the calls files
outputdir=$basedir/data/ByBatch
batchlist=$basedir/QC-Scripts/batchList.txt
sampleExclusions=$basedir/QC-Scripts/referenceList.txt
############

mkdir $outputdir

Batch=`head -n ${SGE_TASK_ID} $batchlist | tail -n 1`
#i=1
#Batch=`head -n ${i} $batchlist | tail -n 1`

echo $Batch

############ FIND LIST OF SNPS FOR EACH BATCH

# NOTE: there is a separate SNP exclude list by batch. NOTE: malefemale test is applied across whole of X chromosome. Don't really want this! FIX LATER...%%%%%%

touch $basedir/QC-Scripts/Sample-QC/SelectSNPs/snpExcludeList-oxfordqc-$Batch-autosome.txt
touch $basedir/QC-Scripts/Sample-QC/SelectSNPs/snpExcludeList-oxfordqc-$Batch-sexchrom.txt

for Test in {"plate_effect","batch_effect","malefemale","hwe"};
do \
echo $Test

awk -v b=$Batch '$1 == b {print $2}' /well/ukbiobank/expt/V2_QCed.SNP-QC/data/SNP-QC_summary/V2_QCed.autosome.$Test.pval_lt_1e-12.UKBiLEVEAX_b1-b11.Batch_b001-b095.txt >> $basedir/QC-Scripts/Sample-QC/SelectSNPs/snpExcludeList-oxfordqc-$Batch-autosome.txt

awk -v b=$Batch '$1 == b {print $2}' /well/ukbiobank/expt/V2_QCed.SNP-QC/data/SNP-QC_summary/V2_QCed.sexchrom.$Test.pval_lt_1e-12.UKBiLEVEAX_b1-b11.Batch_b001-b095.txt >> $basedir/QC-Scripts/Sample-QC/SelectSNPs/snpExcludeList-oxfordqc-$Batch-sexchrom.txt

done;


############ Autosomes

snplist=$basedir/QC-Scripts/Sample-QC/SelectSNPs/snpExcludeList-oxfordqc-$Batch-autosome.txt

echo $plink --bfile $inputdir/$Batch/sorted_bed/$Batch-autosome --exclude $snplist --remove-fam $sampleExclusions --keep-allele-order --make-bed --out $outputdir/$Batch-autosome-oxfordqc

$plink --bfile $inputdir/$Batch/sorted_bed/$Batch-autosome --exclude $snplist --remove-fam $sampleExclusions --keep-allele-order --make-bed --out $outputdir/$Batch-autosome-oxfordqc

rm $outputdir/$Batch-autosome-oxfordqc.log

# Now check that number of snps in plink files is as expected
echo "in plink files..."
wc $outputdir/$Batch-autosome-oxfordqc.bim
echo "snp list..." 
wc $snplist



############ Sex chromosomes

snplist=$basedir/QC-Scripts/Sample-QC/SelectSNPs/snpExcludeList-oxfordqc-$Batch-sexchrom.txt

echo $plink --bfile $inputdir/$Batch/sorted_bed/$Batch-sexchrom --exclude $snplist --remove-fam $sampleExclusions --keep-allele-order --make-bed --out $outputdir/$Batch-sexchrom-oxfordqc

$plink --bfile $inputdir/$Batch/sorted_bed/$Batch-sexchrom --exclude $snplist --remove-fam $sampleExclusions --keep-allele-order --make-bed --out $outputdir/$Batch-sexchrom-oxfordqc

rm $outputdir/$Batch-sexchrom-oxfordqc.log


# Now check that number of snps in plink files is as expected
echo "in plink files..."
wc $outputdir/$Batch-sexchrom-oxfordqc.bim
echo "snp list..."
wc $snplist

