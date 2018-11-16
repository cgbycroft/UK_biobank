#!/bin/bash
#$ -cwd
#$ -N extract-sex-CNV
#$ -o Logs
#$ -j y
#$ -V


########
# shell script to extract sex chromosome CNV values

basedir=/well/ukbiobank/qcoutput.V2_QCed.sample-QC
batchList=$basedir/QC-Scripts/batchList.txt

batch=`head -n $SGE_TASK_ID $batchList | tail -n 1`

echo $SGE_TASK_ID
echo $batch

#batch=Batch_b001
array=`grep Batch <( echo $batch )`
if [ "$array" == "$batch" ]
then
    snporder=/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/Sample-QC/ExtractCNV/Batch_b001-CNV-snp-order.txt
else
    snporder=/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/Sample-QC/ExtractCNV/UKBiLEVEAX_b1-CNV-snp-order.txt
fi;

mkdir $basedir/data/CNV/$batch

# sex chromosomes start at 823236 and end at the bottom of the file
row=`awk '{print $2}' $snporder | grep -w -n X | cut -f1 -d: | sort | head -n 1`


# log2ratios
dataFile=/well/ukbiobank/source/V2_All/$batch/AxiomGT1.cnv.log2ratio.final.txt.gz
outFile=$basedir/data/CNV/$batch/AxiomGT1.cnv.log2ratio.final.sexchrom

echo cat <( zcat $dataFile | head -n 1 ) <( zcat $dataFile | tail -n +$row ) > $outFile.txt
time cat <( zcat $dataFile | head -n 1 ) <( zcat $dataFile | tail -n +$row ) > $outFile.txt


# BAFs
dataFile=/well/ukbiobank/source/V2_All/$batch/AxiomGT1.cnv.baf.final.txt.gz
outFile=$basedir/data/CNV/$batch/AxiomGT1.cnv.baf.final.sexchrom

echo cat <( zcat $dataFile | head -n 1 ) <( zcat $dataFile | tail -n +$row ) > $outFile.txt
time cat <( zcat $dataFile | head -n 1 ) <( zcat $dataFile | tail -n +$row ) > $outFile.txt
