#!/bin/bash
#$ -o Logs
#$ -j y
#$ -cwd
#$ -V



basedir=/well/ukbiobank/qcoutput.V2_QCed.sample-QC
plink=/well/ukbiobank/qcoutput/Src/plink-v1.90b3k-x86_64/plink


kingRunPrefix=b1__b11-b001__b095-pair_batches.filtered

echo $kingRunPrefix

# The genotype data to use
genoTrawData=$basedir/data/Relatedness/${kingRunPrefix}-duplicates-genotypes.traw


mkdir $basedir/data/Relatedness/duplicateDiffs

columnsFile=$basedir/QC-Scripts/Sample-QC/Relatedness/${kingRunPrefix}-duplicates-twins.columnsInGenotypeFile.txt


col1=`head -n $SGE_TASK_ID $columnsFile | tail -n 1 | cut -f3 -d' '`
col2=`head -n $SGE_TASK_ID $columnsFile | tail -n 1 | cut -f4 -d' '`

id1=`head -n $SGE_TASK_ID $columnsFile | tail -n 1 | cut -f1 -d' '`
id2=`head -n $SGE_TASK_ID $columnsFile | tail -n 1 | cut -f2 -d' '`

outfile=$basedir/data/Relatedness/duplicateDiffs/${kingRunPrefix}.genoDiffs.${id1}..${id2}.txt
outfile2=$basedir/data/Relatedness/duplicateDiffs/${kingRunPrefix}.genoDiffs.${id1}..${id2}.allele2.txt

echo $col1
echo $col2
echo $id1
echo $id2
echo $outfile
echo $outfile2

# Just print a list of the rows in the genotypes file that are discordant, or missing in one or two.


#awk -v c1=$col1 -v c2=$col2 '{ if( $c1=="NA" && $c2=="NA" ) print NR, "2"; else if ( $c1=="NA" || $c2=="NA" ) print NR, "1"; else if( $c1!=$c2 ) print NR, "0" }' <( tail -n +2 $genoTrawData ) > $outfile

awk -v c1=$col1 -v c2=$col2 '{ if( $c1=="1" || $c1=="2" || $c2=="1" || $c2=="2" ) print NR }' <( tail -n +2 $genoTrawData ) > $outfile2


head $outfile
wc -l $outfile


echo DONE
