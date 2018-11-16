#!/bin/bash
#$ -cwd
#$ -N interim-compare-count-differences
#$ -o Logs
#$ -j y
#$ -V
#$ -P donnelly.prjc
#$ -q short.qc
#$ -t 1-33

basedir=/well/ukbiobank/qcoutput.V2_QCed.sample-QC
plink=/well/ukbiobank/qcoutput/Src/plink-v1.90b3k-x86_64/plink

interimDir=/well/ukbiobank/expt/check_data_release_1_ctsu/data/calls
newDataDir=/well/ukbiobank/expt/V2_QCed.export/data/calls_export/v3
outDir=$basedir/data/Combined/interimReleaseSamples/Comparison

# chromosomes
#i=22

# batch
batchList=$basedir/QC-Scripts/batchList.txt
#b=1
b=$SGE_TASK_ID

# which array
if [[ $b -le 11 ]]
then
    array=UKBL
    batch=UKBiLEVEAX_b$b;

fi;
   
if [[ $b -ge 12 ]]
then
    b=`expr $b - 11`
    array=UKB
    batch=Batch_b`printf "%03d" $b` 
fi;


# Loop over chromosomes
chroms=`cat <(seq 1 22) <(printf "X\nY\nMT\n")`
#chroms=`printf "X\nY\nXY\nMT\n"`

for i in $chroms;
do \

    # PAR is separated in the new release data, but not in the interim release.
    if [[ "$i" == "XY" ]]; then dataFileInterim=$interimDir/${array}_$batch.chrX.calls.ctsu_basic.txt; fi

    if [[ "$i" != "XY" ]]; then dataFileInterim=$interimDir/${array}_$batch.chr$i.calls.ctsu_basic.txt; fi

    # Count the differences for this chromosome, and these individuals
    Rscript countDifferences.R \
	    $basedir/data/Combined/interimReleaseSamples/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.chr$i.interim_samples_$batch.traw \
	    $outDir/$batch.$i.new.best.array.txt \
	    $dataFileInterim \
	    $outDir/$batch.$i.interim.best.array.txt \
	    $outDir/${array}_$batch.chr$i.calls.ctsu_basic.compareToNew2


done;
    
#rm $batch.$i.interim.best.array.txt
#rm $batch.$i.new.best.array.txt
#rm $batch.$i.interim.snp.order.txt
#rm $batch.$i.new.snp.order.txt
#rm $basedir/data/Combined/interimReleaseSamples/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.chr$i.interim_samples_$batch.traw
