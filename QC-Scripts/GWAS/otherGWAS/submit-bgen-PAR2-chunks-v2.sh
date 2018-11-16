#!/bin/bash
#$ -cwd
#$ -N bgen-PAR2-chunks
#$ -o Logs
#$ -j y
#$ -V

#$ -t 1

basedir=/well/ukbiobank/qcoutput.V2_QCed.sample-QC
#qctool=/well/ukbiobank/qcoutput.V2_QCed.sample-QC/src/qctool_v2.0-beta4-linux-x86_64/qctool
#qctool=/tmp/for_jonathan/qctool_v2.0-beta6
qctool=/well/ukbiobank/qcoutput.V2_QCed.sample-QC/src/qctool_v2.0-beta6
bgenix=/apps/well/bgenix/20160708/bgenix

impDataX=$basedir/data/imputed/hrc+uk10k_sorted_8bit_rsids_chrPAR2.v4.v1_1 # NOTE: this is a symlink made in subset-sex-chroms.sh
outDataX=$basedir/data/imputed/chrPAR2.v2

# values taken from $outDataX.snpInfo.txt
chunks=`seq 154920000 500000 155270000 `
nChunks=`echo $chunks  | tr ' ' '\n' | wc -l`

i=$SGE_TASK_ID

chunk1=`echo $chunks  | tr ' ' '\n' | head -n $i |tail -n 1`
chunk2=`expr $chunk1 + 499999`
chunk=$chunk1-$chunk2


#chunk=2600000-2701000
echo $chunk

echo "$bgenix -g $impDataX.bgen -incl-range X:$chunk > $outDataX.$i.bgen"

$bgenix -g $impDataX.bgen -incl-range PAR2:$chunk > $outDataX.$i.bgen

# get the list of chrom/snp information
$bgenix -g $outDataX.$i.bgen -incl-range PAR2:$chunk -list > $outDataX.$i.snpInfo.txt

# construct new set
nrow=`wc -l $outDataX.$i.snpInfo.txt | cut -f1 -d' '`
nrow=`expr $nrow - 1`

awk '{print $1,$2,$3,$4,$6,$7,$1,$2,"1",$4,$6,$7}' $outDataX.$i.snpInfo.txt | head -n $nrow | tail -n +2 > $outDataX.$i.snpInfoUpdate.txt

# change the chromosome values in the chunks of bgen data.
echo "$qctool -g $outDataX.$i.bgen -map-id-data $outDataX.$i.snpInfoUpdate.txt -og $outDataX.$i.recoded.bgen -ofiletype bgen_v1.1"

$qctool -g $outDataX.$i.bgen -map-id-data $outDataX.$i.snpInfoUpdate.txt -og $outDataX.$i.recoded.bgen -ofiletype bgen_v1.1


echo "DONE!!"
