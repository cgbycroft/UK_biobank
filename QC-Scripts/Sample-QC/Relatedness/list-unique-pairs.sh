# Creates the file with the list of unique batch-pairs required for running submit-pt1-kinship-in-pairs.sh
# Assumes 106 batches UKBiLEVE b1 to UKBiobank b95
# Variables "basedir" and "batchList" must be specified or this won't run

#outfile=b001__b022-pair_batches.txt
outfile=b1__b11-b001__b095-pair_batches.txt

rm $outfile
touch $outfile

for i in {1..106};
do \
    batch1=`head -n $i $batchList | tail -n 1`
    j=`expr $i + 1`

    while [[ $j -le 106 ]];
    do \
	batch2=`head -n $j $batchList | tail -n 1`
	echo $i $j
	
        echo "$basedir/data/ByBatch/$batch1-autosome-sampleqc $basedir/data/ByBatch/$batch2-autosome-sampleqc $basedir/data/Relatedness/$batch1-$batch2" >> $outfile
	((j = j + 1))
    done;
done;


echo "should be 106 choose 2 = 5,565 pairs"

wc -l $outfile
head -n 2 $outfile
