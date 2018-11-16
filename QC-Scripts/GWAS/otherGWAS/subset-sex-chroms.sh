########## 
# script to subset plink files to make chromosome by chromosome files for GWAS. (Clare's script)

basedir=/well/ukbiobank/qcoutput.V2_QCed.sample-QC
plink=/well/ukbiobank/qcoutput/Src/plink-v1.90b3k-x86_64/plink
#qctool=/apps/well/qctool/1.5/qctool
qctool=/well/ukbiobank/qcoutput.V2_QCed.sample-QC/src/qctool_v2.0-beta4-linux-x86_64/qctool
bgenix=/apps/well/bgenix/20160708/bgenix

# based data from Colin
#genoDataAll=/well/ukbiobank/expt/V2_QCed.export/data/calls_export/v2/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095
genoDataAll=/well/ukbiobank/expt/V2_QCed.export/data/calls_export/v3/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095
outdata=$basedir/data/Combined/byChrom/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095-sexchroms_v3


# subset by chromosome and recode for particular Xinactivation model

# separate X chromosome 
$plink --bfile $genoDataAll --chr 23 --split-x hg19 --make-bed --out $outdata-X


# X - non-PAR - note: males must be coded 0 or 2 (this is the case, I believe)
$plink --bfile $outdata-X --chr 23 --make-bed --out $outdata-23
awk 'OFS="\t" {$1=1}1' $outdata-23.bim > $outdata-23-altcode.bim

# X PAR 
$plink --bfile $outdata-X --chr 25 --make-bed --out $outdata-25
#==> recode chromosome to 1
awk 'OFS="\t" {$1=1}1' $outdata-25.bim > $outdata-25-altcode.bim

# Y chromosome (codes are 0,2)
$plink --bfile $genoDataAll --chr 24 --make-bed --out $outdata-24
#==> recode chromosome to 1
awk 'OFS="\t" {$1=1}1' $outdata-24.bim > $outdata-24-altcode.bim

# MT chromosome (codes are 0,2)
$plink --bfile $genoDataAll --chr 26 --make-bed --out $outdata-26
#==> recode chromosome to 1
awk 'OFS="\t" {$1=1}1' $outdata-26.bim > $outdata-26-altcode.bim


# check that males are coded 0,2. Yep.
$plink --bfile $outdata-23 --recode AD --out $outdata-23



# replace FID with IID in column 1 of fam files (in order to match the phenotype file!)
mv $outdata-X.fam $outdata-X.tmp
paste -d' ' <(cut -f2 -d' ' $outdata-X.tmp) <(cut -f2-6 -d' ' $outdata-X.tmp) > $outdata-X.fam

mv $outdata-23.fam $outdata-23.tmp
paste -d' ' <(cut -f2 -d' ' $outdata-23.tmp) <(cut -f2-6 -d' ' $outdata-23.tmp) > $outdata-23.fam

mv $outdata-24.fam $outdata-24.tmp
paste -d' ' <(cut -f2 -d' ' $outdata-24.tmp) <(cut -f2-6 -d' ' $outdata-24.tmp) > $outdata-24.fam

mv $outdata-25.fam $outdata-25.tmp
paste -d' ' <(cut -f2 -d' ' $outdata-25.tmp) <(cut -f2-6 -d' ' $outdata-25.tmp) > $outdata-25.fam

mv $outdata-26.fam $outdata-26.tmp
paste -d' ' <(cut -f2 -d' ' $outdata-26.tmp) <(cut -f2-6 -d' ' $outdata-26.tmp) > $outdata-26.fam



# The above is only relevant for v2 of the input genotype data, which is not separated by chromosome.
# Otherwise just copy over into my folders.
cp $genoDataAll.chrX.bim $outdata-23.bim
cp $genoDataAll.chrX.bed $outdata-23.bed
cp $genoDataAll.chrX.fam $outdata-23.fam
awk 'OFS="\t" {$1=1}1' $outdata-23.bim > $outdata-23-altcode.bim

cp $genoDataAll.chrXY.bim $outdata-25.bim
cp $genoDataAll.chrXY.bed $outdata-25.bed
cp $genoDataAll.chrXY.fam $outdata-25.fam
awk 'OFS="\t" {$1=1}1' $outdata-25.bim > $outdata-25-altcode.bim

cp $genoDataAll.chrY.bim $outdata-24.bim
cp $genoDataAll.chrY.bed $outdata-24.bed
cp $genoDataAll.chrY.fam $outdata-24.fam
awk 'OFS="\t" {$1=1}1' $outdata-24.bim > $outdata-24-altcode.bim

cp $genoDataAll.chrMT.bim $outdata-26.bim
cp $genoDataAll.chrMT.bed $outdata-26.bed
cp $genoDataAll.chrMT.fam $outdata-26.fam
awk 'OFS="\t" {$1=1}1' $outdata-26.bim > $outdata-26-altcode.bim

### ===> note that snp ids are affy id or rsid where possible.



# exclude females from Y chromosome data (so that SNP missing rates aren't artificially low!)
mv $outdata-24.fam $outdata-24-all.fam
mv $outdata-24.bim $outdata-24-all.bim
mv $outdata-24.bed $outdata-24-all.bed
mv $outdata-24-altcode.bim $outdata-24-altcode-all.bim

$plink --bfile $outdata-24-all --filter $outdata-24-all.fam 1 -mfilter 3 --make-bed --out $outdata-24
cp $outdata-24-altcode-all.bim $outdata-24-altcode.bim

wc $outdata-24.fam


################
# Do this for the X chromosome (non-PAR) for the imputed data. <======== first release

cd /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS

impDataX=/well/ukbiobank/imputation/final/full/bgen/hrc+uk10k_sorted_8bit_rsids/chrX
outDataX=$basedir/data/imputed/chrX

# list of all the snps in the X chromosome bgen files
$bgenix -g $outDataX.bgen -incl-range X:0- -list > $outDataX.snpInfo.txt

# Chunk up the bgen files into ~15000 snps
qsub -p donnelly.prjc -q short.qc submit-bgen-X-chunks.sh

# copy over the sample file
cp $impDataX.sample $outDataX.sample

# add sex to the sample file!
#If sex information is in the .sample file, it must be in a column titled 'sex' (capital letters ok) of type 'D' (discrete covariate), and be coded in the usual 1=male/2=female/0=unknown manner, to be loaded.
Rscript addSex.R $outDataX.sample

mv $outDataX.sample $outDataX.sample.raw
mv $outDataX.sample.sex $outDataX.sample



################
# Do this for the X chromosome (PAR and non-PAR) for the imputed data. <======== second release (global v2)

cd /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS

################ For the X
impDataX=/well/ukbiobank/expt/full_release_issues/uk10k_imputation_annotation/hrc+uk10k_sorted_8bit_rsids_chrX.v4.v1_1.bgen
impDataSampleX=/well/ukbiobank/imputation/final/full/bgen/hrc+uk10k_sorted_8bit_rsids/v2/chrX.v2.sample

# Make a sym-link to the real data
ln -s $impDataX $basedir/data/imputed/hrc+uk10k_sorted_8bit_rsids_chrX.v4.v1_1.bgen

outDataX=$basedir/data/imputed/chrX.v2

# list of all the snps in the X chromosome bgen files (this tells us the positions at the start and end of the par regions.)
$bgenix -g $basedir/data/imputed/hrc+uk10k_sorted_8bit_rsids_chrX.v4.v1_1.bgen -incl-range X:0- -list > $outDataX.snpInfo.txt

# Chunk up the bgen files into ~15000 snps
qsub -P donnelly.prjc -q short.qc submit-bgen-X-chunks-v2.sh

# copy over the sample file
cp $impDataSampleX $outDataX.sample

# add sex to the sample file!
#If sex information is in the .sample file, it must be in a column titled 'sex' (capital letters ok) of type 'D' (discrete covariate), and be coded in the usual 1=male/2=female/0=unknown manner, to be loaded.
Rscript addSex.R $outDataX.sample

mv $outDataX.sample $outDataX.sample.raw
mv $outDataX.sample.sex $outDataX.sample


################ For PAR1
impDataX=/well/ukbiobank/expt/full_release_issues/uk10k_imputation_annotation/hrc+uk10k_sorted_8bit_rsids_chrPAR1.v4.v1_1.bgen
impDataSampleX=/well/ukbiobank/imputation/final/full/bgen/hrc+uk10k_sorted_8bit_rsids/v2/chrPAR.sample

# Make a sym-link to the real data
ln -s $impDataX $basedir/data/imputed/hrc+uk10k_sorted_8bit_rsids_chrPAR1.v4.v1_1.bgen

outDataX=$basedir/data/imputed/chrPAR1.v2

# list of all the snps in the X chromosome bgen files (this tells us the positions at the start and end of the par regions.)
$bgenix -g $basedir/data/imputed/hrc+uk10k_sorted_8bit_rsids_chrPAR1.v4.v1_1.bgen -incl-range PAR1:0- -list > $outDataX.snpInfo.txt

# Chunk up the bgen files into ~15000 snps (500Mb)
qsub -P donnelly.prjc -q short.qc submit-bgen-PAR1-chunks-v2.sh

# copy over the sample file
cp $impDataSampleX $outDataX.sample

# add sex to the sample file!
#If sex information is in the .sample file, it must be in a column titled 'sex' (capital letters ok) of type 'D' (discrete covariate), and be coded in the usual 1=male/2=female/0=unknown manner, to be loaded.


################ For PAR2
impDataX=/well/ukbiobank/expt/full_release_issues/uk10k_imputation_annotation/hrc+uk10k_sorted_8bit_rsids_chrPAR2.v4.v1_1.bgen
impDataSampleX=/well/ukbiobank/imputation/final/full/bgen/hrc+uk10k_sorted_8bit_rsids/v2/chrPAR.sample

outDataX=$basedir/data/imputed/chrPAR2.v2

# Make a sym-link to the real data
ln -s $impDataX $basedir/data/imputed/hrc+uk10k_sorted_8bit_rsids_chrPAR2.v4.v1_1.bgen

# list of all the snps in the X chromosome bgen files (this tells us the positions at the start and end of the par regions.)
$bgenix -g $basedir/data/imputed/hrc+uk10k_sorted_8bit_rsids_chrPAR2.v4.v1_1.bgen -incl-range PAR2:0- -list > $outDataX.snpInfo.txt

# Chunk up the bgen files into ~15000 snps (500Mb). JUST ONE CHUNK FOR PAR2
qsub -P donnelly.prjc -q short.qc submit-bgen-PAR2-chunks-v2.sh

# copy over the sample file
cp $impDataSampleX $outDataX.sample



