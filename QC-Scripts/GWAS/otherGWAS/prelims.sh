
basedir=/well/ukbiobank/qcoutput.V2_QCed.sample-QC
plink=/well/ukbiobank/qcoutput/Src/plink-v1.90b3k-x86_64/plink

pheno="Standing.height"
pheno="Place.of.birth.in.UK...north.co.ordinate"
pheno="Place.of.birth.in.UK...east.co.ordinate"
pheno="Body.mass.index..BMI."
pheno="Red.Blood.Cell.Count"

###############
# get latest GWAS catalogue data in the correct build. It will download it again if there is an updated version. For the time-stamp of the file, see 
Rscript get-gwas-catalogue.R hg19 $basedir/QC-Scripts/GWAS/otherGWAS/GWAScatalogue

# get cytobands for plots (build hg19)
Rscript get-cytobands.R hg19 $basedir/QC-Scripts/GWAS/otherGWAS/cytobands

# get gene lists
Rscript get-genes.R hg19 $basedir/QC-Scripts/GWAS/otherGWAS/genes


###############
# All samples (some ad-hoc QC + unrelated); not quantile normalised traits
vers=v1
phenotypes=`echo Standing.height Place.of.birth.in.UK...north.co.ordinate Place.of.birth.in.UK...east.co.ordinate Body.mass.index..BMI. Red.Blood.Cell.Count Haemoglobin.Concentration Haematocrit.percentage Mean.Corpuscular.Volume Mean.Corpuscular.Haemoglobin Mean.Corpuscular.Haemoglobin.Concentration Intra.ocular.pressure..corneal.compensated..right. Intra.ocular.pressure..Goldmann.correlated..right. Intra.ocular.pressure..corneal.compensated..left. Intra.ocular.pressure..Goldmann.correlated..left. Home.area.population.density...urban.or.rural Internal.Pico..ng.uL. CV. dQC Cluster.CR Array.as.binary`
###############
###############
# White British samples, QC + unrelated + quantile normalised traits. STILL NEED TO RUN THIS WITH NEW SET OF WHITE-BRITISH, and A NEW SET OF PCs BASED ON THAT SUBSET
vers=v2  # LMM only
phenotypes=`echo Standing.height Place.of.birth.in.UK...north.co.ordinate Place.of.birth.in.UK...east.co.ordinate Body.mass.index..BMI. Red.Blood.Cell.Count Haemoglobin.Concentration Haematocrit.percentage Mean.Corpuscular.Volume Mean.Corpuscular.Haemoglobin Mean.Corpuscular.Haemoglobin.Concentration Intra.ocular.pressure..corneal.compensated..right. Intra.ocular.pressure..Goldmann.correlated..right. Intra.ocular.pressure..corneal.compensated..left. Intra.ocular.pressure..Goldmann.correlated..left. Intra.ocular.pressure..Goldmann.correlated..mean Intra.ocular.pressure..corneal.compensated..mean Internal.Pico..ng.uL. CV. dQC Cluster.CR Array.as.binary Standing.height.qnorm Body.mass.index..BMI..qnorm Red.Blood.Cell.Count.qnorm Haemoglobin.Concentration.qnorm Haematocrit.percentage.qnorm Mean.Corpuscular.Volume.qnorm Mean.Corpuscular.Haemoglobin.qnorm Mean.Corpuscular.Haemoglobin.Concentration.qnorm Intra.ocular.pressure..corneal.compensated..right..qnorm Intra.ocular.pressure..Goldmann.correlated..right..qnorm Intra.ocular.pressure..corneal.compensated..left..qnorm Intra.ocular.pressure..Goldmann.correlated..left..qnorm Intra.ocular.pressure..Goldmann.correlated..mean.qnorm Intra.ocular.pressure..corneal.compensated..mean.qnorm Internal.Pico..ng.uL..qnorm CV..qnorm dQC.qnorm Cluster.CR.qnorm`
###############
###############
# White British samples, QC + unrelated ===> JUST HEIGHT
phenotypes=`echo Standing.height`
versions=`echo v2 v2a v2b v2c v2d`
vers=v2a # 1000 samples
vers=v2b # 10000 samples
vers=v2c # 50000 samples
vers=v2d # 150000 samples
# 500000 samples (or the white-british subset) are provided in v2
versions=`echo v3 v3a v3b v3c v3d`
###############
vers=v3  # As v2 but LMM + PCs as covariates + white British
# Will also run conditional analysis for Standing.height
###############
###############
vers=v4  # As v3 but dQC + Cluster.CR as covariate, only height & only LREG
###############
vers=v5  # As v3 but using Colin's input data ==> genotype_gwas folder. Contains about 500 more snps than final release (multi-allelics!).
phenotypes=`echo Intra.ocular.pressure..Goldmann.correlated..mean Standing.height Place.of.birth.in.UK...north.co.ordinate Place.of.birth.in.UK...east.co.ordinate Red.Blood.Cell.Count Haemoglobin.Concentration dQC Cluster.CR Array.as.binary Internal.Pico..ng.uL.`
###############
vers=v6  # As v3 but using imputed data chromosome 2
phenotypes=`echo Intra.ocular.pressure..Goldmann.correlated..mean Standing.height Place.of.birth.in.UK...north.co.ordinate Place.of.birth.in.UK...east.co.ordinate Red.Blood.Cell.Count Haemoglobin.Concentration dQC Cluster.CR Array.as.binary Internal.Pico..ng.uL.`
vers=v6a  # As v6 but including LMM option (to get BETAs!)
###############
# v7-10 done by Colin
###############
vers=v11  # As v5 but sex chromosomes only, and v2 of the genotype calls; and remove BOLT's 10% cuttoff
phenotypes=`echo Body.mass.index..BMI. Intra.ocular.pressure..Goldmann.correlated..mean Standing.height Place.of.birth.in.UK...north.co.ordinate Place.of.birth.in.UK...east.co.ordinate Red.Blood.Cell.Count Haemoglobin.Concentration dQC Cluster.CR Array.as.binary Internal.Pico..ng.uL.`
# chr: 23 = X; 24 = Y; 25 = XPAR; 26 = MT.
# NOTE:sex chromosomes were subset in subset-sex-chroms.sh
###############
# v12 is chromosome 20 imputed data (Colin)
###############
vers=v13  # As v6a, but with chromosomes from imputed data from early november 
phenotypes=`echo Intra.ocular.pressure..Goldmann.correlated..mean Standing.height Place.of.birth.in.UK...north.co.ordinate Place.of.birth.in.UK...east.co.ordinate Red.Blood.Cell.Count Haemoglobin.Concentration dQC Cluster.CR Array.as.binary Internal.Pico..ng.uL.`
###############
vers=v14  # As v5, but with release genotype data from calls_export/v2 . Also includes LDscores for lmm mode, and maxmissing per snp set to 1.
phenotypes=`echo Intra.ocular.pressure..Goldmann.correlated..mean Standing.height Place.of.birth.in.UK...north.co.ordinate Place.of.birth.in.UK...east.co.ordinate Red.Blood.Cell.Count Haemoglobin.Concentration dQC Cluster.CR Array.as.binary Internal.Pico..ng.uL.`
vers=v14s  # As v14 but with SEX CHROMOSOMES. Uses same input genotypes as v11, but missing rates filters set to 1. Note lmm option is not possible because of need to rename chromosomes!
###############
vers=v15  # As v13 but with the final *.bgen files as converted by Colin to v.1 Bgen format. Only chromosome 2 is done here. ---> All other autosomes done now actually.
phenotypes=`echo Intra.ocular.pressure..Goldmann.correlated..mean Standing.height Place.of.birth.in.UK...north.co.ordinate Place.of.birth.in.UK...east.co.ordinate Red.Blood.Cell.Count Haemoglobin.Concentration dQC Cluster.CR Array.as.binary Internal.Pico..ng.uL.`
phenotypes2=`echo Intra.ocular.pressure..Goldmann.correlated..mean.qnorm Standing.height.qnorm Red.Blood.Cell.Count.qnorm Haemoglobin.Concentration.qnorm dQC.qnorm Cluster.CR.qnorm Internal.Pico..ng.uL..qnorm`
###############
vers=v15s  # As v15 but X chromosome only.
###############
vers=v16  # As v14 but with the v3 genotype calls files. Autosomes only.
phenotypes=`echo Intra.ocular.pressure..Goldmann.correlated..mean Standing.height Standing.height.qnorm Place.of.birth.in.UK...north.co.ordinate Place.of.birth.in.UK...east.co.ordinate Red.Blood.Cell.Count Haemoglobin.Concentration dQC Cluster.CR Array.as.binary Internal.Pico..ng.uL. Intra.ocular.pressure..Goldmann.correlated..mean.qnorm Standing.height.qnorm Red.Blood.Cell.Count.qnorm Haemoglobin.Concentration.qnorm dQC.qnorm Cluster.CR.qnorm Internal.Pico..ng.uL..qnorm`
###############
vers=v16s  # As v16 but X chromosome only.
###############
vers=v17  # As v15 but a random phenotype ~N(0,1); ~LN(0,1)
phenotypes=`echo Rnorm.0.1 Rlnorm.0.1`
###############
vers=v18  # As v16 but a random phenotype ~N(0,1); ~LN(0,1)
phenotypes=`echo Rnorm.0.1 Rlnorm.0.1`
###############
vers=v19  # As v15 but with the version 2 *.bgen files as converted by Colin to v.1 Bgen format. Only chromosome 2 is done here.  8/8/2018: Clare running other autosomes for Standing Height
phenotypes=`echo Intra.ocular.pressure..Goldmann.correlated..mean Intra.ocular.pressure..Goldmann.correlated..mean.qnorm Standing.height Place.of.birth.in.UK...north.co.ordinate Place.of.birth.in.UK...east.co.ordinate Red.Blood.Cell.Count Haemoglobin.Concentration dQC Cluster.CR Array.as.binary Internal.Pico..ng.uL.`
#phenotypes2=`echo Standing.height.qnorm Red.Blood.Cell.Count.qnorm Haemoglobin.Concentration.qnorm dQC.qnorm Cluster.CR.qnorm Internal.Pico..ng.uL..qnorm`
###############
###############
vers=v19s  # As v19 but X chromosome only. 
###############
###############
vers=v19b  # As v19 but randomly sample 10,000 individuals. <== only chrom2
###############
###############
vers=v16b  # As v16 but using the same random sample as v19b.
###############
###############
vers=v20  # As v19 but using data converted to bgen 1.1 by colin in August 2018
pheno=Standing.height # <==== STILL TO RUN!!!
###############
###############
vers=v20b  # As v20 but randomly sample 10,000 individuals. # <==== STILL TO RUN!!!
###############


# only use this if vers v1,v2, or greater than v7; not v13
#phenotypeFile=$basedir/data/GWAS/otherGWAS/PhenotypesForBOLT-$vers.txt

## Make sure that the phenotypes are actually in the phenotypeFile.
phenotypeFile=$basedir/data/GWAS/otherGWAS/PhenotypesForBOLT-$vers.txt
for pheno in $phenotypes;
do \
grep -w $pheno <(head -n 1 $phenotypeFile) | wc -l
done;

X=X;
X=PAR1
X=PAR2
for pheno in $phenotypes;
do \

#    for vers in $versions;
#    do \

	echo $pheno

	mkdir $basedir/QC-Scripts/GWAS/otherGWAS/$pheno
	mkdir $basedir/QC-Scripts/GWAS/otherGWAS/$pheno/Logs
	mkdir $basedir/QC-Scripts/GWAS/otherGWAS/$pheno/plots

	#cp $basedir/QC-Scripts/GWAS/otherGWAS/runBOLT_$vers-LREG-sex.sh $basedir/QC-Scripts/GWAS/otherGWAS/$pheno

	#sed "s/phenotypeForThisFile/$pheno/g" $basedir/QC-Scripts/GWAS/otherGWAS/$pheno/runBOLT_$vers-LREG-sex.sh | sed "s/versionNumber/$vers/g" > $basedir/QC-Scripts/GWAS/otherGWAS/$pheno/runBOLT_$vers-LREG-sex-$pheno.sh

	#rm $basedir/QC-Scripts/GWAS/otherGWAS/$pheno/runBOLT_$vers-LREG-sex.sh

	cp $basedir/QC-Scripts/GWAS/otherGWAS/runBOLT_$vers-LMM.sh $basedir/QC-Scripts/GWAS/otherGWAS/$pheno
	#cp $basedir/QC-Scripts/GWAS/otherGWAS/runBOLT_$vers-LMM-$X.sh $basedir/QC-Scripts/GWAS/otherGWAS/$pheno

	#cp $basedir/QC-Scripts/GWAS/otherGWAS/runBOLT_v3a-LMM.sh $basedir/QC-Scripts/GWAS/otherGWAS/$pheno/runBOLT_$vers-LMM.sh
	#cp $basedir/QC-Scripts/GWAS/otherGWAS/runBOLT_$vers-LMM-conditional.sh $basedir/QC-Scripts/GWAS/otherGWAS/$pheno

	sed "s/phenotypeForThisFile/$pheno/g" $basedir/QC-Scripts/GWAS/otherGWAS/$pheno/runBOLT_$vers-LMM.sh  | sed "s/versionNumber/$vers/g" > $basedir/QC-Scripts/GWAS/otherGWAS/$pheno/runBOLT_$vers-LMM-$pheno.sh
	#sed "s/phenotypeForThisFile/$pheno/g" $basedir/QC-Scripts/GWAS/otherGWAS/$pheno/runBOLT_$vers-LMM-conditional.sh  | sed "s/versionNumber/$vers/g" > $basedir/QC-Scripts/GWAS/otherGWAS/$pheno/runBOLT_$vers-LMM-conditional-$pheno.sh
	#sed "s/phenotypeForThisFile/$pheno/g" $basedir/QC-Scripts/GWAS/otherGWAS/$pheno/runBOLT_$vers-LMM-$X.sh  | sed "s/versionNumber/$vers/g" > $basedir/QC-Scripts/GWAS/otherGWAS/$pheno/runBOLT_$vers-LMM-$X-$pheno.sh
	
	rm $basedir/QC-Scripts/GWAS/otherGWAS/$pheno/runBOLT_$vers-LMM.sh
	#rm $basedir/QC-Scripts/GWAS/otherGWAS/$pheno/runBOLT_$vers-LMM-conditional.sh
	#rm $basedir/QC-Scripts/GWAS/otherGWAS/$pheno/runBOLT_$vers-LMM-$X.sh

    #done;
done;

#NOTE: make sure to remove Array covariate for the Array.as.binary trait!


##########################################
## Special for IMPUTATION bgen versions v6. create dummy file that has all the right samples 
# create dummy plink file with right list of samples. NOT NEEDED FOR v6a.
if[ "$vers" == "v6" ]
then \
    cp /well/ukbiobank/phasing/final/phased_chunks/chr2.test2.sample chr2.test2.sample
    cut -f1,2 -d' ' chr2.test2.sample | tail -n +3  > toKeep.tmp
    awk '$1==2' /well/ukbiobank/expt/V2_QCed.export/data/genotype_gwas/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.genotype_gwas.bim | head -n 1 > snp.tmp
    $plink --bfile /well/ukbiobank/expt/V2_QCed.export/data/genotype_gwas/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.genotype_gwas --extract snp.tmp --keep-allele-order --keep toKeep.tmp --make-bed --out $basedir/data/GWAS/otherGWAS/ImputedBgenSampleSet.$vers.dummy
fi;

if[ "$vers" == "v6a" ]
then \

    sampleFile=/well/ukbiobank/imputation/final/chr2/out/chr2.test1.r1.hrc.I4.sample
    cp $sampleFile chr2.sample
    cut -f1,2 -d' ' chr2.sample | tail -n +3  > toKeep.2.tmp    
    snpSet=$basedir/data/PCA/b1__b11-b001__b095-autosome-sampleqc-fastpca-White_British-highquality.prune.in  # just use SNPs from original PCA :) NOT ACTUALLY USED IN LREG!!  Note also that this do
    genoData=/well/ukbiobank/expt/V2_QCed.export/data/genotype_gwas/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.genotype_gwas  # <---- this is ~100 SNPs different from the calls_export/v2 files...
    $plink --bfile $genoData --extract $snpSet --keep-allele-order --keep toKeep.2.tmp --make-bed --out $basedir/data/GWAS/otherGWAS/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.genotype_gwas.ImputedBgenSampleSet.modelSnps

fi;

if[ "$vers" == "v13" ]
then \
    # Match genotypes with imputed data (just on subset of snps for LMM model)
# get list of samples in .bgen file
    sampleFile=/well/ukbiobank/imputation/final/full/bgen/chr22.hrc+uk10k.I4.v1.1.sample
    cp $sampleFile chr22.sample
    cut -f1,2 -d' ' chr22.sample | tail -n +3  > toKeep.$vers.tmp
    
# get list of snps for LMM model, in rsid format.
    snpSetbim=$basedir/data/PCA/b1__b11-b001__b095-autosome-sampleqc-fastpca-White_British-highquality-pruned.bim  # just use SNPs from original PCA :) NOT ACTUALLY USED IN LREG!!
    genoData=/well/ukbiobank/expt/V2_QCed.export/data/calls_export/v2/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095  # <---- THis is the release dataset. c.f v6a which used genotype_gwas data
    # get list of snps in rsid format
    snpSet=$basedir/data/GWAS/otherGWAS/b1__b11-b001__b095-autosome-sampleqc-fastpca-White_British-highquality-pruned.snpsetRSIDs.txt
    
    Rscript getReleaseSNPIDs.R $snpSetbim $genoData.bim $snpSet
    
    # NOTE: 2 snps missing from genotype data (missed multi-allelics)
    $plink --bfile $genoData --extract $snpSet --keep-allele-order --keep toKeep.$vers.tmp --make-bed --out $basedir/data/GWAS/otherGWAS/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.ImputedBgenSampleSet.modelSnps

fi;


if[ "$vers" == "v15" ]
then \
    # Match genotypes with imputed data (just on subset of snps for LMM model)
    # get list of samples in .bgen file
    # also applies to v17
    
    sampleFile=/well/ukbiobank/imputation/final/full/bgen/hrc+uk10k_sorted_8bit_rsids/imputed.sample
    cp $sampleFile chr22.sample
    cut -f1,2 -d' ' chr22.sample | tail -n +3  > toKeep.$vers.tmp
    
# get list of snps for LMM model, in rsid format.
    snpSetbim=$basedir/data/PCA/b1__b11-b001__b095-autosome-sampleqc-fastpca-White_British-highquality-pruned.bim  # just use SNPs from original PCA :) NOT ACTUALLY USED IN LREG!!
    genoData=/well/ukbiobank/expt/V2_QCed.export/data/calls_export/v2/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095  # <---- This is the release dataset v2. Actually we released v3 after excluding 33 inds + some snps on X chrom). c.f v6a which used genotype_gwas data
    
    # get list of snps in rsid format
    snpSet=$basedir/data/GWAS/otherGWAS/b1__b11-b001__b095-autosome-sampleqc-fastpca-White_British-highquality-pruned.snpsetRSIDs.txt
    
    Rscript getReleaseSNPIDs.R $snpSetbim $genoData.bim $snpSet
    
    # NOTE: 2 snps missing from genotype data (missed multi-allelics)
    $plink --bfile $genoData --extract $snpSet --keep-allele-order --keep toKeep.$vers.tmp --make-bed --out $basedir/data/GWAS/otherGWAS/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.ImputedBgenSampleSet-for-$vers.modelSnps


    ### NOTE: this file is also used for v19
fi;

if[ "$vers" == "v15s" ]
then \
    # Match genotypes with imputed data (just on subset of snps for LMM model)
    # get list of samples in .bgen file  (486757 samples)
    
    sampleFile=/well/ukbiobank/imputation/final/full/bgen/hrc+uk10k_sorted_8bit_rsids/chrX.sample
    cp $sampleFile chrX.sample
    cut -f1,2 -d' ' chrX.sample | tail -n +3  > toKeep.$vers.tmp
    
    genoData=/well/ukbiobank/expt/V2_QCed.export/data/calls_export/v3/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.chr1  # <---- This is the release dataset v2. Actually we released v3 after excluding 33 inds + some snps on X chrom). c.f v6a which used genotype_gwas data
    
    # Just pick 2 dummy snps
    printf "rs12184325\nrs116390263" > dummy.snps.txt

    $plink --bfile $genoData --extract dummy.snps.txt --keep-allele-order --keep toKeep.$vers.tmp --make-bed --out $basedir/data/GWAS/otherGWAS/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.ImputedBgenSampleSet-for-$vers.modelSnps


fi;

if[ "$vers" == "v19s" ]
then \
    # NOTE: must first run appropriate section of subset-sex-chroms.sh
    # Match genotypes with imputed data (just on subset of snps for LMM model)
    # get list of samples in .bgen file  (486757 samples)

    ####### For the X (486757 samples)
    #    sampleFile=/well/ukbiobank/imputation/final/full/bgen/hrc+uk10k_sorted_8bit_rsids/v2/chrX.sample
    sampleFile=$basedir/data/imputed/chrX.v2.sample
    cp $sampleFile chrX.sample
    cut -f1,2 -d' ' chrX.sample | tail -n +3  > toKeep.$vers.tmp
    
    genoData=/well/ukbiobank/expt/V2_QCed.export/data/calls_export/v3/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.chr1  # <---- This is the release dataset v2. Actually we released v3 after excluding 33 inds + some snps on X chrom). c.f v6a which used genotype_gwas data
    
    # Just pick 2 dummy snps on chromosome 1
    printf "rs12184325\nrs116390263" > dummy.snps.txt

    $plink --bfile $genoData --extract dummy.snps.txt --keep-allele-order --keep toKeep.$vers.tmp --make-bed --out $basedir/data/GWAS/otherGWAS/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.ImputedBgenSampleSet-for-$vers.X.modelSnps

    ####### For the PAR (486443 samples)
    sampleFile=$basedir/data/imputed/chrPAR1.v2.sample
    cp $sampleFile chrPAR.sample
    cut -f1,2 -d' ' chrPAR.sample | tail -n +3  > toKeep.$vers.tmp
    
    genoData=/well/ukbiobank/expt/V2_QCed.export/data/calls_export/v3/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.chr1  # <---- This is the release dataset v2. Actually we released v3 after excluding 33 inds + some snps on X chrom). c.f v6a which used genotype_gwas data
    
    # Just pick 2 dummy snps on chromosome 1
    printf "rs12184325\nrs116390263" > dummy.snps.txt

    $plink --bfile $genoData --extract dummy.snps.txt --keep-allele-order --keep toKeep.$vers.tmp --make-bed --out $basedir/data/GWAS/otherGWAS/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.ImputedBgenSampleSet-for-$vers.PAR.modelSnps

fi;


if[ "$vers" == "v16" ]
then \

    # which samples are not in the imputed data, but are in the genotype data?
    Rscript getSampleExclusionList.R /well/ukbiobank/expt/V2_QCed.export/data/calls_export/v3/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.chr22.fam $basedir/data/GWAS/otherGWAS/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.ImputedBgenSampleSet-for-v15.modelSnps.fam samplesToExclude.$vers.txt
    # 968 samples to exclude from the genotypes so it matches the samples in the imputed data.
fi;

if[ "$vers" == "v16b" ]
then \
    cat $basedir/QC-Scripts/GWAS/otherGWAS/samplesToExclude.v16.txt $basedir/data/GWAS/otherGWAS/SampleExclusionsForBOLT-v19b.txt > $basedir/QC-Scripts/GWAS/otherGWAS/samplesToExclude.$vers.txt
fi;

for pheno in $phenotypes;
do \

    cd $basedir/QC-Scripts/GWAS/otherGWAS/$pheno
#qsub -q short.qa -P donnelly.prja runBOLT_$vers-LREG-sex-$pheno.sh  # run on a-nodes as don't need much memory with sex chromosomes
#head -n 20 runBOLT_$vers-LREG-$pheno.sh | tail -n 10
# genoteyp data (all autosomes):
qsub -q short.qc -P donnelly.prjc runBOLT_$vers-LMM-$pheno.sh 
#for vers in $versions; do qsub -q short.qc -P donnelly.prjc runBOLT_$vers-LMM-$pheno.sh; done;

# by chromosome for imputed data: v13. Run 16-22 November 10th.
#qsub -t 16-22 -q short.qc -P donnelly.prjc runBOLT_$vers-LMM-$pheno.sh 
# for IOP and height run on coolibah without queue
#source ./runBOLT_$vers-LMM-$pheno.sh >> Logs/bolt_lmm_$vers.log 2>&1

# for checks on final dataset. Just chromosome 2.
qsub -t 1-22 -q long.qc -P donnelly.prjc runBOLT_$vers-LMM-$pheno.sh 

# by chromosome for sex chroms
qsub -t 23-26 -q short.qb -P donnelly.prjb runBOLT_$vers-LMM-$pheno.sh 

# by chromosome for sex chroms (bgen v2 in chunks) <====== TO HERE 14/9/2017
qsub -q long.qc -P donnelly.prjc runBOLT_$vers-LMM-X-$pheno.sh 
qsub -q short.qc -P donnelly.prjc runBOLT_$vers-LMM-PAR1-$pheno.sh 
qsub -q short.qc -P donnelly.prjc runBOLT_$vers-LMM-PAR2-$pheno.sh 

# for re-run with new version after setting proper annotations (v19)!
qsub -t 2 -q coolibah.q -P donnelly.prjc runBOLT_$vers-LMM-$pheno.sh
qsub -t 2 -q long.qc -P donnelly.prjc runBOLT_$vers-LMM-$pheno.sh 

# for sub-sampling (genotypes). V16b
qsub -q coolibah.q -P donnelly.prjc runBOLT_$vers-LMM-$pheno.sh

# for v20 and v20b ## <=== queued chroms 13-22 on Friday 10th August, for both v20 and v20b.  Submitted 5-12 on August 11th. Submitted 1-4 on August 12th. screen -r bolt.
qsub -t 13-22 -q long.qc -P donnelly.prjc runBOLT_$vers-LMM-$pheno.sh 
qsub -t 5-12 -q long.qc -P donnelly.prjc runBOLT_$vers-LMM-$pheno.sh 
qsub -t 1-4 -q long.qc -P donnelly.prjc runBOLT_$vers-LMM-$pheno.sh 

done;

######### FOR X CHROMOSOME BGEN (v15s or v19s). RE-COMBINE CHUNKS INTO ONE OUTPUT DATASET (makes plotting easier)


for pheno in $phenotypes2;
do \
    echo $pheno
    chr=X; outchr=23
    #chr=XY; outchr=25
    chr=PAR; outchr=25

    if [ $chr == "X" ]
    then
	ls -v $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-${vers}-chr$chr.*.out | xargs cat > tmp.txt
    fi;
    
    if [ "$chr" == "PAR" ] # COMBINE THE PARS
    then

	ls -v $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-${vers}-chrPAR1.*.out | xargs cat > tmp.txt
	ls -v $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-${vers}-chrPAR2.*.out | xargs cat >> tmp.txt

    fi;
    
    # remove the header lines and convert chromosome back to X
    cat <(head -n 1 tmp.txt) <(awk '!/\<SNP\>/ && NR > 1' tmp.txt | awk -v c=$chr 'OFS="\t" {print $1,c,$3,$4,$5,$6,$7,$8,$9,$10}')  > $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-${vers}-chr$outchr.out    
    
done;
rm tmp.txt


##################################
## Basic manhattan plots -- Genotyped data

for pheno in $phenotypes;
do \
    #for vers in $versions;
    #do \
	cd $basedir/QC-Scripts/GWAS/otherGWAS/$pheno
	
	mkdir plots

	#Rscript ../plot-BOLT-results.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LREG-sex-$vers.out 1,3,4 plots &
	#Rscript ../plot-BOLT-results.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out all plots
	#Rscript ../plot-BOLT-results.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out all plots -qc -ymax 50 -title -QCfiltered
	# no qc (all chroms separately)
	Rscript ../plot-BOLT-results-known-hits.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out all plots -ymax 50 -title -raw &
        # keep rare
	Rscript ../plot-BOLT-results-known-hits.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out all plots -maxmiss 0.05 -ymax 50 -title -QCfiltered &

	# default qc
	Rscript ../plot-BOLT-results-known-hits.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out genome plots -qc -ymax 50 -title -QCfiltered &
	
	Rscript ../plot-BOLT-results-known-hits.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out all plots -qc -ymax 50 -title -QCfiltered &
	Rscript ../plot-BOLT-results-known-hits.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out all plots -qc -ymax 50 -lreg -title -QCfiltered &

	#Rscript ../plot-BOLT-results.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out genome plots -hi -title -QCcolors
	#Rscript ../plot-BOLT-results.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out all plots -linReg -qc -ymax 50 -title -QCfiltered-LREG

	# for sex chroms
	Rscript ../plot-BOLT-results-known-hits.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers-%%.out all plots -qc -ymax 50 -title -QCfiltered -sex
	Rscript ../plot-BOLT-results-known-hits.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers-%%.out all plots -ymax 50 -title -raw -sex
	
    #done;
done;



#########
# Place of Birth particular runs

# North coordinate with all SNPs in LMM (does it make a difference to inflation??)
pheno=Place.of.birth.in.UK...north.co.ordinate
cd $basedir/QC-Scripts/GWAS/otherGWAS/$pheno
qsub -q short.qc -P donnelly.prjc runBOLT_$vers-LMM-all-snps-$pheno.sh # run on Jarrah if you can
Rscript ../plot-BOLT-results.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-all-snps-$vers.out all plots -hi -title -QCcolors

Rscript ../plot-BOLT-results.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-all-snps-$vers.out genome plots -hi -title -QCcolors

#########

###############################################################
# Plots for GWAS example with GWAS hits highlighted: Height + IOP


#pheno=Standing.height
#pheno=Intra.ocular.pressure..Goldmann.correlated..mean

cd $basedir/QC-Scripts/GWAS/otherGWAS/$pheno

# subset GWAS hits
Rscript filter-catalogue.R > Logs/filter-catalogue-$pheno.log

hits=$basedir/QC-Scripts/GWAS/otherGWAS/GWAScatalogue/hg19/gwasCatalog-subset-$pheno.RData

for pheno in $phenotypes;
do \

cd $basedir/QC-Scripts/GWAS/otherGWAS/$pheno

#for vers in $versions;
#do \

Rscript ../plot-BOLT-results-known-hits.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out all plots -ymax 50 -title -raw

Rscript ../plot-BOLT-results-known-hits.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out all plots -ymax 50 -lreg -title -raw

Rscript ../plot-BOLT-results-known-hits.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out all plots -qc -ymax 50 -title -QCfiltered
Rscript ../plot-BOLT-results-known-hits.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out all plots -qc -ymax 50 -lreg -title -QCfiltered

    Rscript ../plot-BOLT-results-known-hits.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out all plots -qc -minmaf 0 -ymax 50 -title -QCfiltered

    Rscript ../plot-BOLT-results-known-hits.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out genome plots -qc -ymax 50 -lreg -hits $hits -title -QCfiltered-Euro-hits &

    Rscript ../plot-BOLT-results-known-hits.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out genome plots -ymax 50 -lreg -hits $hits -title -raw-Euro-hits &

    Rscript ../plot-BOLT-results-known-hits.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out all plots -ymax 50 -hits $hits -title -raw-Euro-hits

#    Rscript ../plot-BOLT-results-known-hits.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out genome plots -qc -ymax 50 -title -QCfiltered

    #Rscript ../plot-BOLT-results-known-hits.R test.out all plots -ymax 50 -hits $hits -title -Euro-hits

## sex chroms
    Rscript ../plot-BOLT-results-known-hits.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers-%%.out all plots -ymax 50 -hits $hits -lreg -title -raw-Euro-hits -sex &
    Rscript ../plot-BOLT-results-known-hits.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers-%%.out all plots -qc -ymax 50 -hits $hits -lreg -title -QCfiltered-Euro-hits -sex &

## X chrom
        Rscript ../plot-BOLT-results-known-hits.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers-23.out 23 plots -ymax 50 -hits $hits -lreg -title -raw-Euro-hits &
        Rscript ../plot-BOLT-results-known-hits.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers-23.out 23 plots -ymax 50 -qc -hits $hits -lreg -title -QCFiltered-Euro-hits &

done;

##################################
## Basic manhattan plots -- Imputed data

# for bgen versions with chromosomes (chr2) --> 16-22 in November 2016.
# NOTE: plotting both LREG and LMM results now
for pheno in $phenotypes2;
do \
    cd $basedir/QC-Scripts/GWAS/otherGWAS/$pheno

    cat ../submit.header.txt <( echo Rscript ../plot-BOLT-results-known-hits.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-${vers}-chr%%.out all plots -ymax 50 -lreg -title -raw ) > $pheno.basic.manhattans.raw.sh
    cat ../submit.header.txt <( echo Rscript ../plot-BOLT-results-known-hits.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-${vers}-chr%%.out all plots -ymax 50 -title -raw ) > $pheno.basic.manhattans.raw.lmm.sh

    cat ../submit.header.txt <( echo Rscript ../plot-BOLT-results-known-hits.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-${vers}-chr%%.out all plots -mininfo 0.3 -ymax 50 -lreg -title -QCfiltered  ) > $pheno.basic.manhattans.qc1.sh
    cat ../submit.header.txt <( echo Rscript ../plot-BOLT-results-known-hits.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-${vers}-chr%%.out all plots -mininfo 0.3 -ymax 50 -title -QCfiltered  ) > $pheno.basic.manhattans.qc1.lmm.sh

    cat ../submit.header.txt <( echo Rscript ../plot-BOLT-results-known-hits.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-${vers}-chr%%.out all plots -mininfo 0.3 -minmaf 0.001 -ymax 50 -lreg -title -QCfiltered  ) > $pheno.basic.manhattans.qc2.sh
    cat ../submit.header.txt <( echo Rscript ../plot-BOLT-results-known-hits.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-${vers}-chr%%.out all plots -mininfo 0.3 -minmaf 0.001 -ymax 50 -title -QCfiltered  ) > $pheno.basic.manhattans.qc2.lmm.sh


    if [[ "$pheno" == "Standing.height" || "$pheno" == "Intra.ocular.pressure..Goldmann.correlated..mean" ]];
    then \
	hits=$basedir/QC-Scripts/GWAS/otherGWAS/GWAScatalogue/hg19/gwasCatalog-subset-$pheno.RData
	
	cat ../submit.header.txt <( echo Rscript ../plot-BOLT-results-known-hits.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-${vers}-chr%%.out all plots -ymax 50 -hits $hits -lreg -title -raw-Euro-hits ) > $pheno.basic.manhattans.raw.sh
	cat ../submit.header.txt <( echo Rscript ../plot-BOLT-results-known-hits.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-${vers}-chr%%.out all plots -ymax 50 -hits $hits -title -raw-Euro-hits ) > $pheno.basic.manhattans.raw.lmm.sh
	cat ../submit.header.txt <( echo Rscript ../plot-BOLT-results-known-hits.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-${vers}-chr%%.out all plots -mininfo 0.3 -minmaf 0.001 -ymax 50 -hits $hits -lreg -title -QCfiltered-Euro-hits ) > $pheno.basic.manhattans.qc2.sh
	cat ../submit.header.txt <( echo Rscript ../plot-BOLT-results-known-hits.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-${vers}-chr%%.out all plots -mininfo 0.3 -minmaf 0.001 -ymax 50 -hits $hits -title -QCfiltered-Euro-hits ) > $pheno.basic.manhattans.qc2.lmm.sh
	

    fi;

done;

# run in cluster
for pheno in $phenotypes2;
do \
    cd $basedir/QC-Scripts/GWAS/otherGWAS/$pheno
    qsub -N $pheno.basic.manhattans.raw -q short.qb -P donnelly.prjb $pheno.basic.manhattans.raw.sh
    qsub -N $pheno.basic.manhattans.raw.lmm -q short.qb -P donnelly.prjb $pheno.basic.manhattans.raw.lmm.sh
    qsub -N $pheno.basic.manhattans.qc1 -q short.qb -P donnelly.prjb $pheno.basic.manhattans.qc1.sh
    qsub -N $pheno.basic.manhattans.qc1.lmm -q short.qb -P donnelly.prjb $pheno.basic.manhattans.qc1.lmm.sh
    qsub -N $pheno.basic.manhattans.qc2 -q short.qb -P donnelly.prjb $pheno.basic.manhattans.qc2.sh
    qsub -N $pheno.basic.manhattans.qc2.lmm -q short.qb -P donnelly.prjb $pheno.basic.manhattans.qc2.lmm.sh
    
done;

# THE ABOVE IS RUNNING ON screen -r gwas. Might take a whole day.  2/1/2017

# Have only run plots for chrom2 and standing height for v15 and v16. 8/3/2017


######################################################
# Manhattans comparing genotypes and imputed data!!!! Chromosome manhattans only.

vers2=v13  # this is the imputed version <=== old, and some chromosomes missing lots of data!
vers=v14   # this is the genotyped version

vers2=v15  # this is the imputed version <=== given to ukbiobank in March. 
vers=v16   # this is the genotyped version <==== given to ukbiobank in March (excluded 33 individuals)

vers2=v17  # this is the imputed version --- random phenotype
vers=v18   # this is the genotyped version --- random phenotype

vers2=v15s  # this is the imputed version <=== X chromosome
vers=v16s   # this is the genotyped version <==== given to ukbiobank in March (excluded 33 individuals)

vers2=v19  # this is the imputed version <=== Following re-do of imputed data (v2)
vers=v16   # this is the genotyped version <==== given to ukbiobank in March (excluded 33 individuals)

vers2=v19s  # this is the imputed version <=== chromsome X Following re-do of imputed data (v2)
vers=v16s   # this is the genotyped version <==== given to ukbiobank in March (excluded 33 individuals)

for pheno in $phenotypes;
do \
    
   cd $basedir/QC-Scripts/GWAS/otherGWAS/$pheno

   # for bgen v2 autosomes and v1 autosomes
   bgenSample=/well/ukbiobank/imputation/final/full/bgen/hrc+uk10k_sorted_8bit_rsids/imputed.sample
   
   # for bgen v2 X chrom
   bgenSample=/well/ukbiobank/imputation/final/full/bgen/hrc+uk10k_sorted_8bit_rsids/imputed.sample/v2/chrX.v2.sample
   # for bgen v2 PAR chrom
   bgenSample=/well/ukbiobank/imputation/final/full/bgen/hrc+uk10k_sorted_8bit_rsids/imputed.sample/v2/chrPAR1.v2.sample
   
   echo Rscript ../plot-GWAS-imputation-comparison.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out \$SGE_TASK_ID plots -ymax 50 -phenoFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt -lreg -title -raw -bgenFile $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers2/$pheno-BOLT-LMM-$vers2-chr\$SGE_TASK_ID.out -sample $bgenSample -dontComputeLD -par 2 > $pheno.raw.impute.geno.plots.sh

   echo Rscript ../plot-GWAS-imputation-comparison.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out \$SGE_TASK_ID plots -ymax 50 -mininfo 0.3 -maxmiss 0.02 -phenoFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt -lreg -title -bgenFile $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers2/$pheno-BOLT-LMM-$vers2-chr\$SGE_TASK_ID.out -sample $bgenSample -dontComputeLD -par 2  > $pheno.mininfo0.3.impute.geno.plots.sh

   echo Rscript ../plot-GWAS-imputation-comparison.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out \$SGE_TASK_ID plots -ymax 50 -qc -phenoFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt -lreg -title -QCFiltered -bgenFile $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers2/$pheno-BOLT-LMM-$vers2-chr\$SGE_TASK_ID.out -sample $bgenSample -dontComputeLD -par 2  > $pheno.QCFiltered.impute.geno.plots.sh

   echo Rscript ../plot-GWAS-imputation-comparison.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out \$SGE_TASK_ID plots -ymax 50 -qc -phenoFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt -title -QCFiltered -bgenFile $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers2/$pheno-BOLT-LMM-$vers2-chr\$SGE_TASK_ID.out -sample $bgenSample -dontComputeLD -par 2  > $pheno.QCFiltered.lmm.impute.geno.plots.sh

   # plot hits on genotype data for Standingheight and IOP
   #echo Rscript ../plot-GWAS-imputation-comparison.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out \$SGE_TASK_ID plots -ymax 50 -phenoFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt -lreg -title -raw -bgenFile $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers2/$pheno-BOLT-LMM-$vers2-chr\$SGE_TASK_ID.out -sample /well/ukbiobank/imputation/final/full/bgen/hrc+uk10k_sorted_v1.1/chr\$SGE_TASK_ID.hrc+uk10k.I4.v1.1.sample -dontComputeLD -hits $hits > $pheno.raw.hits.impute.geno.plots.sh
   
done;

for pheno in $phenotypes;
do \
    cd $basedir/QC-Scripts/GWAS/otherGWAS/$pheno
    
    qsub -t 25 -N $pheno.$vers.raw.impute.geno.plots -pe shmem 2 -o Logs -j y -V -cwd -q short.qb -P donnelly.prjb $pheno.raw.impute.geno.plots.sh
    qsub -t 25 -N $pheno.$vers.mininfo0.3.impute.geno.plots -pe shmem 2 -o Logs -j y -V -cwd -q short.qb -P donnelly.prjb $pheno.mininfo0.3.impute.geno.plots.sh
    qsub -t 25 -N $pheno.$vers.QCFiltered.impute.geno.plots -pe shmem 2 -o Logs -j y -V -cwd -q short.qb -P donnelly.prjb $pheno.QCFiltered.impute.geno.plots.sh
    qsub -t 25 -N $pheno.$vers.QCFiltered.impute.geno.plots -pe shmem 2 -o Logs -j y -V -cwd -q short.qb -P donnelly.prjb $pheno.QCFiltered.lmm.impute.geno.plots.sh
    #qsub -t 23 -N $pheno.$vers.raw.hits.impute.geno.plots -pe shmem 2 -o Logs -j y -V -cwd -q short.qb -P donnelly.prjb $pheno.raw.hits.impute.geno.plots.sh

done;

# THE ABOVE IS RUNNING ON CLUSTER 8/3/2017 <== Done.

# Get p-value orders for random phenotypes
for pheno in $phenotypes;
do \

cd $basedir/QC-Scripts/GWAS/otherGWAS/$pheno    
Rscript ../pvalue-sorted.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers-chr%%.out > Logs/pvalue-sorted.$vers.log &

done;


######################################################
# Compare results from different versions (x/y plots)

for pheno in $phenotypes;
do \
    cd $basedir/QC-Scripts/GWAS/otherGWAS/$pheno
    i=2 # which chromosome?

    # genotypes
   # v1=v14;    v2=v16;
   # Rscript ../plot-GWAS-compare.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$v1/$pheno-BOLT-LMM-$v1.out $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$v2/$pheno-BOLT-LMM-$v2.out $i plots > Logs/plot-GWAS-compare-genotyped-$v1.$v2.log &
    
   # Rscript ../plot-GWAS-compare.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$v2/$pheno-BOLT-LMM-$v2.out $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$v2/$pheno-BOLT-LMM-$v2.out $i plots -lreg1 > Logs/plot-GWAS-compare-genotyped-$v1.lreg.$v2.log &  # v2 lreg vs lmm
    
    # imputed
#    v1=v13;    v2=v15;
#    Rscript ../plot-GWAS-compare.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$v1/$pheno-BOLT-LMM-$v1-chr$i.out $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$v2/$pheno-BOLT-LMM-$v2-chr$i.out $i plots -lreg Logs/plot-GWAS-compare-imputed-lreg.$v1.$v2.log &
    v1=v16
    Rscript ../plot-GWAS-compare.R $basedir/data/GWAS/otherGWAS/$pheno.qnorm/BOLTLMM.$v1/$pheno.qnorm-BOLT-LMM-$v1.out $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$v1/$pheno-BOLT-LMM-$v1.out $i plots -lreg Logs/plot-GWAS-compare-lreg.$v1.qnorm.log &
    v2=v15
    Rscript ../plot-GWAS-compare.R $basedir/data/GWAS/otherGWAS/$pheno.qnorm/BOLTLMM.$v2/$pheno.qnorm-BOLT-LMM-$v2-chr$i.out $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$v2/$pheno-BOLT-LMM-$v2-chr$i.out $i plots -lreg Logs/plot-GWAS-compare-lreg.$v2.qnorm.log &

     
done;
#v1=v15
#v2=v15
#    Rscript ../plot-GWAS-compare.R $basedir/data/GWAS/otherGWAS/dQC/BOLTLMM.$v2/dQC-BOLT-LMM-$v2-chr$i.out $basedir/data/GWAS/otherGWAS/Cluster.CR/BOLTLMM.$v2/Cluster.CR-BOLT-LMM-$v2-chr$i.out $i plots -lreg Logs/plot-GWAS-compare-lreg.$v2.log &

    

######################################################
# Specific region plots for genotype/imputation comparison

# once you have run get.hit.regions.R you can use these results to plot the regions separately, based on the set of regions in regions.to.plot.txt
# NOTE: will get hit regions based on the genotype data (i.e v14) rather than the imputed data!


pheno=Array.as.binary
pheno=Standing.height
pheno=Place.of.birth.in.UK...north.co.ordinate
pheno=Internal.Pico..ng.uL.
pheno=Intra.ocular.pressure..Goldmann.correlated..mean
pheno=Intra.ocular.pressure..Goldmann.correlated..mean.qnorm

cd $basedir/QC-Scripts/GWAS/otherGWAS/$pheno


vers=v16
vers2=v19
vers=v16s
vers2=v19s

if [ "$pheno" == "Intra.ocular.pressure..Goldmann.correlated..mean" ];   
then \
    
    hits=$basedir/QC-Scripts/GWAS/otherGWAS/GWAScatalogue/hg19/gwasCatalog-subset-$pheno.RData
    
    Rscript ../Internal.Pico..ng.uL./plot-known-hits.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out all plots -ymax 50 -hits $hits -title -raw-hits & # genotypes
    Rscript ../Internal.Pico..ng.uL./plot-known-hits.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out genome plots -ymax 50 -hits $hits -title -raw-hits & # genot    
    Rscript ../Internal.Pico..ng.uL./plot-known-hits.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers2/$pheno-BOLT-LMM-$vers2-chr%%.out all plots -ymax 50 -hits $hits -title -raw-hits & # imputed

fi;
if [ "$pheno" == "Intra.ocular.pressure..Goldmann.correlated..mean.qnorm" ];   
then \
    
    hits=$basedir/QC-Scripts/GWAS/otherGWAS/GWAScatalogue/hg19/gwasCatalog-subset-Intra.ocular.pressure..Goldmann.correlated..mean.RData
    
    Rscript ../Internal.Pico..ng.uL./plot-known-hits.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out all plots -ymax 50 -hits $hits -title -raw-hits & # genotypes
    Rscript ../Internal.Pico..ng.uL./plot-known-hits.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out genome plots -ymax 50 -hits $hits -title -raw-hits & # genot     
    Rscript ../Internal.Pico..ng.uL./plot-known-hits.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers2/$pheno-BOLT-LMM-$vers2-chr%%.out all plots -ymax 50 -hits $hits -title -raw-hits & # imputed

fi;

# extract specific regions for Array.as.binary (create regions dataset manually based on looking at manhattan plots!):
if [ "$pheno" == "Array.as.binary" ];
then \
        
    Rscript ../get.hit.regions.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out 1-22 plots -ymax 50 -phenoFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt -lreg -title -raw -regions $basedir/QC-Scripts/GWAS/otherGWAS/$pheno/regions.to.plot.$vers.txt &
    
    echo Rscript ../plot-GWAS-imputation-comparison.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out \$SGE_TASK_ID plots -ymax 50 -phenoFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt -lreg -title -raw -bgenFile $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers2/$pheno-BOLT-LMM-$vers2-chr\$SGE_TASK_ID.out -sample /well/ukbiobank/imputation/final/full/bgen/chr\$SGE_TASK_ID.hrc+uk10k.I4.v1.1.sample -dontPlotManhattans -ldRData $pheno-BOLT-LMM-$vers.out.chrgenome.maf0.miss1.pruned-raw-lreg-hitsAndLDcalc.RData > $pheno.$vers.with.$vers2.raw.impute.geno.plots.regions.sh
    
fi;

if [ "$pheno" == "Standing.height" ];
then \
        hits=$basedir/QC-Scripts/GWAS/otherGWAS/GWAScatalogue/hg19/gwasCatalog-subset-$pheno.RData

    # extract regions for standing height based on (default) qc-filtered genotype data! BUT: plot all data around hit region (i.e no qc)
    # NOTE: Females only used in ld for one of the versions (has female in the name) 
    Rscript ../region.plots.new.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out 23 plots -ymax 50 -qc -phenoFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt -lreg -title -QCFiltered -dontPlotRegions & # <===== LREG regions

    Rscript ../region.plots.new.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out 1-22 plots -ymax 50 -qc -phenoFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt -title -QCFiltered -dontPlotRegions & # <===== LMM regions (autosomes only)

    # Get imputation LD stats for chromosome 2 only, based on genotype hit regions. But first run $basedir/data/imputed/bgenixIndex.sh. TAKES A WHILE!
    i=23; i2=X;# callT=0.3   # <=== Which chromosome to do?
    i=2; i2=2
    i=25; i2=PAR

    # first is for v13, which used indexed files from bgen v1.1
    #Rscript ../get.hit.regions.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers2/$pheno-BOLT-LMM-$vers2-chr$i.out $i plots -ymax 50 -phenoFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt -lreg -title -raw -bgen /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/imputed/chr$i.hrc+uk10k.I4.v1.1.bgen -sample /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/imputed/chr$i.hrc+uk10k.I4.v1.1.sample -regions $pheno-BOLT-LMM-$vers.out.chrgenome.maf0.001.miss0.05.pruned-QCFiltered-lreg-hitsAndLDcalc.RData

    # Second is for v15, which used indexed files from bgen v1.2 that colin converted to v1.1    
    # Compute LD around top hit snps for imputed data, based on regions from QC'd genotype data.

    #bgenFile=/well/ukbiobank/expt/V2_QCed.imputation.sanity_check/data/chr$i2.hrc+uk10k_sorted_8bit_rsids.v1.1.bgen
    # Make a sym-link to the real data
    ln -s /well/ukbiobank/expt/full_release_issues/uk10k_imputation_annotation/hrc+uk10k_sorted_8bit_rsids_chr$i2.v4.v1_1.bgen $basedir/data/imputed/hrc+uk10k_sorted_8bit_rsids_chr$i2.v4.v1_1.bgen

    bgenFile=$basedir/data/imputed/hrc+uk10k_sorted_8bit_rsids_chr$i2.v4.v1_1.bgen
    bgenSample=/well/ukbiobank/imputation/final/full/bgen/hrc+uk10k_sorted_8bit_rsids/imputed.sample
    # for X chromosome
    bgenSample=$basedir/data/imputed/chrX.v2.sample
    # for PAR
    bgenSample=$basedir/data/imputed/chr$i2.v2.sample # same for both pars

     
    ### Using a different callthreshold (only necessary in X chrom, maybe)
    #Rscript ../get.hit.regions.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers2/$pheno-BOLT-LMM-$vers2-chr$i.out $i plots -ymax 50 -phenoFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt -lreg -title -raw -bgen $bgenFile -sample /well/ukbiobank/imputation/final/full/bgen/hrc+uk10k_sorted_8bit_rsids/imputed.sample -regions $pheno-BOLT-LMM-$vers-$i.out.chr$i.maf0.001.miss0.05.pruned-QCFiltered-lreg-hitsAndLDcalc.RData -callthreshold $callT > Logs/get.hit.regions-$vers2-lreg-raw.log &

    # NOTE: top snp is found with qcfiltered imputed data, but ld computed with all snps!
    # autosomes
    theRegions=$pheno-BOLT-LMM-$vers.out.chrgenome.maf0.001.miss0.05.pruned-QCFiltered-lreg-hitsAndLDcalc.RData
    # X chrom
    theRegions=$pheno-BOLT-LMM-$vers-$i.out.chr$i.maf0.001.miss0.05.pruned-QCFiltered-lreg-hitsAndLDcalc.RData
    
    Rscript ../get.hit.regions.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers2/$pheno-BOLT-LMM-$vers2-chr$i.out $i plots -ymax 50 -phenoFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt -qc -lreg -title -QCFiltered -bgen $bgenFile -sample $bgenSample -regions $theRegions -callthreshold 0.1 > Logs/get.hit.regions-$vers2-lreg-QCFiltered.log &

	# The above two running on coolibah. 10/5/2017
	
    rm topVariants.*.tmp
    #rm genotypes*tmp*
    
    # actually do the plotting. bgen v1
    echo Rscript ../plot-GWAS-imputation-comparison.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out $i plots -ymax 50 -phenoFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt -lreg -title -raw.$vers2 -bgenFile $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers2/$pheno-BOLT-LMM-$vers2-chr$i.out -sample $bgenSample -dontPlotManhattans -ldRData $pheno-BOLT-LMM-$vers.out.chrgenome.maf0.001.miss0.05.pruned-QCFiltered-lreg-hitsAndLDcalc.RData -ldRDataImp $pheno-BOLT-LMM-$vers2-chr$i.out.chr$i2.maf0.info0.pruned-raw-lreg-hitsAndLDcalc.RData > $pheno.$vers.with.$vers2.raw.impute.geno.plots.regions.sh

    # Plot special region where GIANT has a hit... ===> do this in plot-GWAS-imputation-comparison.R
    # Chromosome X === raw bgen (ct 0.3)
    echo Rscript ../plot-GWAS-imputation-comparison.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out $i plots -ymax 50 -phenoFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt -lreg -title -raw -bgenFile $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers2/$pheno-BOLT-LMM-$vers2-chr$i.out -sample /well/ukbiobank/imputation/final/full/bgen/chrX.hrc+uk10k.I4.v1.1.sample -dontPlotManhattans -ldRData $pheno-BOLT-LMM-$vers-$i.out.chr$i2.maf0.001.miss0.05.pruned-QCFiltered-lreg-hitsAndLDcalc.RData -ldRDataImp $pheno-BOLT-LMM-$vers2-chr$i.out.chr$i2.maf0.info0.pruned-raw-lreg-ct0.1-hitsAndLDcalc.RData -hits $hits > $pheno.$vers.with.$vers2.raw.impute.geno.plots.regions.sh

        # Chromosome X === raw bgen (ct 0.1)
    echo Rscript ../plot-GWAS-imputation-comparison.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out 23 plots -ymax 50 -phenoFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt -lreg -title -rawCt0.1 -bgenFile $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers2/$pheno-BOLT-LMM-$vers2-chr23.out -sample $bgenSample -dontPlotManhattans -ldRData $pheno-BOLT-LMM-$vers-23.out.chr23.maf0.001.miss0.05.pruned-QCFiltered-lreg-hitsAndLDcalc.RData -ldRDataImp $pheno-BOLT-LMM-$vers2-chr23.out.chrX.maf0.info0.pruned-raw-lreg-ct0.1-hitsAndLDcalc.RData -hits $hits > $pheno.$vers.with.$vers2.rawct0.1.impute.geno.plots.regions.sh
    
        # Chromosome X and PAR === qcd bgen v2
    echo Rscript ../plot-GWAS-imputation-comparison.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out $i plots -ymax 50 -phenoFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt -lreg -title -rawQChits.$vers2 -bgenFile $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers2/$pheno-BOLT-LMM-$vers2-chr$i.out -sample $bgenSample -ldRData $pheno-BOLT-LMM-$vers-$i.out.chr$i.maf0.001.miss0.05.pruned-QCFiltered-lreg-hitsAndLDcalc.RData -ldRDataImp $pheno-BOLT-LMM-$vers2-chr$i.out.chr$i2.maf0.001.info0.3.pruned-QCFiltered-lreg-ct0.1-hitsAndLDcalc.RData -hits $hits -dontPlotManhattans > $pheno.$vers.with.$vers2.QCFiltered.impute.geno.plots.regions.sh
    
            # Autosome === qcd bgen v2
    echo Rscript ../plot-GWAS-imputation-comparison.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out $i plots -ymax 50 -phenoFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt -lreg -title -rawQChits-$vers2 -bgenFile $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers2/$pheno-BOLT-LMM-$vers2-chr$i.out -sample $bgenSample -ldRData $pheno-BOLT-LMM-$vers.out.chrgenome.maf0.001.miss0.05.pruned-QCFiltered-lreg-hitsAndLDcalc.RData -ldRDataImp $pheno-BOLT-LMM-$vers2-chr$i.out.chr$i2.maf0.001.info0.3.pruned-QCFiltered-lreg-ct0.1-hitsAndLDcalc.RData -hits $hits > $pheno.$vers.with.$vers2.QCFiltered.impute.geno.plots.regions.sh

    source ./$pheno.$vers.with.$vers2.raw.impute.geno.plots.regions.sh > $pheno.$vers.with.$vers2.raw.impute.geno.plots.regions.log & # used for v1 bgen
    source ./$pheno.$vers.with.$vers2.QCFiltered.impute.geno.plots.regions.sh > Logs/$pheno.$vers.with.$vers2.QCFiltered.impute.geno.plots.regions.log &  # <==== just to this one for v2 bgen and X
    source ./$pheno.$vers.with.$vers2.rawct0.1.impute.geno.plots.regions.sh > pheno.$vers.with.$vers2.rawct0.1.impute.geno.plots.regions.log &

# Then run combine-region-plots.R manually
    
if [ "$pheno" == "Internal.Pico..ng.uL..qnorm" ];
then \

    # get hit regions based raw data for qnorm
    cd $basedir/QC-Scripts/GWAS/otherGWAS/$pheno

    # genotyped data hit regions and LD
    Rscript ../region.plots.new.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out all plots -ymax 50 -phenoFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt -lreg -title -raw -dontPlotRegions &

    # imputed data LD (based on genotyped regions)
    for i in `seq 1 22`;
    do \
	Rscript ../get.hit.regions.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers2/$pheno-BOLT-LMM-$vers2-chr$i.out $i plots -ymax 50 -phenoFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt -lreg -title -raw -bgen /well/ukbiobank/expt/V2_QCed.imputation.sanity_check/data/chr$i.hrc+uk10k_sorted_8bit_rsids.v1.1.bgen -sample /well/ukbiobank/imputation/final/full/bgen/hrc+uk10k_sorted_8bit_rsids/imputed.sample -regions $pheno-BOLT-LMM-$vers.out.chrgenome.maf0.miss1.pruned-raw-lreg-hitsAndLDcalc.RData > Logs/get.hit.regions-$vers2-lreg-raw-chr$i.log &
    done;     

    echo Rscript ../plot-GWAS-imputation-comparison.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out \$SGE_TASK_ID plots -ymax 50 -phenoFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt -lreg -title -raw -bgenFile $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers2/$pheno-BOLT-LMM-$vers2-chr\$SGE_TASK_ID.out -sample /well/ukbiobank/imputation/final/full/bgen/chr\$SGE_TASK_ID.hrc+uk10k.I4.v1.1.sample -dontPlotManhattans -ldRData $pheno-BOLT-LMM-$vers.out.chrgenome.maf0.miss1.pruned-raw-lreg-hitsAndLDcalc.RData -ldRDataImp $pheno-BOLT-LMM-$vers2-chr\$SGE_TASK_ID.out.chr\$SGE_TASK_ID.maf0.info0.pruned-raw-lreg-hitsAndLDcalc.RData > $pheno.$vers.with.$vers2.raw.impute.geno.plots.regions.sh

    # plot chromosome manhattans with catalogue hits.
    hits=$basedir/QC-Scripts/GWAS/otherGWAS/GWAScatalogue/hg19/gwasCatalog-subset-Internal.Pico..ng.uL..RData

    Rscript ../Internal.Pico..ng.uL./plot-known-hits.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out all plots -ymax 50 -hits $hits -title -raw-hits & # genotypes
    Rscript ../Internal.Pico..ng.uL./plot-known-hits.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers2/$pheno-BOLT-LMM-$vers2-chr%%.out all plots -ymax 50 -hits $hits -title -raw-hits & # imputed


fi;
    
if [ "$pheno" == "Internal.Pico..ng.uL." ];
then \
       
    cd $basedir/QC-Scripts/GWAS/otherGWAS/$pheno
    
    # now get regions based on QC'd data from non-qnorm version, but use regions from the qnorm version.
    Rscript ../get.hit.regions.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out 1-22 plots -ymax 50 -qc -phenoFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt -lreg -title -QCFiltered -regions ../$pheno.qnorm/$pheno.qnorm-BOLT-LMM-$vers.out.chrgenome.maf0.miss1.pruned-raw-lreg-hitsAndLDcalc.RData &

    # imputed data LD (based on genotyped regions); again use qc'd data (otherwise rare variants dominate)
    for i in `seq 1 22`;
    do \
	Rscript ../get.hit.regions.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers2/$pheno-BOLT-LMM-$vers2-chr$i.out $i plots -ymax 50 -phenoFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt -qc -lreg -title -QCFiltered -bgen /well/ukbiobank/expt/V2_QCed.imputation.sanity_check/data/chr$i.hrc+uk10k_sorted_8bit_rsids.v1.1.bgen -sample /well/ukbiobank/imputation/final/full/bgen/hrc+uk10k_sorted_8bit_rsids/imputed.sample -regions $pheno-BOLT-LMM-$vers.out.chrgenome.maf0.001.miss0.05.pruned-QCFiltered-lreg-hitsAndLDcalc.RData > Logs/get.hit.regions-$vers2-lreg-qc-chr$i.log &
    done;
    # THE ABOVE IS RUNNING.

    # plot raw data but only at qc;d regions.
    echo Rscript ../plot-GWAS-imputation-comparison.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out \$SGE_TASK_ID plots -ymax 50 -phenoFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt -qc -lreg -title -QCFiltered -bgenFile $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers2/$pheno-BOLT-LMM-$vers2-chr\$SGE_TASK_ID.out -sample /well/ukbiobank/imputation/final/full/bgen/chr\$SGE_TASK_ID.hrc+uk10k.I4.v1.1.sample -dontPlotManhattans -ldRData $pheno-BOLT-LMM-$vers.out.chrgenome.maf0.001.miss0.05.pruned-QCFiltered-lreg-hitsAndLDcalc.RData -ldRDataImp $pheno-BOLT-LMM-$vers2-chr\$SGE_TASK_ID.out.chr\$SGE_TASK_ID.maf0.001.info0.3.pruned-QCFiltered-lreg-hitsAndLDcalc.RData > $pheno.$vers.with.$vers2.qc.impute.geno.plots.regions.sh    
        # TO HERE, but doesn't work well with raw data as ld input because they aren't in LD with much!

    # plot chromosome manhattans with catalogue hits.
    Rscript plot-known-hits.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out all plots -ymax 50 -hits $hits -title -raw-hits & # genotypes
    Rscript plot-known-hits.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers2/$pheno-BOLT-LMM-$vers2-chr%%.out all plots -ymax 50 -hits $hits -title -raw-hits & # imputed

fi;



if [ "$pheno" == "Place.of.birth.in.UK...north.co.ordinate" ];
then \
    cd $basedir/QC-Scripts/GWAS/otherGWAS/$pheno
    # find regions using qc-filgered genotype data
    Rscript ../region.plots.new.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out all plots -ymax 50 -qc -phenoFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt -lreg -title -QCFiltered -dontPlotRegions &

    
    # imputed data LD (based on genotyped regions). Better to use qc'd version of data for this.
    for i in `seq 1 22`;
    do \
	Rscript ../get.hit.regions.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers2/$pheno-BOLT-LMM-$vers2-chr$i.out $i plots -ymax 50 -phenoFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt -qc -lreg -title -QCFiltered -bgen /well/ukbiobank/expt/V2_QCed.imputation.sanity_check/data/chr$i.hrc+uk10k_sorted_8bit_rsids.v1.1.bgen -sample /well/ukbiobank/imputation/final/full/bgen/hrc+uk10k_sorted_8bit_rsids/imputed.sample -regions $pheno-BOLT-LMM-$vers.out.chrgenome.maf0.001.miss0.05.pruned-QCFiltered-lreg-hitsAndLDcalc.RData > Logs/get.hit.regions-$vers2-lreg-qc-chr$i.log &
    done;     

    
    # plot them (plot raw, but only ld with qc'd snps)
    echo Rscript ../plot-GWAS-imputation-comparison.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out \$SGE_TASK_ID plots -ymax 50 -phenoFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt -qc -lreg -title -QCFiltered -bgenFile $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers2/$pheno-BOLT-LMM-$vers2-chr\$SGE_TASK_ID.out -dontPlotManhattans -ldRData $pheno-BOLT-LMM-$vers.out.chrgenome.maf0.001.miss0.05.pruned-QCFiltered-lreg-hitsAndLDcalc.RData -ldRDataImp $pheno-BOLT-LMM-$vers2-chr\$SGE_TASK_ID.out.chr\$SGE_TASK_ID.maf0.001.info0.3.pruned-QCFiltered-lreg-hitsAndLDcalc.RData > $pheno.$vers.with.$vers2.qc.impute.geno.plots.regions.sh

fi;
    
qsub -t 1-22 -N $pheno.$vers.with.$vers2.raw.impute.geno.plots.regions -pe shmem 2 -o Logs -j y -V -cwd -q short.qb -P donnelly.prjb $pheno.$vers.with.$vers2.raw.impute.geno.plots.regions.sh
qsub -t 1-22 -N $pheno.$vers.with.$vers2.qc.impute.geno.plots.regions -pe shmem 2 -o Logs -j y -V -cwd -q short.qb -P donnelly.prjb $pheno.$vers.with.$vers2.qc.impute.geno.plots.regions.sh






####################################
# Basic region plots for standing height and IOP (just the current version plotted - i.e imputation not on top etc.)

vers=v5
pheno=Standing.height
pheno=Intra.ocular.pressure..Goldmann.correlated..mean

hits=$basedir/QC-Scripts/GWAS/otherGWAS/GWAScatalogue/hg19/gwasCatalog-subset-$pheno.RData
cd $basedir/QC-Scripts/GWAS/otherGWAS/$pheno


####### Find hit regions and plot them. Other scripts use the output from this!
# use -dontPlotRegions do not plot each region separately
# use -dontComputeLD if don't want to compute LD

# autosomes
Rscript ../region.plots.new.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out genome plots -ymax 50 -qc -hits $hits -phenoFile $basedir/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt -lreg -title -recombRegions -dontComputeLD > $basedir/QC-Scripts/GWAS/otherGWAS/$pheno/$pheno.qc.region.plots.new.sh



for pheno in $phenotypes;
do \
    echo $pheno
    cd $basedir/QC-Scripts/GWAS/otherGWAS/$pheno

    # X chrom
    Rscript ../region.plots.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers-23.out 23 plots -ymax 50 -qc -hits $hits -phenoFile $basedir/data/GWAS/otherGWAS/PhenotypesForBOLT-v11-23.txt -title -recombRegions -sex

    Rscript ../region.plots.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers-23.out 23 plots -ymax 50 -hits $hits -phenoFile $basedir/data/GWAS/otherGWAS/PhenotypesForBOLT-v11-23.txt -title -recombRegions-raw -sex

    Rscript ../region.plots.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers-25.out 25 plots -ymax 50 -qc -hits $hits -phenoFile $basedir/data/GWAS/otherGWAS/PhenotypesForBOLT-v11-25.txt -title -recombRegions -sex

    # other non-autosomes chroms
    Rscript ../region.plots.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers-24.out 24 plots -ymax 50 -qc -hits $hits -phenoFile $basedir/data/GWAS/otherGWAS/PhenotypesForBOLT-v11-24.txt -title -recombRegions -sex
    Rscript ../region.plots.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers-26.out 26 plots -ymax 50 -qc -hits $hits -phenoFile $basedir/data/GWAS/otherGWAS/PhenotypesForBOLT-v11-26.txt -title -recombRegions -sex

done;


# for bgen versions
chr=2 
# first apply no qc
Rscript ../region.plots.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-${vers}-chr$chr.out genome plots -ymax 50 -hits $hits -phenoFile $basedir/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt -title -raw-recombRegions -dontComputeLD

# apply standard qc
Rscript ../region.plots.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-${vers}-chr$chr.out genome plots -ymax 50 -qc -hits $hits -phenoFile $basedir/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt -title -recombRegions -dontComputeLD



################################################################
# Conditional analyses
################################################################


########
vers=v3  # have only done conditional stuff for v3 on height.
########
vers=v16  # have only done conditional stuff for genotypes on height.
pheno=Standing.height
########

cd $basedir/QC-Scripts/GWAS/otherGWAS/$pheno

# set up for conditional analysis (have to run region.plots.R first)
# v3
Rscript ../conditional-analysis-setup.R $basedir/data/Combined/b1__b11-b001__b095-autosome-oxfordqc -vers $vers -snps $basedir/QC-Scripts/GWAS/otherGWAS/$pheno/$pheno-BOLT-LMM-$vers.out.chrgenome.maf0.001.miss0.05.pruned-recombRegions-hitsAndLDcalc.RData
# v16
Rscript ../conditional-analysis-setup.R /well/ukbiobank/expt/V2_QCed.export/data/calls_export/v3/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.chr%% -vers $vers -snps $basedir/QC-Scripts/GWAS/otherGWAS/$pheno/$pheno-BOLT-LMM-$vers.out.chrgenome.maf0.001.miss0.05.pruned-QCFiltered-lreg-hitsAndLDcalc.RData


# run conditional analysis
cp $basedir/QC-Scripts/GWAS/otherGWAS/runBOLT_$vers-LMM-conditional.sh $basedir/QC-Scripts/GWAS/otherGWAS/$pheno
sed "s/phenotypeForThisFile/$pheno/g" $basedir/QC-Scripts/GWAS/otherGWAS/$pheno/runBOLT_$vers-LMM-conditional.sh  | sed "s/versionNumber/$vers/g" > $basedir/QC-Scripts/GWAS/otherGWAS/$pheno/runBOLT_$vers-LMM-conditional-$pheno.sh

#nRegions=`wc -l $pheno-BOLT-LMM-$vers.out.chrgenome.maf0.001.miss0.05.pruned-recombRegions-snpRegionList.txt | cut -f1 -d' '`
nRegions=`wc -l $pheno-BOLT-LMM-$vers.out.chrgenome.maf0.001.miss0.05.pruned-QCFiltered-lreg-snpRegionList.txt | cut -f1 -d' '`
qsub -t 1-$nRegions -P donnelly.prjc -q short.qc $basedir/QC-Scripts/GWAS/otherGWAS/$pheno/runBOLT_$vers-LMM-conditional-$pheno.sh


##### Conditonal analysis for v16 started but aborted by clare.... need to re-run jobs.


# plot conditional analysis
Rscript ../plot-conditional.R $pheno-BOLT-LMM-$vers.out.chrgenome.maf0.001.miss0.05.pruned-recombRegions-snpRegionList.txt -readRaw -hits $hits -chr all -qc -ymax 50 -resultsPrefix $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/Standing.height-BOLT-LMM-$vers-cond -phenoFile $basedir/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt

##### plot certain snps (autosomes)
# IOP
Rscript ../cluster-region-plots.R $pheno -vers $vers -outdir clusterplots_snpsToCheck -snps Affx-27465929,Affx-37324711,Affx-3476614,Affx-12307198
# Affx-27465299 not a snp. Affx-27465929 is. chrom 6

# Height
Rscript ../cluster-region-plots.R $pheno -vers $vers -outdir clusterplots_snpsToCheck -snps Affx-22403083,Affx-25620729,Affx-50129154,Affx-34427587,Affx-34429175,Affx-16615086,Affx-16640493

# Height: SNP which is in GWAS catalogue but low in UKBiobank
Rscript ../cluster-region-plots.R $pheno -vers $vers -outdir clusterplots_snpsToCheck -snps Affx-17862683

# HLA 
Rscript ../cluster-region-plots.R $pheno -vers $vers -outdir clusterplots_snpsToCheck -snps Affx-28469679,Affx-28501022,Affx-28434783

# Affx-501291154 is not a snp.

# Array (v14) hits on chrom19 in genotype data
###########
Rscript ../cluster-region-plots.R $pheno -vers $vers -rsids -outdir clusterplots_snpsToCheck -snps Affx-89022245,Affx-79516215



##### plot certain snps (sex chroms)

# North coordinate (Y chrom + MT + PAR)
pheno=Place.of.birth.in.UK...north.co.ordinate
Rscript ../cluster-region-plots.R $pheno -vers $vers -outdir clusterplots_snpsToCheck -snps Affx-89022310,Affx-89025669,Affx-34738499

# red blood cell counts (X chrom)
pheno=Red.Blood.Cell.Count
Rscript ../cluster-region-plots.R $pheno -vers $vers -outdir clusterplots_snpsToCheck -snps $pheno-BOLT-LMM-v11-23.out.chr23.maf0.001.miss0.05.pruned-recombRegions-topVariants.txt

# Haemoglobin.Concentration
pheno=Haemoglobin.Concentration
Rscript ../cluster-region-plots.R $pheno -vers $vers -outdir clusterplots_snpsToCheck -snps $pheno-BOLT-LMM-v11-23.out.chr23.maf0.001.miss0.05.pruned-recombRegions-topVariants.txt

# Body.mass.index..BMI.
pheno=Body.mass.index..BMI.
Rscript ../cluster-region-plots.R $pheno -vers $vers -outdir clusterplots_snpsToCheck -snps $pheno-BOLT-LMM-v11-23.out.chr23.maf0.001.miss0.05.pruned-recombRegions-topVariants.txt

# standing height
pheno=Standing.height
Rscript ../cluster-region-plots.R $pheno -vers $vers -outdir clusterplots_snpsToCheck -snps Affx-34918238,Affx-34868823,Affx-52308892,Affx-37660483,Affx-35022227,Affx-37665723
Rscript ../cluster-region-plots.R $pheno -vers $vers -outdir clusterplots_snpsToCheck -snps $pheno-BOLT-LMM-v11-25.out.chr25.maf0.001.miss0.05.pruned-recombRegions-topVariants.txt

# array (X chrom)
pheno=Array.as.binary
# first region.plots by snp (for ones that aren't in the 'top hits' list)
Rscript ../region.plots.by.snp.R $pheno-BOLT-LMM-$vers-23.out.chr23.maf0.miss1.pruned-recombRegions-raw-hitsAndLDcalc.RData allsnps 23 -range 1.53e8,1.534e8 # print all snps within range of weird thing.
Rscript ../cluster-region-plots.R $pheno -vers $vers -outdir clusterplots_snpsToCheck -snps Affx-89019147,Affx-89011845,Affx-89011807,Affx-89023616,Affx-89023210,Affx-80267005



######### testing

#Rscript $basedir/QC-Scripts/GWAS/otherGWAS/plot-BOLT-results-known-hits.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out 22 plots -ymax 50 -hits $hits -title -raw-Euro-hits


#pheno=Intra.ocular.pressure..Goldmann.correlated..mean
#hits=$basedir/QC-Scripts/GWAS/otherGWAS/GWAScatalogue/hg19/gwasCatalog-subset-$pheno.RData

#Rscript $basedir/QC-Scripts/GWAS/otherGWAS/plot-BOLT-results-known-hits.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out 2 plots -ymax 50 -hits $hits -title -raw-Euro-hits


#Rscript $basedir/QC-Scripts/GWAS/otherGWAS/region.plots.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out genome plots -ymax 50 -qc -title -test -dontPlotRegions -dontComputeLD
