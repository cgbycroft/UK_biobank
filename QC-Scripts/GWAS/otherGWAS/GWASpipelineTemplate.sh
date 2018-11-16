########################
# Template for GWAS pipeline. Don't run this script itself, but run steps separately by copy/pasting into commandline
########################
# Authorship: Clare Bycroft; July 2016


##############################
# 0. Initial bash variables

# Latest version of plink
plink=/well/ukbiobank/qcoutput/Src/plink-v1.90b3k-x86_64/plink

# This is where data and scripts for UKBiobank are stored.
basedir=/well/ukbiobank/qcoutput.V2_QCed.sample-QC


##############################
# 1. Set up version parameters

#######
vers=v7
phenotypes=`echo Intra.ocular.pressure..Goldmann.correlated..mean Standing.height`
phenotypeFile=$basedir/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt # v3 is file with white british samples, unrelated samples only
####### 


##############################
# 2. Make sure that the phenotypes are actually in the phenotypeFile.
for pheno in $phenotypes;
do \
grep -w $pheno <(head -n 1 $phenotypeFile) | wc -l
# values should all be 1.
done;


##############################
# 3. MANUAL STEP:  Write template file with desired bolt input data and input parameters. Save in $basedir/QC-Scripts/GWAS/otherGWAS

# Call this file runBOLT_$vers-LMM.sh so that the next step works
# Example using genotype data: runBOLT_v5-LMM.sh
# Example using imputed data: runBOLT_v6-LMM.sh


##############################
# 4. Create qsub scripts for each phenotype

for pheno in $phenotypes;
do \
    echo $pheno

    mkdir $basedir/QC-Scripts/GWAS/otherGWAS/$pheno
    mkdir $basedir/QC-Scripts/GWAS/otherGWAS/$pheno/Logs
    mkdir $basedir/QC-Scripts/GWAS/otherGWAS/$pheno/plots

    # copy template
    cp $basedir/QC-Scripts/GWAS/otherGWAS/ $basedir/QC-Scripts/GWAS/otherGWAS/$pheno
    
    # update values in template for this phenotype
    sed "s/phenotypeForThisFile/$pheno/g" $basedir/QC-Scripts/GWAS/otherGWAS/$pheno/runBOLT_$vers-LMM.sh  | sed "s/versionNumber/$vers/g" > $basedir/QC-Scripts/GWAS/otherGWAS/$pheno/runBOLT_$vers-LMM-$pheno.sh

    # remove template from phenotype folder
    rm $basedir/QC-Scripts/GWAS/otherGWAS/$pheno/runBOLT_$vers-LMM.sh

done;
    


##############################
# 5. IMPUTATION BGEN VERSIONS ONLY: Create dummy plink file that has all the right samples. Make sure this matches the parameter --bfile in the qsub scripts

if[ "$vers" == "v6" ]
then \
    cp /well/ukbiobank/phasing/final/phased_chunks/chr2.test2.sample chr2.test2.sample
    cut -f1,2 -d' ' chr2.test2.sample | tail -n +3  > toKeep.tmp
    awk '$1==2' /well/ukbiobank/expt/V2_QCed.export/data/genotype_gwas/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.genotype_gwas.bim | head -n 1 > snp.tmp
    $plink --bfile /well/ukbiobank/expt/V2_QCed.export/data/genotype_gwas/V2_QCed.export.UKBiLEVEAX_b1-b11.Batch_b001_b095.genotype_gwas --extract snp.tmp --keep-allele-order --keep toKeep.tmp --make-bed --out $basedir/data/GWAS/otherGWAS/ImputedBgenSampleSet.$vers.dummy
fi;




##############################
# 6. Submit jobs to cluster (or source)

for pheno in $phenotypes;
do \

cd $basedir/QC-Scripts/GWAS/otherGWAS/$pheno

# FOR EXAMPLE 1: all chromosomes if using genotyped data. Eg. v5
# qsub --q short.qc -P donnelly.prjc runBOLT_$vers-LMM-$pheno.sh

# FOR EXAMPLE 2: by chromosome for imputed data. Eg. v6
# qsub -t 1-22 -q short.qc -P donnelly.prjc runBOLT_$vers-LMM-$pheno.sh

done;


##############################
# 7. Basic manhattan plots

# USAGE: Rscript $basedir/QC-Scripts/GWAS/otherGWAS/plot-BOLT-results-known-hits.R <BOLT output GWAS data> <chromosomes> <output directory for plots> [-qc] [-ymax <max value>] [-hits <RData file with known hits>] [-title <text added to plot filename name>]

# -qc option exclude snps with maf < 0.001, missing > 0.05, info < 0.3 (see the R script). Filter values are always shown in the name of the output file.
# -hits option only works if <chromosome> is not 'genome' and hits file exists

for pheno in $phenotypes;
do \

    cd $basedir/QC-Scripts/GWAS/otherGWAS/$pheno # only necessary if <output directory for plots> is set relative to current directory

    # R data file with known hits (only exist for Standing.height and IOP)
    hits=$basedir/QC-Scripts/GWAS/otherGWAS/GWAScatalogue/hg19/gwasCatalog-subset-$pheno.RData
    
    # EXAMPLES:

    # plot all chromosomes that exist in output data, but in separate plots, and qc-filtered
    Rscript $basedir/QC-Scripts/GWAS/otherGWAS/plot-BOLT-results-known-hits.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out all plots -qc -ymax 50

    # plot all chromosomes that exist in output data, but in one single plot, and qc-filtered
    Rscript $basedir/QC-Scripts/GWAS/otherGWAS/plot-BOLT-results-known-hits.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out genome plots -qc -ymax 50

    # plot all chromosomes that exist in output data, but in one single plot, not qc-filtered
    Rscript $basedir/QC-Scripts/GWAS/otherGWAS/plot-BOLT-results-known-hits.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out genome plots -ymax 50

    # plot all chromosomes that exist in output data, but in separate plots, not qc-filtered
    Rscript $basedir/QC-Scripts/GWAS/otherGWAS/plot-BOLT-results-known-hits.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out all plots -ymax 50 -hits $hits -title -raw-Euro-hits

done;


##############################
# 8. Region plots - find regions of hits and save lists of top hits.

# Same usage as $basedir/QC-Scripts/GWAS/otherGWAS/plot-BOLT-results-known-hits.R except the following extra options:
# -phenoFile <phenotypeFile> is required, but not if using the -dontComputeLD flag
# -dontPlotRegions to not plot each region separately. Only plot inferred regions on each chromsome
# -dontComputeLD if you don't want to compute LD stats (otherwise takes half an hour or so)


for pheno in $phenotypes;
do \

    cd $basedir/QC-Scripts/GWAS/otherGWAS/$pheno  # only necessary if <output directory for plots> is set relative to current directory

    # Find regions and only plot chromosome-wide manhattan plot
    Rscript $basedir/QC-Scripts/GWAS/otherGWAS/region.plots.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out genome plots -ymax 50 -qc -title -recombRegions -dontPlotRegions -dontComputeLD


    # Find regions and show known hits, qc-filtered. Compute LD and plot regions separately
    Rscript $basedir/QC-Scripts/GWAS/otherGWAS/region.plots.R $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$vers/$pheno-BOLT-LMM-$vers.out genome plots -ymax 50 -qc -hits $hits -phenoFile $phenotypeFile -title -recombRegions


done;



##############################
# 9. MANUAL STEP: Clap your hands and say 'yeah, we're awesome!'
