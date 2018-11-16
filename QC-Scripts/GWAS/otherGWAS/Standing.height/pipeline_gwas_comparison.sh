################
# Pipeline for running analysis of credible regions, hit overlaps, etc. between GIANT and standing height.
################

basedir=/well/ukbiobank/qcoutput.V2_QCed.sample-QC
pheno=Standing.height
cd $basedir/QC-Scripts/GWAS/otherGWAS/$pheno


########################################################################
## 0: Define versions of GWAS results to use.
########################################################################


#versImp=v19
versImp=v20
versImpb=v20b
versGen=v16
#chr=2

# Genotyped UKBiobank
UKBBGenotypedGWAS=$basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$versGen/$pheno-BOLT-LMM-$versGen.out
# Imputed UKBiobank
UKBBImputedGWAS=$basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$versImp/$pheno-BOLT-LMM-${versImp}-chr%%.out
UKBBImputedGWAS_sub=$basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$versImpb/$pheno-BOLT-LMM-${versImpb}-chr%%.out
# GIANT
GIANTGWAS=$basedir/QC-Scripts/GWAS/otherGWAS/Standing.height/GIANT/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq-with-positions.txt.gz

# Sanity check of numbers of SNPs in imputed data
for chr in {1..4}; # <== Still to check 1-4 once done!
do \
    echo $chr
    grep 'BGEN snpBlocks' Logs/bolt_lmm_$versImp.o*.$chr
    wc -l $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$versImp/$pheno-BOLT-LMM-${versImp}-chr$chr.out
done;
# ====> All correct for v20 5-22; and v20b 1-22

# Subset these for significant SNPs only (do once. convenience for reading quickly later)

for chr in {1..22};
do \
    echo $chr
    #file=$basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$versImp/$pheno-BOLT-LMM-${versImp}-chr$chr.out
    #cat <(head -n 1 $file) <(awk '$14 < 5e-8'  $file ) > $file.signif
    
    file=$basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$versImpb/$pheno-BOLT-LMM-${versImpb}-chr$chr.out
    cat <(head -n 1 $file) <(awk '$14 < 5e-8'  $file ) > $file.signif

done;

cat <(head -n 1 $UKBBGenotypedGWAS ) <( awk '$14 < 5e-8' $UKBBGenotypedGWAS ) > $UKBBGenotypedGWAS.signif
zcat $GIANTGWAS | awk -F' ' '$11 < 5e-8' > $GIANTGWAS.signif

######
# Change permissions on GWAS output data so they can't be accidentally overwritten (but only once jobs have finished running!)

chmod -w $UKBBGenotypedGWAS

for chr in {1..22};
do \
    echo $chr
    chmod -w $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$versImp/$pheno-BOLT-LMM-${versImp}-chr$chr.out
    chmod -w $basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$versImp/$pheno-BOLT-LMM-${versImp}-chr$chr.out.geno
done;


########################################################################
# 0b: Count snps using different filter criteria
########################################################################

rm $outdir/$versImp.filtered-mininfo0.3-minmaf0.001.counts
touch $outdir/$versImp.filtered-mininfo0.3-minmaf0.001.counts
for chr in {1..22};
do \
    echo $chr
    file=$basedir/data/GWAS/otherGWAS/$pheno/BOLTLMM.$versImp/$pheno-BOLT-LMM-${versImp}-chr$chr.out
    #    8==INFO
    #    7 == A1FREQ
    awk ' $8 > 0.3 && $7 > 0.001 && $7 < 0.999 ' $file | wc -l >> $outdir/$versImp.filtered-mininfo0.3-minmaf0.001.counts
    
done;

zcat $GIANTGWAS | awk ' $5 > 0.001 && $5 < 0.999 ' | wc -l > $outdir/GIANT.filtered-mininfo0.3-minmaf0.001.counts

awk ' $7 > 0.001 && $7 < 0.999 ' $UKBBGenotypedGWAS  | wc -l  > $outdir/GENO.filtered-mininfo0.3-minmaf0.001.counts



########################################################################
# 1: Make a venn diagram of hit regions across different datasets. 
########################################################################

## First define a set of regions (simple definition using 1Mb windows).

Rscript define-regions-windows.R -ws 1 # 1Mb window sizes
# /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/regions-hg19-window-1Mb.txt


regions=/well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/GWAS/otherGWAS/regions-hg19-window-1Mb.txt
title=1MbWin

outdir=$basedir/data/GWAS/otherGWAS/$pheno/hitCounts
mkdir $outdir



## Use this set of regions to count the number of hits in each region for three GWAS sets
cd $basedir/QC-Scripts/GWAS/otherGWAS/$pheno

cat ../submit.header.txt <( echo Rscript ../count-hit-regions-basic.R $UKBBGenotypedGWAS $regions -title .$title.GENO -chr '$SGE_TASK_ID' -outdir $outdir ) > submit.$pheno.count-hit-regions-basic.$title.geno.sh

cat ../submit.header.txt <( echo Rscript ../count-hit-regions-basic.R $UKBBImputedGWAS $regions -title .$title.IMP -chr '$SGE_TASK_ID' -outdir $outdir ) > submit.$pheno.count-hit-regions-basic.$title.imp.sh

# For v20b
cat ../submit.header.txt <( echo Rscript ../count-hit-regions-basic.R $UKBBImputedGWAS_sub $regions -title .$title.IMPSUB -chr '$SGE_TASK_ID' -outdir $outdir ) > submit.$pheno.count-hit-regions-basic.$title.impsub.sh

cat ../submit.header.txt <( echo Rscript ../count-hit-regions-basic.R $GIANTGWAS $regions -title .$title.GIANT -chr '$SGE_TASK_ID' -outdir $outdir ) > submit.$pheno.count-hit-regions-basic.$title.giant.sh

# submit to cluster (still to run other chromosomes fomr v120 imputation)

qsub -t 1-22 -N count-hit-regions-basic.$title.geno -q short.qc -P donnelly.prjc submit.$pheno.count-hit-regions-basic.$title.geno.sh
qsub -t 1-22 -N count-hit-regions-basic.$title.imp -q short.qc -P donnelly.prjc submit.$pheno.count-hit-regions-basic.$title.imp.sh
qsub -t 1-22 -N count-hit-regions-basic.$title.impsub -q short.qc -P donnelly.prjc submit.$pheno.count-hit-regions-basic.$title.impsub.sh
qsub -t 1-22 -N count-hit-regions-basic.$title.giant -q short.qc -P donnelly.prjc submit.$pheno.count-hit-regions-basic.$title.giant.sh


## Plot the counts in these regions.  Set up for imputation/genotyped/giant

##########
genoCounts=$outdir/$pheno-BOLT-LMM-$versGen.out.chr%%.maf0.miss1.pruned.$title.GENO-5e-08.txt
impCounts=$outdir/$pheno-BOLT-LMM-$versImp-chr%%.out.chr%%.maf0.info0.pruned.$title.IMP-5e-08.txt
giantCounts=$outdir/GIANT.chr%%.$title.GIANT-5e-08.txt
##########

Rscript ../plot-count-hit-regions-basic.R \
	-plotdir $basedir/QC-Scripts/GWAS/otherGWAS/Standing.height \
	-counts $genoCounts,$impCounts,$giantCounts \
	-chr 1-22 -title hit-counts-Imp.$versImp.Geno.$versGen.GIANT-$title \
	-outdir $outdir \
    &> Logs/plot-count-hit-regions-basic.noqc.log &

Rscript ../plot-count-hit-regions-basic.R \
	-plotdir $basedir/QC-Scripts/GWAS/otherGWAS/Standing.height \
	-counts $genoCounts,$impCounts,$giantCounts \
	-chr 1-22 -title hit-counts-Imp.$versImp.Geno.$versGen.GIANT-$title \
	-outdir $outdir \
	-mininfo 0.3 \
    &> Logs/plot-count-hit-regions-basic.info0.3.log &

Rscript ../plot-count-hit-regions-basic.R \
	-plotdir $basedir/QC-Scripts/GWAS/otherGWAS/Standing.height \
	-counts $genoCounts,$impCounts,$giantCounts \
	-chr 1-22 -title hit-counts-Imp.$versImp.Geno.$versGen.GIANT-$title \
	-outdir $outdir \
	-minmaf 0.001 \
    &> Logs/plot-count-hit-regions-basic.maf0.001.log &

Rscript ../plot-count-hit-regions-basic.R \
	-plotdir $basedir/QC-Scripts/GWAS/otherGWAS/Standing.height \
	-counts $genoCounts,$impCounts,$giantCounts \
	-chr 1-22 -title hit-counts-Imp.$versImp.Geno.$versGen.GIANT-$title \
	-outdir $outdir \
	-mininfo 0.3 \
	-minmaf 0.001 \
    &> Logs/plot-count-hit-regions-basic.info0.3.maf0.001.log &


Rscript ../plot-count-hit-regions-basic.R \
	-plotdir $basedir/QC-Scripts/GWAS/otherGWAS/Standing.height \
	-counts $genoCounts,$impCounts,$giantCounts \
	-chr 1-22 -title hit-counts-Imp.$versImp.Geno.$versGen.GIANT-$title \
	-outdir $outdir \
	-mininfo 0.3 \
	-minmaf 0.001 \
	-maxmiss 0.05 \
    &> Logs/plot-count-hit-regions-basic.info0.3.maf0.001.miss0.05.log &



############################################################
# 2: Credible set analysis
############################################################

#### Compute bayes factors for all snps (relatively slow)
myPrior=0.2

# For imputed data go over chromosomes.
cat ../submit.header.txt <( echo Rscript ../compute-bayes-factors.R -gwas $UKBBImputedGWAS -chr '$SGE_TASK_ID' -prior $myPrior ) > submit.compute-bayes-factors.imp.$versImp.sh

qsub -t 1-22 -N compute-bayes-factors.imp.$versImp.$myPrior -q short.qc -P donnelly.prjc submit.compute-bayes-factors.imp.$versImp.sh

# For imputed data (subset)  go over chromosomes.
cat ../submit.header.txt <( echo Rscript ../compute-bayes-factors.R -gwas $UKBBImputedGWAS_sub -chr '$SGE_TASK_ID' -prior $myPrior ) > submit.compute-bayes-factors.impsub.$versImpb.sh

qsub -t 1-22 -N compute-bayes-factors.impsub.$versImpb.$myPrior -q short.qc -P donnelly.prjc submit.compute-bayes-factors.impsub.$versImpb.sh

# Genotypes
Rscript ../compute-bayes-factors.R -gwas $UKBBGenotypedGWAS -prior $myPrior &> Logs/compute-bayes-factors.geno.$versGen.$myPrior.log &

# Giant
Rscript ../compute-bayes-factors.R -gwas $GIANTGWAS -prior $myPrior &> Logs/compute-bayes-factors.giant.$myPrior.log &


#### Then define the regions (use GIANT)

Rscript ../define-regions-based-on-hits.R -gwas $GIANTGWAS.signif -chr 1-22 -outdir $outdir &> Logs/define-regions-based-on-hits.GIANT.log &

Rscript ../define-regions-based-on-hits.R -gwas $UKBBGenotypedGWAS.signif -chr 1-22 -outdir $outdir &> Logs/define-regions-based-on-hits.GENO.log &



###########################
## Subset by region, and compute posteriors for these SNPs 

regionsBF=$outdir/regions-0.125cM-25KB-GIANT.chrgenome.txt
#chr=5

cat ../submit.header.txt <( echo Rscript ../compute-posteriors.R -gwas $UKBBImputedGWAS -prior $myPrior -regions $regionsBF -outdir $outdir -chr '$SGE_TASK_ID' -title .giantRegs ) > submit.compute-posteriors.imp.sh

cat ../submit.header.txt <( echo Rscript ../compute-posteriors.R -gwas $UKBBImputedGWAS_sub -prior $myPrior -regions $regionsBF -outdir $outdir -chr '$SGE_TASK_ID' -title .giantRegs ) > submit.compute-posteriors.impsub.sh

cat ../submit.header.txt <( echo Rscript ../compute-posteriors.R -gwas $UKBBGenotypedGWAS -prior $myPrior -regions $regionsBF -outdir $outdir -chr '$SGE_TASK_ID' -title .giantRegs ) > submit.compute-posteriors.geno.sh

cat ../submit.header.txt <( echo Rscript ../compute-posteriors.R -gwas $GIANTGWAS -prior $myPrior -regions $regionsBF -outdir $outdir -chr '$SGE_TASK_ID' -title .giantRegs ) > submit.compute-posteriors.giant.sh


qsub -t 1-22 -N compute-posteriors.giantRegs.impsub.$myPrior -q short.qc -P donnelly.prjc submit.compute-posteriors.impsub.sh
qsub -t 1-22 -N compute-posteriors.giantRegs.imp.$myPrior -q short.qc -P donnelly.prjc submit.compute-posteriors.imp.sh
qsub -t 1-22 -N compute-posteriors.giantRegs.geno.$myPrior -q short.qc -P donnelly.prjc submit.compute-posteriors.geno.sh
qsub -t 1-22 -N compute-posteriors.giantRegs.giant.$myPrior -q short.qc -P donnelly.prjc submit.compute-posteriors.giant.sh


# Check the logs
grep filtered Logs/compute-posteriors.giantRegs*.o*.1


#################################
## Merge the regions

########
v=version1
#myPrior=0.2
# NO QC ON THE SNPS; Prior = 0.2
impPosteriors=Posteriors-Standing.height-BOLT-LMM-$versImp-chr%%.out.chr%%.maf0.info0.pruned.giantRegs.$myPrior.RData
impsubPosteriors=Posteriors-Standing.height-BOLT-LMM-$versImpb-chr%%.out.chr%%.maf0.info0.pruned.giantRegs.$myPrior.RData
genoPosteriors=Posteriors-Standing.height-BOLT-LMM-$versGen.out.chr%%.maf0.miss1.pruned.giantRegs.$myPrior.RData
giantPosteriors=Posteriors-GIANT.chr%%.giantRegs.$myPrior.RData
########


Rscript ../merge-posteriors.R \
	-posteriors $impPosteriors,$genoPosteriors,$giantPosteriors \
	-outdir $outdir \
	-chr 1-22 \
	-title $v.$myPrior \
    &> Logs/merge-posteriors.$v.$myPrior.log &

grep 'ERROR|WARNING' Logs/merge-posteriors.$v.$myPrior.log


###########
## Calculate credible sets



############################################
## Download Neale's data for comparison!

mkdir $basedir/QC-Scripts/GWAS/otherGWAS/Standing.height/Neale2018
cd $basedir/QC-Scripts/GWAS/otherGWAS/Standing.height/Neale2018


wget https://www.dropbox.com/s/ou12jm89v74k55e/50_irnt.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 50_irnt.gwas.imputed_v3.both_sexes.tsv.bgz

# THE MANIFEST FOR THIS DATA IS HERE:


