

###########
# v2 : Y l2r mean (instead of Affy's NP probes)
pheno=Y.mean.l2r
vers=v2

###########
# v3 : X l2r mean (instead of Affy's NP probes)
pheno=X.mean.l2r
vers=v3

###########
# v4 : X l2r mean (instead of Affy's NP probes) ==> sex chroms
pheno=X.mean.l2r
vers=v4

###########
# v5 : Y l2r mean (instead of Affy's NP probes) ==> sex chroms
pheno=Y.mean.l2r
vers=v5

###########
# v6 : X-Y l2r mean (females) ==> sex chroms and autosomes in same version. No 24, obviously.
pheno=XYdiff.mean.l2r
vers=v6

###########
# v7 : X0 or XX (only 135 samples!)
pheno=X0
vers=v7


# for autosomes
Rscript ../otherGWAS/plot-BOLT-results-known-hits.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/YXintensity/BOLTLMM.$vers/Ychrom-BOLT-LMM-quant-Age-$vers.out all plots -qc -ymax 50 -title -$pheno-QCfiltered
Rscript ../otherGWAS/plot-BOLT-results-known-hits.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/YXintensity/BOLTLMM.$vers/Ychrom-BOLT-LMM-quant-Age-$vers.out all plots -qc -ymax 50 -minmaf 0 -maxmiss 0.05 -title -$pheno-QCfiltered2
Rscript ../otherGWAS/plot-BOLT-results-known-hits.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/YXintensity/BOLTLMM.$vers/Ychrom-BOLT-LMM-quant-Age-$vers.out all plots -ymax 50 -title -$pheno-raw


# for sex chroms
Rscript ../otherGWAS/plot-BOLT-results-known-hits.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/YXintensity/BOLTLMM.$vers/Ychrom-BOLT-LMM-quant-Age-$vers-%%.out 23-26 plots -qc -ymax 50 -title -$pheno-QCfiltered -sex
Rscript ../otherGWAS/plot-BOLT-results-known-hits.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/YXintensity/BOLTLMM.$vers/Ychrom-BOLT-LMM-quant-Age-$vers-%%.out 23-26 plots -ymax 50 -title -$pheno-raw -sex
