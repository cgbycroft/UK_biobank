Rscript ../plot-GWAS-imputation-comparison.R /data/GWAS/otherGWAS/Place.of.birth.in.UK...east.co.ordinate/BOLTLMM.v14/Place.of.birth.in.UK...east.co.ordinate-BOLT-LMM-v14.out $SGE_TASK_ID plots -ymax 50 -phenoFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt -lreg -title -raw -bgenFile /data/GWAS/otherGWAS/Place.of.birth.in.UK...east.co.ordinate/BOLTLMM./Place.of.birth.in.UK...east.co.ordinate-BOLT-LMM--chr$SGE_TASK_ID.out -sample /well/ukbiobank/imputation/final/full/bgen/chr$SGE_TASK_ID.hrc+uk10k.I4.v1.1.sample -dontComputeLD