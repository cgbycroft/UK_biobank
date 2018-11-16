Rscript ../plot-GWAS-imputation-comparison.R /data/GWAS/otherGWAS/Array.as.binary/BOLTLMM.v14/Array.as.binary-BOLT-LMM-v14.out $SGE_TASK_ID plots -ymax 50 -phenoFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt -lreg -title -raw -bgenFile /data/GWAS/otherGWAS/Array.as.binary/BOLTLMM./Array.as.binary-BOLT-LMM--chr$SGE_TASK_ID.out -sample /well/ukbiobank/imputation/final/full/bgen/chr$SGE_TASK_ID.hrc+uk10k.I4.v1.1.sample -dontComputeLD