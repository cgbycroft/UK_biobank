Rscript ../plot-GWAS-imputation-comparison.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/dQC.qnorm/BOLTLMM.v16s/dQC.qnorm-BOLT-LMM-v16s.out $SGE_TASK_ID plots -ymax 50 -phenoFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt -lreg -title -raw -bgenFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/dQC.qnorm/BOLTLMM.v15s/dQC.qnorm-BOLT-LMM-v15s-chr$SGE_TASK_ID.out -sample /well/ukbiobank/imputation/final/full/bgen/hrc+uk10k_sorted_v1.1/chr$SGE_TASK_ID.hrc+uk10k.I4.v1.1.sample -dontComputeLD
