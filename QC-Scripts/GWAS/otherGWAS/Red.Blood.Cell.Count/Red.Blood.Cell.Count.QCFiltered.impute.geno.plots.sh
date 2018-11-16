Rscript ../plot-GWAS-imputation-comparison.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Red.Blood.Cell.Count/BOLTLMM.v16s/Red.Blood.Cell.Count-BOLT-LMM-v16s.out $SGE_TASK_ID plots -ymax 50 -qc -phenoFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/PhenotypesForBOLT-v3.txt -lreg -title -QCFiltered -bgenFile /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Red.Blood.Cell.Count/BOLTLMM.v15s/Red.Blood.Cell.Count-BOLT-LMM-v15s-chr$SGE_TASK_ID.out -sample /well/ukbiobank/imputation/final/full/bgen/hrc+uk10k_sorted_v1.1/chr$SGE_TASK_ID.hrc+uk10k.I4.v1.1.sample -dontComputeLD
