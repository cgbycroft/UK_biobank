#!/bin/bash
#$ -cwd
#$ -o Logs
#$ -j y
#$ -V

Rscript ../plot-BOLT-results-known-hits.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Haemoglobin.Concentration.qnorm/BOLTLMM.v15/Haemoglobin.Concentration.qnorm-BOLT-LMM-v15-chr%%.out all plots -mininfo 0.3 -minmaf 0.001 -ymax 50 -title -QCfiltered
