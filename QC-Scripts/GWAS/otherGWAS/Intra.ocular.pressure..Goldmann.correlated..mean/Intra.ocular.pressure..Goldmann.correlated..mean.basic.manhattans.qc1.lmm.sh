#!/bin/bash
#$ -cwd
#$ -o Logs
#$ -j y
#$ -V

Rscript ../plot-BOLT-results-known-hits.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Intra.ocular.pressure..Goldmann.correlated..mean/BOLTLMM.v15/Intra.ocular.pressure..Goldmann.correlated..mean-BOLT-LMM-v15-chr%%.out all plots -mininfo 0.3 -ymax 50 -title -QCfiltered
