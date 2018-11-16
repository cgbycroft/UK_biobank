#!/bin/bash
#$ -cwd
#$ -o Logs
#$ -j y
#$ -V

Rscript ../plot-BOLT-results-known-hits.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/Internal.Pico..ng.uL./BOLTLMM.v15/Internal.Pico..ng.uL.-BOLT-LMM-v15-chr%%.out all plots -mininfo 0.3 -ymax 50 -lreg -title -QCfiltered
