#!/bin/bash
#$ -cwd
#$ -o Logs
#$ -j y
#$ -V

Rscript ../plot-BOLT-results-known-hits.R /well/ukbiobank/qcoutput.V2_QCed.sample-QC/data/GWAS/otherGWAS/dQC.qnorm/BOLTLMM.v15/dQC.qnorm-BOLT-LMM-v15-chr%%.out all plots -ymax 50 -title -raw
