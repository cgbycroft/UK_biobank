#!/bin/bash
#$ -N submit-moran-10Km-20neibrs
#$ -o Logs
#$ -j y
#$ -cwd
#$ -V
#$ -t 1-40
Rscript /well/ukbiobank/qcoutput.V2_QCed.sample-QC/QC-Scripts/PCA/pca-UKBio/run-Moran-Null.R $SGE_TASK_ID 10000 200 20

