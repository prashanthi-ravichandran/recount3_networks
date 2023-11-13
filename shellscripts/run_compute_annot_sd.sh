#!/bin/bash -l
#SBATCH
#SBATCH --job-name=compute_SD
#SBATCH --time=0:15:0
#SBATCH --partition=defq
#SBATCH --mem=20G

module load r

Rscript /data/abattle4/prashanthi/recount3_paper/src/08_stratified_LDSC/32-compute_annot_SD.R $1 $2 $3

