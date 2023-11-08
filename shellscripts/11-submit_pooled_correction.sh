#!/bin/bash -l
#SBATCH
#SBATCH --job-name=submit_pooled_correction
#SBATCH --mail-type=END
#SBATCH --mail-user=pravich2@jhu.edu
#SBATCH --time=0:10:0
#SBATCH --partition=express
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

sbatch run_correct_pooled.sh GTEx
sbatch run_correct_pooled.sh sra_normal

