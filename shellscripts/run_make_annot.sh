#!/bin/bash
#SBATCH
#SBATCH --job-name=make_annot
#SBATCH --time=0:15:0
#SBATCH --partition=defq
#SBATCH --mem=8G

module load anaconda
conda activate ldsc 

STUDY=$1
CENTRALITY=$2
CHR=$3

echo $STUDY
echo $CENTRALITY
echo $CHR

python make_annot_cont.py $1 $2 $3