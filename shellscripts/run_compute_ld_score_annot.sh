#!/bin/bash
#SBATCH
#SBATCH --job-name=compute_ld
#SBATCH --time=2:0:0
#SBATCH --partition=defq
#SBATCH --mem=4G

module load anaconda
conda activate ldsc 

STUDY=$1
CENTRALITY=$2
CHR=$3

echo $STUDY
echo $CENTRALITY
echo $CHR

python /data/abattle4/prashanthi/recount3_paper/src/08_stratified_LDSC/ldsc.py --l2 --bfile /data/abattle4/prashanthi/recount3_paper/data/s_LDSC/1000G_EUR_Phase3_plink/1000G.EUR.QC.$CHR --ld-wind-cm 1 --annot "/data/abattle4/prashanthi/recount3_paper/results/s_LDSC/"$STUDY"/"$CENTRALITY"/ldscore/"$STUDY"."$CHR".annot.gz" --out "/data/abattle4/prashanthi/recount3_paper/results/s_LDSC/"$STUDY"/"$CENTRALITY"/ldscore/"$STUDY"."$CHR --print-snps /data/abattle4/prashanthi/recount3_paper/data/s_LDSC/list.txt
