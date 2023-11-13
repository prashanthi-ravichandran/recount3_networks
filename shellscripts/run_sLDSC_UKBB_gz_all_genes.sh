#!/bin/bash
#SBATCH
#SBATCH --job-name=run_ldsc_UKBB_gz_all_genes
#SBATCH --time=2:0:0
#SBATCH --partition=defq
#SBATCH --mem=60G

module load anaconda
conda activate ldsc 

STUDY=$1
CENTRALITY=$2
TRAIT=$3

echo $STUDY
echo $CENTRALITY
echo $TRAIT

python /data/abattle4/prashanthi/recount3_paper/src/08_stratified_LDSC/ldsc.py --h2 "/data/abattle4/prashanthi/recount3_paper/data/s_LDSC/UKBB/"$TRAIT".imputed_v3.ldsc.both_sexes.tsv.gz" --w-ld-chr /data/abattle4/prashanthi/recount3_paper/data/s_LDSC/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. --ref-ld-chr /data/abattle4/prashanthi/recount3_paper/data/s_LDSC/all-genes/SAMGENIC.,"/data/abattle4/prashanthi/recount3_paper/results/s_LDSC/"$STUDY"/"$CENTRALITY"/ldscore/"$STUDY"." --overlap-annot --frqfile-chr /data/abattle4/prashanthi/recount3_paper/data/s_LDSC/1000G_Phase3_frq/1000G.EUR.QC. --out "/data/abattle4/prashanthi/recount3_paper/results/s_LDSC/"$STUDY"/"$CENTRALITY"/results/UKBB/"$TRAIT"_all_genes" --print-coefficients