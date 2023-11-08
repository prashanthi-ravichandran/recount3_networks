#!/bin/bash -l
#SBATCH
#SBATCH --job-name=submit_cov_aggregate
#SBATCH --mail-type=END
#SBATCH --mail-user=pravich2@jhu.edu
#SBATCH --time=0:10:0
#SBATCH --partition=express
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

sbatch run_aggregate_covariances.sh all_consensus
sbatch run_aggregate_covariances.sh normal_consensus
sbatch run_aggregate_covariances.sh cancer_consensus
sbatch run_aggregate_covariances.sh GTEx
sbatch run_aggregate_covariances.sh sra_normal
sbatch run_aggregate_covariances.sh adipose/all
sbatch run_aggregate_covariances.sh adipose/GTEx
sbatch run_aggregate_covariances.sh airway/all
sbatch run_aggregate_covariances.sh B_cells/all
sbatch run_aggregate_covariances.sh B_cells/GTEx
sbatch run_aggregate_covariances.sh blood/all
sbatch run_aggregate_covariances.sh blood/GTEx
sbatch run_aggregate_covariances.sh breast/all
sbatch run_aggregate_covariances.sh breast/GTEx
sbatch run_aggregate_covariances.sh cardiac/all
sbatch run_aggregate_covariances.sh cardiac/GTEx
sbatch run_aggregate_covariances.sh central_nervous_system/all
sbatch run_aggregate_covariances.sh central_nervous_sytem/GTEx
sbatch run_aggregate_covariances.sh colon/all
sbatch run_aggregate_covariances.sh colon/GTEx
sbatch run_aggregate_covariances.sh esophagus/all
sbatch run_aggregate_covariances.sh esophagus/GTEx
sbatch run_aggregate_covariances.sh eye/all
sbatch run_aggregate_covariances.sh fibroblasts/all
sbatch run_aggregate_covariances.sh fibroblasts/GTEx
sbatch run_aggregate_covariances.sh hescs/all
sbatch run_aggregate_covariances.sh intestine/all
sbatch run_aggregate_covariances.sh intestine/GTEx
sbatch run_aggregate_covariances.sh ipscs/all
sbatch run_aggregate_covariances.sh kidney/all
sbatch run_aggregate_covariances.sh kidney/GTEx
sbatch run_aggregate_covariances.sh liver/all
sbatch run_aggregate_covariances.sh liver/GTEx
sbatch run_aggregate_covariances.sh lung/all
sbatch run_aggregate_covariances.sh lung/GTEx
sbatch run_aggregate_covariances.sh multipotent_cells/all
sbatch run_aggregate_covariances.sh myeloid_cells/all
sbatch run_aggregate_covariances.sh nervous_system/all
sbatch run_aggregate_covariances.sh nervous_system/GTEx
sbatch run_aggregate_covariances.sh pancreas/all
sbatch run_aggregate_covariances.sh pancreas/GTEx
sbatch run_aggregate_covariances.sh pbmcs_t_cells/all
sbatch run_aggregate_covariances.sh prostate/all
sbatch run_aggregate_covariances.sh prostate/GTEx
sbatch run_aggregate_covariances.sh skeletal_muscle/all
sbatch run_aggregate_covariances.sh skeletal_muscle/GTEx
sbatch run_aggregate_covariances.sh skin/all
sbatch run_aggregate_covariances.sh skin/GTEx
sbatch run_aggregate_covariances.sh stomach/all
sbatch run_aggregate_covariances.sh stomach/GTEx
sbatch run_aggregate_covariances.sh vascular/all
sbatch run_aggregate_covariances.sh vascular/GTEx

