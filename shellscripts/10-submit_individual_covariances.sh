#!/bin/bash -l
#SBATCH
#SBATCH --job-name=submit_ind_covariance
#SBATCH --mail-type=END
#SBATCH --mail-user=pravich2@jhu.edu
#SBATCH --time=0:10:0
#SBATCH --partition=express
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

sbatch run_covariance_individual.sh all_consensus
sbatch run_covariance_individual.sh normal_consensus
sbatch run_covariance_individual.sh cancer_consensus
sbatch run_covariance_individual.sh GTEx
sbatch run_covariance_individual.sh sra_normal
sbatch run_covariance_individual.sh adipose/all
sbatch run_covariance_individual.sh adipose/GTEx
sbatch run_covariance_individual.sh airway/all
sbatch run_covariance_individual.sh B_cells/all
sbatch run_covariance_individual.sh B_cells/GTEx
sbatch run_covariance_individual.sh blood/all
sbatch run_covariance_individual.sh blood/GTEx
sbatch run_covariance_individual.sh breast/all
sbatch run_covariance_individual.sh breast/GTEx
sbatch run_covariance_individual.sh cardiac/all
sbatch run_covariance_individual.sh cardiac/GTEx
sbatch run_covariance_individual.sh central_nervous_system/all
sbatch run_covariance_individual.sh central_nervous_sytem/GTEx
sbatch run_covariance_individual.sh colon/all
sbatch run_covariance_individual.sh colon/GTEx
sbatch run_covariance_individual.sh esophagus/all
sbatch run_covariance_individual.sh esophagus/GTEx
sbatch run_covariance_individual.sh eye/all
sbatch run_covariance_individual.sh fibroblasts/all
sbatch run_covariance_individual.sh fibroblasts/GTEx
sbatch run_covariance_individual.sh hescs/all
sbatch run_covariance_individual.sh intestine/all
sbatch run_covariance_individual.sh intestine/GTEx
sbatch run_covariance_individual.sh ipscs/all
sbatch run_covariance_individual.sh kidney/all
sbatch run_covariance_individual.sh kidney/GTEx
sbatch run_covariance_individual.sh liver/all
sbatch run_covariance_individual.sh liver/GTEx
sbatch run_covariance_individual.sh lung/all
sbatch run_covariance_individual.sh lung/GTEx
sbatch run_covariance_individual.sh multipotent_cells/all
sbatch run_covariance_individual.sh myeloid_cells/all
sbatch run_covariance_individual.sh nervous_system/all
sbatch run_covariance_individual.sh nervous_system/GTEx
sbatch run_covariance_individual.sh pancreas/all
sbatch run_covariance_individual.sh pancreas/GTEx
sbatch run_covariance_individual.sh pbmcs_t_cells/all
sbatch run_covariance_individual.sh prostate/all
sbatch run_covariance_individual.sh prostate/GTEx
sbatch run_covariance_individual.sh skeletal_muscle/all
sbatch run_covariance_individual.sh skeletal_muscle/GTEx
sbatch run_covariance_individual.sh skin/all
sbatch run_covariance_individual.sh skin/GTEx
sbatch run_covariance_individual.sh stomach/all
sbatch run_covariance_individual.sh stomach/GTEx
sbatch run_covariance_individual.sh vascular/all
sbatch run_covariance_individual.sh vascular/GTEx


