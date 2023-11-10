#!/bin/bash -l
#SBATCH
#SBATCH --job-name=submit_scale_free
#SBATCH --mail-type=END
#SBATCH --mail-user=pravich2@jhu.edu
#SBATCH --time=0:20:0
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

sbatch run_evaluate_scale_free.sh all_consensus
sbatch run_evaluate_scale_free.sh normal_consensus
sbatch run_evaluate_scale_free.sh cancer_consensus
sbatch run_evaluate_scale_free.sh GTEx
sbatch run_evaluate_scale_free.sh sra_normal

sbatch run_evaluate_scale_free.sh airway/all 
sbatch run_evaluate_scale_free.sh blood/all 
sbatch run_evaluate_scale_free.sh B_cells/all 
sbatch run_evaluate_scale_free.sh breast/all 
sbatch run_evaluate_scale_free.sh cardiac/all 
sbatch run_evaluate_scale_free.sh central_nervous_system/all 
sbatch run_evaluate_scale_free.sh colon/all 
sbatch run_evaluate_scale_free.sh eye/all 
sbatch run_evaluate_scale_free.sh esophagus/all 
sbatch run_evaluate_scale_free.sh fibroblasts/all 
sbatch run_evaluate_scale_free.sh hescs/all 
sbatch run_evaluate_scale_free.sh intestine/all 
sbatch run_evaluate_scale_free.sh ipscs/all 
sbatch run_evaluate_scale_free.sh kidney/all 
sbatch run_evaluate_scale_free.sh liver/all 
sbatch run_evaluate_scale_free.sh lung/all 
sbatch run_evaluate_scale_free.sh multipotent_cells/all 
sbatch run_evaluate_scale_free.sh myeloid_cells/all 
sbatch run_evaluate_scale_free.sh pancreas/all 
sbatch run_evaluate_scale_free.sh pbmcs_t_cells/all 
sbatch run_evaluate_scale_free.sh prostate/all 
sbatch run_evaluate_scale_free.sh skeletal_muscle/all 
sbatch run_evaluate_scale_free.sh skin/all 
sbatch run_evaluate_scale_free.sh stomach/all 
sbatch run_evaluate_scale_free.sh vascular/all 

sbatch run_evaluate_scale_free.sh adipose/GTEx 
sbatch run_evaluate_scale_free.sh blood/GTEx 
sbatch run_evaluate_scale_free.sh B_cells/GTEx 
sbatch run_evaluate_scale_free.sh breast/GTEx 
sbatch run_evaluate_scale_free.sh cardiac/GTEx 
sbatch run_evaluate_scale_free.sh central_nervous_system/GTEx 
sbatch run_evaluate_scale_free.sh colon/GTEx 
sbatch run_evaluate_scale_free.sh esophagus/GTEx 
sbatch run_evaluate_scale_free.sh fibroblasts/GTEx 
sbatch run_evaluate_scale_free.sh intestine/GTEx 
sbatch run_evaluate_scale_free.sh kidney/GTEx 
sbatch run_evaluate_scale_free.sh liver/GTEx 
sbatch run_evaluate_scale_free.sh lung/GTEx 
sbatch run_evaluate_scale_free.sh pancreas/GTEx 
sbatch run_evaluate_scale_free.sh prostate/GTEx 
sbatch run_evaluate_scale_free.sh skeletal_muscle/GTEx 
sbatch run_evaluate_scale_free.sh skin/GTEx 
sbatch run_evaluate_scale_free.sh stomach/GTEx 
sbatch run_evaluate_scale_free.sh vascular/GTEx 


