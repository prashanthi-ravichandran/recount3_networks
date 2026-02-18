#!/bin/bash
#SBATCH
#SBATCH --job-name=run_TF_enr
#SBATCH --time=2:0:0
#SBATCH --partition=shared
#SBATCH --mem=60G

module load r/4.0.2
homeDir=`pwd | sed -e 's/\/shellscripts//g'`
echo $homeDir
scriptDir=$homeDir/src/04_evaluate_data_aggregation

Rscript $scriptDir/evaluate_TF_target_enrichment.R $1 $2 
