#!/bin/bash
#SBATCH
#SBATCH --job-name=run_glasso_single
#SBATCH --time=24:0:0
#SBATCH --partition=shared
#SBATCH --mem=20G


module load r/4.0.2 
homeDir=`pwd | sed -e 's/\/shellscripts//g'`
echo $homeDir
scriptDir=$homeDir/src/03_network_inference

echo $1
Rscript $scriptDir/18-infer_network_single_study.R $1 



