#!/bin/bash
#SBATCH
#SBATCH --job-name=run_SNAP_enrichment
#SBATCH --time=1:0:0
#SBATCH --partition=defq
#SBATCH --mem=40G


module load r 
homeDir=`pwd | sed -e 's/\/shellscripts//g'`
echo $homeDir
scriptDir=$homeDir/src/04_evaluate_data_aggregation
logDir=$homeDir/logs


Rscript $scriptDir/22-compute_SNAP_encrichment.R $1 $2 $3
