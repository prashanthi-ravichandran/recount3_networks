#!/bin/bash
#SBATCH
#SBATCH --job-name=run_GO_enrichment
#SBATCH --time=1:0:0
#SBATCH --partition=defq
#SBATCH --mem=40G


module load r 
homeDir=`pwd | sed -e 's/\/shellscripts//g'`
echo $homeDir
scriptDir=$homeDir/src/05_evaluate_biological_properties
logDir=$homeDir/logs


Rscript $scriptDir/24-enrichment_specific_pathways.R 
