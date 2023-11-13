#!/bin/bash
#SBATCH
#SBATCH --job-name=run_excess_overlap
#SBATCH --time=2:0:0
#SBATCH --partition=defq
#SBATCH --mem=50G


module load r 
homeDir=`pwd | sed -e 's/\/shellscripts//g'`
echo $homeDir
scriptDir=$homeDir/src/05_evaluate_biological_properties
logDir=$homeDir/logs


Rscript $scriptDir/25-excess_overlap_known_geneSets.R 
