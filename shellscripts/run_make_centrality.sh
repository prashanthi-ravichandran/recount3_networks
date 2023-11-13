#!/bin/bash
#SBATCH
#SBATCH --job-name=run_glasso
#SBATCH --time=1:0:0
#SBATCH --partition=defq
#SBATCH --mem=40G


module load r 
homeDir=`pwd | sed -e 's/\/shellscripts//g'`
echo $homeDir
scriptDir=$homeDir/src/08_stratified_LDSC
logDir=$homeDir/logs

Rscript $scriptDir/29-make_centrality_bed.R $1 $2 $3