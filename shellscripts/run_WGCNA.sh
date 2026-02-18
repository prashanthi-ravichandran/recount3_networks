#!/bin/bash
#SBATCH
#SBATCH --job-name=consensus_WGCNA
#SBATCH --time=24:0:0
#SBATCH --partition=shared
#SBATCH --mem=60G

module load r/4.2.0
homeDir=`pwd | sed -e 's/\/shellscripts//g'`
scriptDir=$homeDir/src/10_WGCNA
echo $1
Rscript $scriptDir/02-infer_networks.R $1 $2
