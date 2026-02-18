#!/bin/bash
#SBATCH
#SBATCH --job-name=consensus_WGCNA
#SBATCH --time=24:0:0
#SBATCH --partition=bigmem
#SBATCH -A abattle44_bigmem
#SBATCH --mem=250G

module load r/4.2.0
homeDir=`pwd | sed -e 's/\/shellscripts//g'`
scriptDir=$homeDir/src/10_WGCNA
echo $1
echo $2
Rscript $scriptDir/05-consensus_network_inference.R $1 $2
