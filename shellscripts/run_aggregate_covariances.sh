#!/bin/bash -l
#SBATCH
#SBATCH --job-name=aggregateCov
#SBATCH --mail-type=END
#SBATCH --mail-user=pravich2@jhu.edu
#SBATCH --time=72:0:0
#SBATCH --partition=defq
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=2

module load r-sva

homeDir=`pwd | sed -e 's/\/shellscripts//g'`
echo $homeDir
logDir=$homeDir/log
srcDir=$homeDir/src/03_network_inference


Rscript $srcDir/17-aggregating_individual_covariances.R $1 > $logDir"/aggregate_covariances_"$1".log"
