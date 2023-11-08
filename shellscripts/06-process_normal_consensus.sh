#!/bin/bash -l
#SBATCH
#SBATCH --job-name=processNormalConsensus
#SBATCH --mail-type=END
#SBATCH --mail-user=pravich2@jhu.edu
#SBATCH --time=4:0:0
#SBATCH --partition=lrgmem
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=2

module load gcc/5.5.0
module load R/4.0.2

homeDir=`pwd | sed -e 's/\/shellscripts//g'`
echo $homeDir
datDir=$homeDir/data
logDir=$homeDir/log
srcDir=$homeDir/src/01_data_preprocessing


Rscript $srcDir/process_expr_normal_consensus.R > $logDir/normal_consensus.log

