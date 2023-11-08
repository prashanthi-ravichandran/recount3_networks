#!/bin/bash -l
#SBATCH
#SBATCH --job-name=processReplicates
#SBATCH --mail-type=END
#SBATCH --mail-user=pravich2@jhu.edu
#SBATCH --time=24:0:0
#SBATCH --partition=lrgmem
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=8

module load gcc/5.5.0
module load R/4.0.2

homeDir=`pwd | sed -e 's/\/shellscripts//g'`
echo $homeDir
logDir=$homeDir/log
srcDir=$homeDir/src/01_data_preprocessing


Rscript $srcDir/processing_smallrna_replicates.R > $logDir/process_replicates.log

