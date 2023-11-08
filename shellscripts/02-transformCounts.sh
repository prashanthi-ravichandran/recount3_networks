#!/bin/bash -l
#SBATCH
#SBATCH --job-name=transformCounts
#SBATCH --mail-type=END
#SBATCH --mail-user=pravich2@jhu.edu
#SBATCH --time=10:0:0
#SBATCH --partition=lrgmem
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=2

module load gcc/5.5.0
module load R/4.0.2

homeDir=`pwd | sed -e 's/\/shellscripts//g'`
echo $homeDir
scriptDir=$homeDir/src
datDir=$homeDir/data
logDir=$homeDir/log
rpkm=$homeDir/data/rpkm
auc=$homeDir/data/auc_scaled
raw=$homeDir/data/raw
srcDir=$homeDir/src/01_data_preprocessing
srametDir=$datDir/sra_metadata/




Rscript $srcDir/transformCounts.R > $logDir/transformCounts.log

