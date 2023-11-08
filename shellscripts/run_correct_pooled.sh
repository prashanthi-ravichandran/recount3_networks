#!/bin/bash -l
#SBATCH
#SBATCH --job-name=correctPooled
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
srcDir=$homeDir/src/02_data_correction


Rscript $srcDir/16-correcting_pooled_studies.R $1 > $logDir"/correct_pooled_expression_"$1".log"
