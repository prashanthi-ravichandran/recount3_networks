#!/bin/bash
#SBATCH
#SBATCH --job-name=run_SNAP_enrichment_perm
#SBATCH --time=24:0:0
#SBATCH --partition=shared
#SBATCH --mem=40G


module load r/4.2.0
homeDir=`pwd | sed -e 's/\/shellscripts//g'`
echo $homeDir
scriptDir=$homeDir/src/04_evaluate_data_aggregation/

Rscript $scriptDir/SNAP_permutation_test_2.R $1 $2 $3
