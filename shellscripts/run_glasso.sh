#!/bin/bash
#SBATCH
#SBATCH --job-name=run_glasso
#SBATCH --time=1:0:0
#SBATCH --partition=defq
#SBATCH --mem=40G


module load r 
homeDir=`pwd | sed -e 's/\/shellscripts//g'`
echo $homeDir
scriptDir=$homeDir/src/03_network_inference
logDir=$homeDir/logs

echo $1
echo $2
echo $3
echo $4
if [ "$3" = "weighted" ]
then 
  echo "weighted covariance aggregation"
  mkdir -p $logDir"/weighted_cov_networks/"$1"/glasso_"$2"/"
  Rscript $scriptDir/18-infer_network.R $1 $2 $3 $4 > $logDir"/weighted_cov_networks/"$1"/glasso_"$2"/"$4"_lambda.log"
fi

if [ "$3" = "unweighted" ]
then
  echo "unweighted covariance aggregation"
  mkdir -p $logDir"/unweighted_cov_networks/"$1"/glasso_"$2"/"
  Rscript $scriptDir/18-infer_network.R $1 $2 $3 $4 > "$logDir/unweighted_cov_networks/"$1"/glasso_"$2"/"$4"_lambda.log"
fi

if [ "$3" = "pooling_first" ]
then 
  echo "correcting after pooling"
  mkdir -p $logDir"/correcting_after_pooling/"$1"/glasso_"$2"/"
  Rscript $scriptDir/18-infer_network.R $1 $2 $3 $4 > "$logDir/correcting_after_pooling/"$1"/glasso_"$2"/"$4"_lambda.log"
fi

if [ "$3" = "correcting_first" ]
then
  echo "pooling after correcting"
  mkdir -p $logDir"/pooling_after_correcting/"$1"/glasso_"$2"/"
  Rscript $scriptDir/18-infer_network.R $1 $2 $3 $4 > "$logDir/pooling_after_correcting/"$1"/glasso_"$2"/"$4"_lambda.log"
fi


