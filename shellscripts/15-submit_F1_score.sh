#!/bin/bash -l
#SBATCH
#SBATCH --job-name=submit_odds
#SBATCH --mail-type=END
#SBATCH --mail-user=pravich2@jhu.edu
#SBATCH --time=0:20:0
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

for K in $(seq 0 10 50);
do
	if [ $K == 0 ]
		then K=1
	fi
	if [ $K == 50 ]
    then K=49
  fi
  sbatch run_odds_ratio.sh GTEx weighted_cov_networks $K
  sbatch run_odds_ratio.sh GTEx unweighted_cov_networks $K
  sbatch run_odds_ratio.sh GTEx correcting_after_pooling $K
  sbatch run_odds_ratio.sh GTEx pooling_after_correcting $K   
done

for K in $(seq 0 100 600);
do
	if [ $K == 0 ]
		then K=1
	fi
	if [ $K == 600 ]
      then K=566
  fi
  sbatch run_F1_score.sh sra_normal weighted_cov_networks $K
  sbatch run_F1_score.sh sra_normal unweighted_cov_networks $K
  sbatch run_F1_score.sh sra_normal correcting_after_pooling $K
  sbatch run_F1_score.sh sra_normal pooling_after_correcting $K  
done

sbatch run_F1_score.sh all_consensus weighted_cov_networks 966
sbatch run_F1_score.sh normal_consensus weighted_cov_networks 629
sbatch run_F1_score.sh cancer_consensus weighted_cov_networks 386
