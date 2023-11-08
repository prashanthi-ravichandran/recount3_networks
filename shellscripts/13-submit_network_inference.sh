#!/bin/bash -l
#SBATCH
#SBATCH --job-name=submit_network_inference
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
	for i in $(seq 0.04 0.02 1.00);
	do
		sbatch run_glasso.sh GTEx $K weighted $i
		sbatch run_glasso.sh GTEx $K unweighted $i
		sbatch run_glasso.sh GTEx $K pooling_first $i
		sbatch run_glasso.sh GTEx $K correcting_first $i
	done
done

for K in $(seq 0 10 50);
do
	if [ $K == 0 ]
		then K=1
	fi
	if [ $K == 50 ]
                then K=49
        fi
	for i in $(seq 0.04 0.02 1.00);
	do
		sbatch run_glasso.sh GTEx $K weighted $i
		sbatch run_glasso.sh GTEx $K unweighted $i
		sbatch run_glasso.sh GTEx $K pooling_first $i
		sbatch run_glasso.sh GTEx $K correcting_first $i
	done
done

for K in $(seq 0 100 600);
do
	if [ $K == 0 ]
		then K=1
	fi
	if [ $K == 600 ]
                then K=566
        fi
	for i in $(seq 0.04 0.02 1.00);
	do
		sbatch run_glasso.sh sra_normal $K weighted $i
		sbatch run_glasso.sh sra_normal $K unweighted $i
		sbatch run_glasso.sh sra_normal $K pooling_first $i
		sbatch run_glasso.sh sra_normal $K correcting_first $i
	done
done

for i in $(seq 0.04 0.02 1.00);
do
  sbatch run_glasso.sh all_consensus 966 weighted $i
  sbatch run_glasso.sh normal_consensus 629 weighted $i
  sbatch run_glasso.sh cancer_consensus 386 weighted $i
  sbatch run_glasso.sh adipose/all 11 weighted $i
  sbatch run_glasso.sh adipose/GTEx 2 weighted $i
  sbatch run_glasso.sh airway/all 8 weighted $i
  sbatch run_glasso.sh B_cells/all 17 weighted $i
  sbatch run_glasso.sh B_cells/GTEx 1 weighted $i
  sbatch run_glasso.sh blood/all 65 weighted $i
  sbatch run_glasso.sh blood/GTEx 1 weighted $i
  sbatch run_glasso.sh breast/all 12 weighted $i
  sbatch run_glasso.sh breast/GTEx 1 weighted $i
  sbatch run_glasso.sh cardiac/all 13 weighted $i
  sbatch run_glasso.sh cardiac/GTEx 2 weighted $i
  sbatch run_glasso.sh central_nervous_system/all 53 weighted $i
  sbatch run_glasso.sh central_nervous_system/GTEx 13 weighted $i
  sbatch run_glasso.sh colon/all 10 weighted $i
  sbatch run_glasso.sh colon/GTEx 2 weighted $i
  sbatch run_glasso.sh esophagus/all 3 weighted $i
  sbatch run_glasso.sh esophagus/GTEx 3 weighted $i
  sbatch run_glasso.sh eye/all 2 weighted $i
  sbatch run_glasso.sh fibroblasts/all 13 weighted $i
  sbatch run_glasso.sh fibroblasts/GTEx 1 weighted $i
  sbatch run_glasso.sh hescs/all 21 weighted $i
  sbatch run_glasso.sh intestine/all 9 weighted $i
  sbatch run_glasso.sh intestine/GTEx 1 weighted $i
  sbatch run_glasso.sh ipscs/all 14 weighted $i
  sbatch run_glasso.sh kidney/all 26 weighted $i
  sbatch run_glasso.sh kidney/GTEx 1 weighted $i
  sbatch run_glasso.sh liver/all 28 weighted $i
  sbatch run_glasso.sh liver/GTEx 1 weighted $i
  sbatch run_glasso.sh lung/all 10 weighted $i
  sbatch run_glasso.sh lung/GTEx 1 weighted $i
  sbatch run_glasso.sh multipotent_cells/all 6 weighted $i
  sbatch run_glasso.sh myeloid_cells/all 29 weighted $i
  sbatch run_glasso.sh nervous_system/all 18 weighted $i
  sbatch run_glasso.sh nervous_system/GTEx 1 weighted $i
  sbatch run_glasso.sh pancreas/all 3 weighted $i
  sbatch run_glasso.sh pancreas/GTEx 1 weighted $i
  sbatch run_glasso.sh pbmcs_t_cells/all 58 weighted $i
  sbatch run_glasso.sh prostate/all 6 weighted $i
  sbatch run_glasso.sh prostate/GTEx 1 weighted $i
  sbatch run_glasso.sh skeletal_muscle/all 8 weighted $i
  sbatch run_glasso.sh skeletal_muscle/GTEx 1 weighted $i
  sbatch run_glasso.sh skin/all 20 weighted $i
  sbatch run_glasso.sh skin/GTEx 2 weighted $i
  sbatch run_glasso.sh stomach/all 4 weighted $i
  sbatch run_glasso.sh stomach/GTEx 1 weighted $i
  sbatch run_glasso.sh vascular/all 4 weighted $i
  sbatch run_glasso.sh vascular/GTEx 3 weighted $i
done

