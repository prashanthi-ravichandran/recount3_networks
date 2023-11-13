#!/bin/bash -l
#SBATCH
#SBATCH --job-name=submit_scale_free
#SBATCH --mail-type=END
#SBATCH --mail-user=pravich2@jhu.edu
#SBATCH --time=0:20:0
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

sbatch run_make_centrality.sh all_consensus 966 0.14
sbatch run_make_centrality.sh all_consensus 966 0.16
sbatch run_make_centrality.sh all_consensus 966 0.18
sbatch run_make_centrality.sh all_consensus 966 0.20

sbatch run_make_centrality.sh normal_consensus 629 0.14
sbatch run_make_centrality.sh normal_consensus 629 0.16
sbatch run_make_centrality.sh normal_consensus 629 0.18
sbatch run_make_centrality.sh normal_consensus 629 0.20

sbatch run_make_centrality.sh blood/all 65 0.18
sbatch run_make_centrality.sh blood/all 65 0.20
sbatch run_make_centrality.sh blood/all 65 0.22
sbatch run_make_centrality.sh blood/all 65 0.24
sbatch run_make_centrality.sh blood/all 65 0.26

sbatch run_make_centrality.sh blood/GTEx 1 0.24
sbatch run_make_centrality.sh blood/GTEx 1 0.26
sbatch run_make_centrality.sh blood/GTEx 1 0.28
sbatch run_make_centrality.sh blood/GTEx 1 0.30
sbatch run_make_centrality.sh blood/GTEx 1 0.32

sbatch run_make_centrality.sh central_nervous_system/all 53 0.20
sbatch run_make_centrality.sh central_nervous_system/all 53 0.22
sbatch run_make_centrality.sh central_nervous_system/all 53 0.24
sbatch run_make_centrality.sh central_nervous_system/all 53 0.26
sbatch run_make_centrality.sh central_nervous_system/all 53 0.28
sbatch run_make_centrality.sh central_nervous_system/all 53 0.30
sbatch run_make_centrality.sh central_nervous_system/all 53 0.32

sbatch run_make_centrality.sh central_nervous_system/GTEx 1 0.22
sbatch run_make_centrality.sh central_nervous_system/GTEx 1 0.24
sbatch run_make_centrality.sh central_nervous_system/GTEx 1 0.26
sbatch run_make_centrality.sh central_nervous_system/GTEx 1 0.28
sbatch run_make_centrality.sh central_nervous_system/GTEx 1 0.30
sbatch run_make_centrality.sh central_nervous_system/GTEx 1 0.32
sbatch run_make_centrality.sh central_nervous_system/GTEx 1 0.34
sbatch run_make_centrality.sh central_nervous_system/GTEx 1 0.36




