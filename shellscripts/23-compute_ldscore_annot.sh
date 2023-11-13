#!/bin/bash -l
#SBATCH
#SBATCH --job-name=submit_make_annot
#SBATCH --mail-type=END
#SBATCH --mail-user=pravich2@jhu.edu
#SBATCH --time=0:20:0
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

for chr in $(seq 1 1 22);
do
  sbatch run_compute_ld_score_annot.sh all_consensus_0.14 betweeness $chr
  sbatch run_compute_ld_score_annot.sh all_consensus_0.14 closeness $chr
  sbatch run_compute_ld_score_annot.sh all_consensus_0.14 degree $chr
  sbatch run_compute_ld_score_annot.sh all_consensus_0.14 eigen_centrality $chr
  sbatch run_compute_ld_score_annot.sh all_consensus_0.14 maximum_weight $chr
  sbatch run_compute_ld_score_annot.sh all_consensus_0.14 strength $chr
  sbatch run_compute_ld_score_annot.sh all_consensus_0.14 page_rank $chr
  
  sbatch run_compute_ld_score_annot.sh all_consensus_0.16 betweeness $chr
  sbatch run_compute_ld_score_annot.sh all_consensus_0.16 closeness $chr
  sbatch run_compute_ld_score_annot.sh all_consensus_0.16 degree $chr
  sbatch run_compute_ld_score_annot.sh all_consensus_0.16 eigen_centrality $chr
  sbatch run_compute_ld_score_annot.sh all_consensus_0.16 maximum_weight $chr
  sbatch run_compute_ld_score_annot.sh all_consensus_0.16 strength $chr
  sbatch run_compute_ld_score_annot.sh all_consensus_0.16 page_rank $chr
  
  sbatch run_compute_ld_score_annot.sh all_consensus_0.18 betweeness $chr
  sbatch run_compute_ld_score_annot.sh all_consensus_0.18 closeness $chr
  sbatch run_compute_ld_score_annot.sh all_consensus_0.18 degree $chr
  sbatch run_compute_ld_score_annot.sh all_consensus_0.18 eigen_centrality $chr
  sbatch run_compute_ld_score_annot.sh all_consensus_0.18 maximum_weight $chr
  sbatch run_compute_ld_score_annot.sh all_consensus_0.18 strength $chr
  sbatch run_compute_ld_score_annot.sh all_consensus_0.18 page_rank $chr
  
  sbatch run_compute_ld_score_annot.sh all_consensus_0.20 betweeness $chr
  sbatch run_compute_ld_score_annot.sh all_consensus_0.20 closeness $chr
  sbatch run_compute_ld_score_annot.sh all_consensus_0.20 degree $chr
  sbatch run_compute_ld_score_annot.sh all_consensus_0.20 eigen_centrality $chr
  sbatch run_compute_ld_score_annot.sh all_consensus_0.20 maximum_weight $chr
  sbatch run_compute_ld_score_annot.sh all_consensus_0.20 strength $chr
  sbatch run_compute_ld_score_annot.sh all_consensus_0.20 page_rank $chr
  
  sbatch run_compute_ld_score_annot.sh normal_consensus_0.14 betweeness $chr
  sbatch run_compute_ld_score_annot.sh normal_consensus_0.14 closeness $chr
  sbatch run_compute_ld_score_annot.sh normal_consensus_0.14 degree $chr
  sbatch run_compute_ld_score_annot.sh normal_consensus_0.14 eigen_centrality $chr
  sbatch run_compute_ld_score_annot.sh normal_consensus_0.14 maximum_weight $chr
  sbatch run_compute_ld_score_annot.sh normal_consensus_0.14 strength $chr
  sbatch run_compute_ld_score_annot.sh normal_consensus_0.14 page_rank $chr
  
  sbatch run_compute_ld_score_annot.sh normal_consensus_0.16 betweeness $chr
  sbatch run_compute_ld_score_annot.sh normal_consensus_0.16 closeness $chr
  sbatch run_compute_ld_score_annot.sh normal_consensus_0.16 degree $chr
  sbatch run_compute_ld_score_annot.sh normal_consensus_0.16 eigen_centrality $chr
  sbatch run_compute_ld_score_annot.sh normal_consensus_0.16 maximum_weight $chr
  sbatch run_compute_ld_score_annot.sh normal_consensus_0.16 strength $chr
  sbatch run_compute_ld_score_annot.sh normal_consensus_0.16 page_rank $chr
  
  sbatch run_compute_ld_score_annot.sh normal_consensus_0.18 betweeness $chr
  sbatch run_compute_ld_score_annot.sh normal_consensus_0.18 closeness $chr
  sbatch run_compute_ld_score_annot.sh normal_consensus_0.18 degree $chr
  sbatch run_compute_ld_score_annot.sh normal_consensus_0.18 eigen_centrality $chr
  sbatch run_compute_ld_score_annot.sh normal_consensus_0.18 maximum_weight $chr
  sbatch run_compute_ld_score_annot.sh normal_consensus_0.18 strength $chr
  sbatch run_compute_ld_score_annot.sh normal_consensus_0.18 page_rank $chr
  
  sbatch run_compute_ld_score_annot.sh normal_consensus_0.20 betweeness $chr
  sbatch run_compute_ld_score_annot.sh normal_consensus_0.20 closeness $chr
  sbatch run_compute_ld_score_annot.sh normal_consensus_0.20 degree $chr
  sbatch run_compute_ld_score_annot.sh normal_consensus_0.20 eigen_centrality $chr
  sbatch run_compute_ld_score_annot.sh normal_consensus_0.20 maximum_weight $chr
  sbatch run_compute_ld_score_annot.sh normal_consensus_0.20 strength $chr
  sbatch run_compute_ld_score_annot.sh normal_consensus_0.20 page_rank $chr
  
  sbatch run_compute_ld_score_annot.sh blood_consensus_0.18 betweeness $chr
  sbatch run_compute_ld_score_annot.sh blood_consensus_0.18 closeness $chr
  sbatch run_compute_ld_score_annot.sh blood_consensus_0.18 degree $chr
  sbatch run_compute_ld_score_annot.sh blood_consensus_0.18 eigen_centrality $chr
  sbatch run_compute_ld_score_annot.sh blood_consensus_0.18 maximum_weight $chr
  sbatch run_compute_ld_score_annot.sh blood_consensus_0.18 strength $chr
  sbatch run_compute_ld_score_annot.sh blood_consensus_0.18 page_rank $chr

  sbatch run_compute_ld_score_annot.sh blood_consensus_0.20 betweeness $chr
  sbatch run_compute_ld_score_annot.sh blood_consensus_0.20 closeness $chr
  sbatch run_compute_ld_score_annot.sh blood_consensus_0.20 degree $chr
  sbatch run_compute_ld_score_annot.sh blood_consensus_0.20 eigen_centrality $chr
  sbatch run_compute_ld_score_annot.sh blood_consensus_0.20 maximum_weight $chr
  sbatch run_compute_ld_score_annot.sh blood_consensus_0.20 strength $chr
  sbatch run_compute_ld_score_annot.sh blood_consensus_0.20 page_rank $chr

  sbatch run_compute_ld_score_annot.sh blood_consensus_0.22 betweeness $chr
  sbatch run_compute_ld_score_annot.sh blood_consensus_0.22 closeness $chr
  sbatch run_compute_ld_score_annot.sh blood_consensus_0.22 degree $chr
  sbatch run_compute_ld_score_annot.sh blood_consensus_0.22 eigen_centrality $chr
  sbatch run_compute_ld_score_annot.sh blood_consensus_0.22 maximum_weight $chr
  sbatch run_compute_ld_score_annot.sh blood_consensus_0.22 strength $chr
  sbatch run_compute_ld_score_annot.sh blood_consensus_0.22 page_rank $chr
  
  sbatch run_compute_ld_score_annot.sh blood_consensus_0.24 betweeness $chr
  sbatch run_compute_ld_score_annot.sh blood_consensus_0.24 closeness $chr
  sbatch run_compute_ld_score_annot.sh blood_consensus_0.24 degree $chr
  sbatch run_compute_ld_score_annot.sh blood_consensus_0.24 eigen_centrality $chr
  sbatch run_compute_ld_score_annot.sh blood_consensus_0.24 maximum_weight $chr
  sbatch run_compute_ld_score_annot.sh blood_consensus_0.24 strength $chr
  sbatch run_compute_ld_score_annot.sh blood_consensus_0.24 page_rank $chr
  
  sbatch run_compute_ld_score_annot.sh blood_consensus_0.26 betweeness $chr
  sbatch run_compute_ld_score_annot.sh blood_consensus_0.26 closeness $chr
  sbatch run_compute_ld_score_annot.sh blood_consensus_0.26 degree $chr
  sbatch run_compute_ld_score_annot.sh blood_consensus_0.26 eigen_centrality $chr
  sbatch run_compute_ld_score_annot.sh blood_consensus_0.26 maximum_weight $chr
  sbatch run_compute_ld_score_annot.sh blood_consensus_0.26 strength $chr
  sbatch run_compute_ld_score_annot.sh blood_consensus_0.26 page_rank $chr
  
  sbatch run_compute_ld_score_annot.sh blood_GTEx_0.24 betweeness $chr
  sbatch run_compute_ld_score_annot.sh blood_GTEx_0.24 closeness $chr
  sbatch run_compute_ld_score_annot.sh blood_GTEx_0.24 degree $chr
  sbatch run_compute_ld_score_annot.sh blood_GTEx_0.24 eigen_centrality $chr
  sbatch run_compute_ld_score_annot.sh blood_GTEx_0.24 maximum_weight $chr
  sbatch run_compute_ld_score_annot.sh blood_GTEx_0.24 strength $chr
  sbatch run_compute_ld_score_annot.sh blood_GTEx_0.24 page_rank $chr

  sbatch run_compute_ld_score_annot.sh blood_GTEx_0.26 betweeness $chr
  sbatch run_compute_ld_score_annot.sh blood_GTEx_0.26 closeness $chr
  sbatch run_compute_ld_score_annot.sh blood_GTEx_0.26 degree $chr
  sbatch run_compute_ld_score_annot.sh blood_GTEx_0.26 eigen_centrality $chr
  sbatch run_compute_ld_score_annot.sh blood_GTEx_0.26 maximum_weight $chr
  sbatch run_compute_ld_score_annot.sh blood_GTEx_0.26 strength $chr
  sbatch run_compute_ld_score_annot.sh blood_GTEx_0.26 page_rank $chr
  
  sbatch run_compute_ld_score_annot.sh blood_GTEx_0.28 betweeness $chr
  sbatch run_compute_ld_score_annot.sh blood_GTEx_0.28 closeness $chr
  sbatch run_compute_ld_score_annot.sh blood_GTEx_0.28 degree $chr
  sbatch run_compute_ld_score_annot.sh blood_GTEx_0.28 eigen_centrality $chr
  sbatch run_compute_ld_score_annot.sh blood_GTEx_0.28 maximum_weight $chr
  sbatch run_compute_ld_score_annot.sh blood_GTEx_0.28 strength $chr
  sbatch run_compute_ld_score_annot.sh blood_GTEx_0.28 page_rank $chr

  sbatch run_compute_ld_score_annot.sh blood_GTEx_0.30 betweeness $chr
  sbatch run_compute_ld_score_annot.sh blood_GTEx_0.30 closeness $chr
  sbatch run_compute_ld_score_annot.sh blood_GTEx_0.30 degree $chr
  sbatch run_compute_ld_score_annot.sh blood_GTEx_0.30 eigen_centrality $chr
  sbatch run_compute_ld_score_annot.sh blood_GTEx_0.30 maximum_weight $chr
  sbatch run_compute_ld_score_annot.sh blood_GTEx_0.30 strength $chr
  sbatch run_compute_ld_score_annot.sh blood_GTEx_0.30 page_rank $chr
  
  sbatch run_compute_ld_score_annot.sh blood_GTEx_0.32 betweeness $chr
  sbatch run_compute_ld_score_annot.sh blood_GTEx_0.32 closeness $chr
  sbatch run_compute_ld_score_annot.sh blood_GTEx_0.32 degree $chr
  sbatch run_compute_ld_score_annot.sh blood_GTEx_0.32 eigen_centrality $chr
  sbatch run_compute_ld_score_annot.sh blood_GTEx_0.32 maximum_weight $chr
  sbatch run_compute_ld_score_annot.sh blood_GTEx_0.32 strength $chr
  sbatch run_compute_ld_score_annot.sh blood_GTEx_0.32 page_rank $chr
  
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.20 betweeness $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.20 closeness $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.20 degree $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.20 eigen_centrality $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.20 maximum_weight $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.20 strength $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.20 page_rank $chr

  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.22 betweeness $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.22 closeness $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.22 degree $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.22 eigen_centrality $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.22 maximum_weight $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.22 strength $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.22 page_rank $chr

  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.24 betweeness $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.24 closeness $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.24 degree $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.24 eigen_centrality $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.24 maximum_weight $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.24 strength $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.24 page_rank $chr

  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.26 betweeness $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.26 closeness $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.26 degree $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.26 eigen_centrality $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.26 maximum_weight $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.26 strength $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.26 page_rank $chr

  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.28 betweeness $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.28 closeness $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.28 degree $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.28 eigen_centrality $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.28 maximum_weight $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.28 strength $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.28 page_rank $chr

  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.30 betweeness $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.30 closeness $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.30 degree $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.30 eigen_centrality $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.30 maximum_weight $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.30 strength $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.30 page_rank $chr

  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.32 betweeness $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.32 closeness $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.32 degree $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.32 eigen_centrality $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.32 maximum_weight $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.32 strength $chr
  sbatch run_compute_ld_score_annot.sh CNS_consensus_0.32 page_rank $chr

  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.22 betweeness $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.22 closeness $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.22 degree $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.22 eigen_centrality $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.22 maximum_weight $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.22 strength $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.22 page_rank $chr

  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.24 betweeness $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.24 closeness $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.24 degree $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.24 eigen_centrality $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.24 maximum_weight $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.24 strength $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.24 page_rank $chr

  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.26 betweeness $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.26 closeness $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.26 degree $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.26 eigen_centrality $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.26 maximum_weight $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.26 strength $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.26 page_rank $chr

  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.28 betweeness $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.28 closeness $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.28 degree $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.28 eigen_centrality $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.28 maximum_weight $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.28 strength $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.28 page_rank $chr

  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.30 betweeness $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.30 closeness $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.30 degree $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.30 eigen_centrality $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.30 maximum_weight $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.30 strength $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.30 page_rank $chr

  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.32 betweeness $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.32 closeness $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.32 degree $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.32 eigen_centrality $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.32 maximum_weight $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.32 strength $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.32 page_rank $chr
  
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.34 betweeness $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.34 closeness $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.34 degree $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.34 eigen_centrality $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.34 maximum_weight $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.34 strength $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.34 page_rank $chr

  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.36 betweeness $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.36 closeness $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.36 degree $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.36 eigen_centrality $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.36 maximum_weight $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.36 strength $chr
  sbatch run_compute_ld_score_annot.sh CNS_GTEx_0.36 page_rank $chr
done
