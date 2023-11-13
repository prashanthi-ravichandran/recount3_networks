#!/bin/bash -l
#SBATCH
#SBATCH --job-name=submit_make_annot
#SBATCH --mail-type=END
#SBATCH --mail-user=pravich2@jhu.edu
#SBATCH --time=0:20:0
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

sbatch run_compute_annot_sd.sh all_consensus_0.14 betweeness all_genes
sbatch run_compute_annot_sd.sh all_consensus_0.14 closeness all_genes
sbatch run_compute_annot_sd.sh all_consensus_0.14 degree all_genes
sbatch run_compute_annot_sd.sh all_consensus_0.14 eigen_centrality all_genes
sbatch run_compute_annot_sd.sh all_consensus_0.14 maximum_weight all_genes
sbatch run_compute_annot_sd.sh all_consensus_0.14 strength all_genes
sbatch run_compute_annot_sd.sh all_consensus_0.14 page_rank all_genes

sbatch run_compute_annot_sd.sh all_consensus_0.16 betweeness all_genes
sbatch run_compute_annot_sd.sh all_consensus_0.16 closeness all_genes
sbatch run_compute_annot_sd.sh all_consensus_0.16 degree all_genes
sbatch run_compute_annot_sd.sh all_consensus_0.16 eigen_centrality all_genes
sbatch run_compute_annot_sd.sh all_consensus_0.16 maximum_weight all_genes
sbatch run_compute_annot_sd.sh all_consensus_0.16 strength all_genes
sbatch run_compute_annot_sd.sh all_consensus_0.16 page_rank all_genes

sbatch run_compute_annot_sd.sh all_consensus_0.18 betweeness all_genes
sbatch run_compute_annot_sd.sh all_consensus_0.18 closeness all_genes
sbatch run_compute_annot_sd.sh all_consensus_0.18 degree all_genes
sbatch run_compute_annot_sd.sh all_consensus_0.18 eigen_centrality all_genes
sbatch run_compute_annot_sd.sh all_consensus_0.18 maximum_weight all_genes
sbatch run_compute_annot_sd.sh all_consensus_0.18 strength all_genes
sbatch run_compute_annot_sd.sh all_consensus_0.18 page_rank all_genes

sbatch run_compute_annot_sd.sh all_consensus_0.20 betweeness all_genes
sbatch run_compute_annot_sd.sh all_consensus_0.20 closeness all_genes
sbatch run_compute_annot_sd.sh all_consensus_0.20 degree all_genes
sbatch run_compute_annot_sd.sh all_consensus_0.20 eigen_centrality all_genes
sbatch run_compute_annot_sd.sh all_consensus_0.20 maximum_weight all_genes
sbatch run_compute_annot_sd.sh all_consensus_0.20 strength all_genes
sbatch run_compute_annot_sd.sh all_consensus_0.20 page_rank all_genes

sbatch run_compute_annot_sd.sh normal_consensus_0.14 betweeness all_genes
sbatch run_compute_annot_sd.sh normal_consensus_0.14 closeness all_genes
sbatch run_compute_annot_sd.sh normal_consensus_0.14 degree all_genes
sbatch run_compute_annot_sd.sh normal_consensus_0.14 eigen_centrality all_genes
sbatch run_compute_annot_sd.sh normal_consensus_0.14 maximum_weight all_genes
sbatch run_compute_annot_sd.sh normal_consensus_0.14 strength all_genes
sbatch run_compute_annot_sd.sh normal_consensus_0.14 page_rank all_genes

sbatch run_compute_annot_sd.sh normal_consensus_0.16 betweeness all_genes
sbatch run_compute_annot_sd.sh normal_consensus_0.16 closeness all_genes
sbatch run_compute_annot_sd.sh normal_consensus_0.16 degree all_genes
sbatch run_compute_annot_sd.sh normal_consensus_0.16 eigen_centrality all_genes
sbatch run_compute_annot_sd.sh normal_consensus_0.16 maximum_weight all_genes
sbatch run_compute_annot_sd.sh normal_consensus_0.16 strength all_genes
sbatch run_compute_annot_sd.sh normal_consensus_0.16 page_rank all_genes

sbatch run_compute_annot_sd.sh normal_consensus_0.18 betweeness all_genes
sbatch run_compute_annot_sd.sh normal_consensus_0.18 closeness all_genes
sbatch run_compute_annot_sd.sh normal_consensus_0.18 degree all_genes
sbatch run_compute_annot_sd.sh normal_consensus_0.18 eigen_centrality all_genes
sbatch run_compute_annot_sd.sh normal_consensus_0.18 maximum_weight all_genes
sbatch run_compute_annot_sd.sh normal_consensus_0.18 strength all_genes
sbatch run_compute_annot_sd.sh normal_consensus_0.18 page_rank all_genes

sbatch run_compute_annot_sd.sh normal_consensus_0.20 betweeness all_genes
sbatch run_compute_annot_sd.sh normal_consensus_0.20 closeness all_genes
sbatch run_compute_annot_sd.sh normal_consensus_0.20 degree all_genes
sbatch run_compute_annot_sd.sh normal_consensus_0.20 eigen_centrality all_genes
sbatch run_compute_annot_sd.sh normal_consensus_0.20 maximum_weight all_genes
sbatch run_compute_annot_sd.sh normal_consensus_0.20 strength all_genes
sbatch run_compute_annot_sd.sh normal_consensus_0.20 page_rank all_genes

sbatch run_compute_annot_sd.sh blood_consensus_0.18 betweeness all_genes
sbatch run_compute_annot_sd.sh blood_consensus_0.18 closeness all_genes
sbatch run_compute_annot_sd.sh blood_consensus_0.18 degree all_genes
sbatch run_compute_annot_sd.sh blood_consensus_0.18 eigen_centrality all_genes
sbatch run_compute_annot_sd.sh blood_consensus_0.18 maximum_weight all_genes
sbatch run_compute_annot_sd.sh blood_consensus_0.18 strength all_genes
sbatch run_compute_annot_sd.sh blood_consensus_0.18 page_rank all_genes

sbatch run_compute_annot_sd.sh blood_consensus_0.20 betweeness all_genes
sbatch run_compute_annot_sd.sh blood_consensus_0.20 closeness all_genes
sbatch run_compute_annot_sd.sh blood_consensus_0.20 degree all_genes
sbatch run_compute_annot_sd.sh blood_consensus_0.20 eigen_centrality all_genes
sbatch run_compute_annot_sd.sh blood_consensus_0.20 maximum_weight all_genes
sbatch run_compute_annot_sd.sh blood_consensus_0.20 strength all_genes
sbatch run_compute_annot_sd.sh blood_consensus_0.20 page_rank all_genes

sbatch run_compute_annot_sd.sh blood_consensus_0.22 betweeness all_genes
sbatch run_compute_annot_sd.sh blood_consensus_0.22 closeness all_genes
sbatch run_compute_annot_sd.sh blood_consensus_0.22 degree all_genes
sbatch run_compute_annot_sd.sh blood_consensus_0.22 eigen_centrality all_genes
sbatch run_compute_annot_sd.sh blood_consensus_0.22 maximum_weight all_genes
sbatch run_compute_annot_sd.sh blood_consensus_0.22 strength all_genes
sbatch run_compute_annot_sd.sh blood_consensus_0.22 page_rank all_genes

sbatch run_compute_annot_sd.sh blood_consensus_0.24 betweeness all_genes
sbatch run_compute_annot_sd.sh blood_consensus_0.24 closeness all_genes
sbatch run_compute_annot_sd.sh blood_consensus_0.24 degree all_genes
sbatch run_compute_annot_sd.sh blood_consensus_0.24 eigen_centrality all_genes
sbatch run_compute_annot_sd.sh blood_consensus_0.24 maximum_weight all_genes
sbatch run_compute_annot_sd.sh blood_consensus_0.24 strength all_genes
sbatch run_compute_annot_sd.sh blood_consensus_0.24 page_rank all_genes

sbatch run_compute_annot_sd.sh blood_consensus_0.26 betweeness all_genes
sbatch run_compute_annot_sd.sh blood_consensus_0.26 closeness all_genes
sbatch run_compute_annot_sd.sh blood_consensus_0.26 degree all_genes
sbatch run_compute_annot_sd.sh blood_consensus_0.26 eigen_centrality all_genes
sbatch run_compute_annot_sd.sh blood_consensus_0.26 maximum_weight all_genes
sbatch run_compute_annot_sd.sh blood_consensus_0.26 strength all_genes
sbatch run_compute_annot_sd.sh blood_consensus_0.26 page_rank all_genes

sbatch run_compute_annot_sd.sh blood_GTEx_0.24 betweeness all_genes
sbatch run_compute_annot_sd.sh blood_GTEx_0.24 closeness all_genes
sbatch run_compute_annot_sd.sh blood_GTEx_0.24 degree all_genes
sbatch run_compute_annot_sd.sh blood_GTEx_0.24 eigen_centrality all_genes
sbatch run_compute_annot_sd.sh blood_GTEx_0.24 maximum_weight all_genes
sbatch run_compute_annot_sd.sh blood_GTEx_0.24 strength all_genes
sbatch run_compute_annot_sd.sh blood_GTEx_0.24 page_rank all_genes

sbatch run_compute_annot_sd.sh blood_GTEx_0.26 betweeness all_genes
sbatch run_compute_annot_sd.sh blood_GTEx_0.26 closeness all_genes
sbatch run_compute_annot_sd.sh blood_GTEx_0.26 degree all_genes
sbatch run_compute_annot_sd.sh blood_GTEx_0.26 eigen_centrality all_genes
sbatch run_compute_annot_sd.sh blood_GTEx_0.26 maximum_weight all_genes
sbatch run_compute_annot_sd.sh blood_GTEx_0.26 strength all_genes
sbatch run_compute_annot_sd.sh blood_GTEx_0.26 page_rank all_genes

sbatch run_compute_annot_sd.sh blood_GTEx_0.28 betweeness all_genes
sbatch run_compute_annot_sd.sh blood_GTEx_0.28 closeness all_genes
sbatch run_compute_annot_sd.sh blood_GTEx_0.28 degree all_genes
sbatch run_compute_annot_sd.sh blood_GTEx_0.28 eigen_centrality all_genes
sbatch run_compute_annot_sd.sh blood_GTEx_0.28 maximum_weight all_genes
sbatch run_compute_annot_sd.sh blood_GTEx_0.28 strength all_genes
sbatch run_compute_annot_sd.sh blood_GTEx_0.28 page_rank all_genes

sbatch run_compute_annot_sd.sh blood_GTEx_0.30 betweeness all_genes
sbatch run_compute_annot_sd.sh blood_GTEx_0.30 closeness all_genes
sbatch run_compute_annot_sd.sh blood_GTEx_0.30 degree all_genes
sbatch run_compute_annot_sd.sh blood_GTEx_0.30 eigen_centrality all_genes
sbatch run_compute_annot_sd.sh blood_GTEx_0.30 maximum_weight all_genes
sbatch run_compute_annot_sd.sh blood_GTEx_0.30 strength all_genes
sbatch run_compute_annot_sd.sh blood_GTEx_0.30 page_rank all_genes

sbatch run_compute_annot_sd.sh blood_GTEx_0.32 betweeness all_genes
sbatch run_compute_annot_sd.sh blood_GTEx_0.32 closeness all_genes
sbatch run_compute_annot_sd.sh blood_GTEx_0.32 degree all_genes
sbatch run_compute_annot_sd.sh blood_GTEx_0.32 eigen_centrality all_genes
sbatch run_compute_annot_sd.sh blood_GTEx_0.32 maximum_weight all_genes
sbatch run_compute_annot_sd.sh blood_GTEx_0.32 strength all_genes
sbatch run_compute_annot_sd.sh blood_GTEx_0.32 page_rank all_genes

sbatch run_compute_annot_sd.sh CNS_consensus_0.20 betweeness all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.20 closeness all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.20 degree all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.20 eigen_centrality all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.20 maximum_weight all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.20 strength all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.20 page_rank all_genes

sbatch run_compute_annot_sd.sh CNS_consensus_0.22 betweeness all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.22 closeness all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.22 degree all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.22 eigen_centrality all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.22 maximum_weight all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.22 strength all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.22 page_rank all_genes

sbatch run_compute_annot_sd.sh CNS_consensus_0.24 betweeness all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.24 closeness all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.24 degree all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.24 eigen_centrality all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.24 maximum_weight all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.24 strength all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.24 page_rank all_genes

sbatch run_compute_annot_sd.sh CNS_consensus_0.26 betweeness all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.26 closeness all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.26 degree all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.26 eigen_centrality all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.26 maximum_weight all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.26 strength all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.26 page_rank all_genes

sbatch run_compute_annot_sd.sh CNS_consensus_0.28 betweeness all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.28 closeness all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.28 degree all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.28 eigen_centrality all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.28 maximum_weight all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.28 strength all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.28 page_rank all_genes

sbatch run_compute_annot_sd.sh CNS_consensus_0.30 betweeness all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.30 closeness all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.30 degree all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.30 eigen_centrality all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.30 maximum_weight all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.30 strength all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.30 page_rank all_genes

sbatch run_compute_annot_sd.sh CNS_consensus_0.32 betweeness all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.32 closeness all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.32 degree all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.32 eigen_centrality all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.32 maximum_weight all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.32 strength all_genes
sbatch run_compute_annot_sd.sh CNS_consensus_0.32 page_rank all_genes

sbatch run_compute_annot_sd.sh CNS_GTEx_0.22 betweeness all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.22 closeness all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.22 degree all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.22 eigen_centrality all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.22 maximum_weight all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.22 strength all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.22 page_rank all_genes

sbatch run_compute_annot_sd.sh CNS_GTEx_0.24 betweeness all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.24 closeness all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.24 degree all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.24 eigen_centrality all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.24 maximum_weight all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.24 strength all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.24 page_rank all_genes

sbatch run_compute_annot_sd.sh CNS_GTEx_0.26 betweeness all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.26 closeness all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.26 degree all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.26 eigen_centrality all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.26 maximum_weight all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.26 strength all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.26 page_rank all_genes

sbatch run_compute_annot_sd.sh CNS_GTEx_0.28 betweeness all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.28 closeness all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.28 degree all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.28 eigen_centrality all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.28 maximum_weight all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.28 strength all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.28 page_rank all_genes

sbatch run_compute_annot_sd.sh CNS_GTEx_0.30 betweeness all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.30 closeness all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.30 degree all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.30 eigen_centrality all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.30 maximum_weight all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.30 strength all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.30 page_rank all_genes

sbatch run_compute_annot_sd.sh CNS_GTEx_0.32 betweeness all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.32 closeness all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.32 degree all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.32 eigen_centrality all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.32 maximum_weight all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.32 strength all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.32 page_rank all_genes

sbatch run_compute_annot_sd.sh CNS_GTEx_0.34 betweeness all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.34 closeness all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.34 degree all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.34 eigen_centrality all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.34 maximum_weight all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.34 strength all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.34 page_rank all_genes

sbatch run_compute_annot_sd.sh CNS_GTEx_0.36 betweeness all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.36 closeness all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.36 degree all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.36 eigen_centrality all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.36 maximum_weight all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.36 strength all_genes
sbatch run_compute_annot_sd.sh CNS_GTEx_0.36 page_rank all_genes

sbatch run_compute_annot_sd.sh all_consensus_0.14 betweeness baseline
sbatch run_compute_annot_sd.sh all_consensus_0.14 closeness baseline
sbatch run_compute_annot_sd.sh all_consensus_0.14 degree baseline
sbatch run_compute_annot_sd.sh all_consensus_0.14 eigen_centrality baseline
sbatch run_compute_annot_sd.sh all_consensus_0.14 maximum_weight baseline
sbatch run_compute_annot_sd.sh all_consensus_0.14 strength baseline
sbatch run_compute_annot_sd.sh all_consensus_0.14 page_rank baseline

sbatch run_compute_annot_sd.sh all_consensus_0.16 betweeness baseline
sbatch run_compute_annot_sd.sh all_consensus_0.16 closeness baseline
sbatch run_compute_annot_sd.sh all_consensus_0.16 degree baseline
sbatch run_compute_annot_sd.sh all_consensus_0.16 eigen_centrality baseline
sbatch run_compute_annot_sd.sh all_consensus_0.16 maximum_weight baseline
sbatch run_compute_annot_sd.sh all_consensus_0.16 strength baseline
sbatch run_compute_annot_sd.sh all_consensus_0.16 page_rank baseline

sbatch run_compute_annot_sd.sh all_consensus_0.18 betweeness baseline
sbatch run_compute_annot_sd.sh all_consensus_0.18 closeness baseline
sbatch run_compute_annot_sd.sh all_consensus_0.18 degree baseline
sbatch run_compute_annot_sd.sh all_consensus_0.18 eigen_centrality baseline
sbatch run_compute_annot_sd.sh all_consensus_0.18 maximum_weight baseline
sbatch run_compute_annot_sd.sh all_consensus_0.18 strength baseline
sbatch run_compute_annot_sd.sh all_consensus_0.18 page_rank baseline

sbatch run_compute_annot_sd.sh all_consensus_0.20 betweeness baseline
sbatch run_compute_annot_sd.sh all_consensus_0.20 closeness baseline
sbatch run_compute_annot_sd.sh all_consensus_0.20 degree baseline
sbatch run_compute_annot_sd.sh all_consensus_0.20 eigen_centrality baseline
sbatch run_compute_annot_sd.sh all_consensus_0.20 maximum_weight baseline
sbatch run_compute_annot_sd.sh all_consensus_0.20 strength baseline
sbatch run_compute_annot_sd.sh all_consensus_0.20 page_rank baseline

sbatch run_compute_annot_sd.sh normal_consensus_0.14 betweeness baseline
sbatch run_compute_annot_sd.sh normal_consensus_0.14 closeness baseline
sbatch run_compute_annot_sd.sh normal_consensus_0.14 degree baseline
sbatch run_compute_annot_sd.sh normal_consensus_0.14 eigen_centrality baseline
sbatch run_compute_annot_sd.sh normal_consensus_0.14 maximum_weight baseline
sbatch run_compute_annot_sd.sh normal_consensus_0.14 strength baseline
sbatch run_compute_annot_sd.sh normal_consensus_0.14 page_rank baseline

sbatch run_compute_annot_sd.sh normal_consensus_0.16 betweeness baseline
sbatch run_compute_annot_sd.sh normal_consensus_0.16 closeness baseline
sbatch run_compute_annot_sd.sh normal_consensus_0.16 degree baseline
sbatch run_compute_annot_sd.sh normal_consensus_0.16 eigen_centrality baseline
sbatch run_compute_annot_sd.sh normal_consensus_0.16 maximum_weight baseline
sbatch run_compute_annot_sd.sh normal_consensus_0.16 strength baseline
sbatch run_compute_annot_sd.sh normal_consensus_0.16 page_rank baseline

sbatch run_compute_annot_sd.sh normal_consensus_0.18 betweeness baseline
sbatch run_compute_annot_sd.sh normal_consensus_0.18 closeness baseline
sbatch run_compute_annot_sd.sh normal_consensus_0.18 degree baseline
sbatch run_compute_annot_sd.sh normal_consensus_0.18 eigen_centrality baseline
sbatch run_compute_annot_sd.sh normal_consensus_0.18 maximum_weight baseline
sbatch run_compute_annot_sd.sh normal_consensus_0.18 strength baseline
sbatch run_compute_annot_sd.sh normal_consensus_0.18 page_rank baseline

sbatch run_compute_annot_sd.sh normal_consensus_0.20 betweeness baseline
sbatch run_compute_annot_sd.sh normal_consensus_0.20 closeness baseline
sbatch run_compute_annot_sd.sh normal_consensus_0.20 degree baseline
sbatch run_compute_annot_sd.sh normal_consensus_0.20 eigen_centrality baseline
sbatch run_compute_annot_sd.sh normal_consensus_0.20 maximum_weight baseline
sbatch run_compute_annot_sd.sh normal_consensus_0.20 strength baseline
sbatch run_compute_annot_sd.sh normal_consensus_0.20 page_rank baseline

sbatch run_compute_annot_sd.sh blood_consensus_0.18 betweeness baseline
sbatch run_compute_annot_sd.sh blood_consensus_0.18 closeness baseline
sbatch run_compute_annot_sd.sh blood_consensus_0.18 degree baseline
sbatch run_compute_annot_sd.sh blood_consensus_0.18 eigen_centrality baseline
sbatch run_compute_annot_sd.sh blood_consensus_0.18 maximum_weight baseline
sbatch run_compute_annot_sd.sh blood_consensus_0.18 strength baseline
sbatch run_compute_annot_sd.sh blood_consensus_0.18 page_rank baseline

sbatch run_compute_annot_sd.sh blood_consensus_0.20 betweeness baseline
sbatch run_compute_annot_sd.sh blood_consensus_0.20 closeness baseline
sbatch run_compute_annot_sd.sh blood_consensus_0.20 degree baseline
sbatch run_compute_annot_sd.sh blood_consensus_0.20 eigen_centrality baseline
sbatch run_compute_annot_sd.sh blood_consensus_0.20 maximum_weight baseline
sbatch run_compute_annot_sd.sh blood_consensus_0.20 strength baseline
sbatch run_compute_annot_sd.sh blood_consensus_0.20 page_rank baseline

sbatch run_compute_annot_sd.sh blood_consensus_0.22 betweeness baseline
sbatch run_compute_annot_sd.sh blood_consensus_0.22 closeness baseline
sbatch run_compute_annot_sd.sh blood_consensus_0.22 degree baseline
sbatch run_compute_annot_sd.sh blood_consensus_0.22 eigen_centrality baseline
sbatch run_compute_annot_sd.sh blood_consensus_0.22 maximum_weight baseline
sbatch run_compute_annot_sd.sh blood_consensus_0.22 strength baseline
sbatch run_compute_annot_sd.sh blood_consensus_0.22 page_rank baseline

sbatch run_compute_annot_sd.sh blood_consensus_0.24 betweeness baseline
sbatch run_compute_annot_sd.sh blood_consensus_0.24 closeness baseline
sbatch run_compute_annot_sd.sh blood_consensus_0.24 degree baseline
sbatch run_compute_annot_sd.sh blood_consensus_0.24 eigen_centrality baseline
sbatch run_compute_annot_sd.sh blood_consensus_0.24 maximum_weight baseline
sbatch run_compute_annot_sd.sh blood_consensus_0.24 strength baseline
sbatch run_compute_annot_sd.sh blood_consensus_0.24 page_rank baseline

sbatch run_compute_annot_sd.sh blood_consensus_0.26 betweeness baseline
sbatch run_compute_annot_sd.sh blood_consensus_0.26 closeness baseline
sbatch run_compute_annot_sd.sh blood_consensus_0.26 degree baseline
sbatch run_compute_annot_sd.sh blood_consensus_0.26 eigen_centrality baseline
sbatch run_compute_annot_sd.sh blood_consensus_0.26 maximum_weight baseline
sbatch run_compute_annot_sd.sh blood_consensus_0.26 strength baseline
sbatch run_compute_annot_sd.sh blood_consensus_0.26 page_rank baseline

sbatch run_compute_annot_sd.sh blood_GTEx_0.24 betweeness baseline
sbatch run_compute_annot_sd.sh blood_GTEx_0.24 closeness baseline
sbatch run_compute_annot_sd.sh blood_GTEx_0.24 degree baseline
sbatch run_compute_annot_sd.sh blood_GTEx_0.24 eigen_centrality baseline
sbatch run_compute_annot_sd.sh blood_GTEx_0.24 maximum_weight baseline
sbatch run_compute_annot_sd.sh blood_GTEx_0.24 strength baseline
sbatch run_compute_annot_sd.sh blood_GTEx_0.24 page_rank baseline

sbatch run_compute_annot_sd.sh blood_GTEx_0.26 betweeness baseline
sbatch run_compute_annot_sd.sh blood_GTEx_0.26 closeness baseline
sbatch run_compute_annot_sd.sh blood_GTEx_0.26 degree baseline
sbatch run_compute_annot_sd.sh blood_GTEx_0.26 eigen_centrality baseline
sbatch run_compute_annot_sd.sh blood_GTEx_0.26 maximum_weight baseline
sbatch run_compute_annot_sd.sh blood_GTEx_0.26 strength baseline
sbatch run_compute_annot_sd.sh blood_GTEx_0.26 page_rank baseline

sbatch run_compute_annot_sd.sh blood_GTEx_0.28 betweeness baseline
sbatch run_compute_annot_sd.sh blood_GTEx_0.28 closeness baseline
sbatch run_compute_annot_sd.sh blood_GTEx_0.28 degree baseline
sbatch run_compute_annot_sd.sh blood_GTEx_0.28 eigen_centrality baseline
sbatch run_compute_annot_sd.sh blood_GTEx_0.28 maximum_weight baseline
sbatch run_compute_annot_sd.sh blood_GTEx_0.28 strength baseline
sbatch run_compute_annot_sd.sh blood_GTEx_0.28 page_rank baseline

sbatch run_compute_annot_sd.sh blood_GTEx_0.30 betweeness baseline
sbatch run_compute_annot_sd.sh blood_GTEx_0.30 closeness baseline
sbatch run_compute_annot_sd.sh blood_GTEx_0.30 degree baseline
sbatch run_compute_annot_sd.sh blood_GTEx_0.30 eigen_centrality baseline
sbatch run_compute_annot_sd.sh blood_GTEx_0.30 maximum_weight baseline
sbatch run_compute_annot_sd.sh blood_GTEx_0.30 strength baseline
sbatch run_compute_annot_sd.sh blood_GTEx_0.30 page_rank baseline

sbatch run_compute_annot_sd.sh blood_GTEx_0.32 betweeness baseline
sbatch run_compute_annot_sd.sh blood_GTEx_0.32 closeness baseline
sbatch run_compute_annot_sd.sh blood_GTEx_0.32 degree baseline
sbatch run_compute_annot_sd.sh blood_GTEx_0.32 eigen_centrality baseline
sbatch run_compute_annot_sd.sh blood_GTEx_0.32 maximum_weight baseline
sbatch run_compute_annot_sd.sh blood_GTEx_0.32 strength baseline
sbatch run_compute_annot_sd.sh blood_GTEx_0.32 page_rank baseline

sbatch run_compute_annot_sd.sh CNS_consensus_0.20 betweeness baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.20 closeness baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.20 degree baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.20 eigen_centrality baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.20 maximum_weight baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.20 strength baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.20 page_rank baseline

sbatch run_compute_annot_sd.sh CNS_consensus_0.22 betweeness baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.22 closeness baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.22 degree baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.22 eigen_centrality baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.22 maximum_weight baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.22 strength baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.22 page_rank baseline

sbatch run_compute_annot_sd.sh CNS_consensus_0.24 betweeness baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.24 closeness baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.24 degree baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.24 eigen_centrality baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.24 maximum_weight baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.24 strength baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.24 page_rank baseline

sbatch run_compute_annot_sd.sh CNS_consensus_0.26 betweeness baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.26 closeness baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.26 degree baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.26 eigen_centrality baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.26 maximum_weight baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.26 strength baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.26 page_rank baseline

sbatch run_compute_annot_sd.sh CNS_consensus_0.28 betweeness baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.28 closeness baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.28 degree baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.28 eigen_centrality baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.28 maximum_weight baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.28 strength baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.28 page_rank baseline

sbatch run_compute_annot_sd.sh CNS_consensus_0.30 betweeness baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.30 closeness baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.30 degree baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.30 eigen_centrality baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.30 maximum_weight baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.30 strength baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.30 page_rank baseline

sbatch run_compute_annot_sd.sh CNS_consensus_0.32 betweeness baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.32 closeness baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.32 degree baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.32 eigen_centrality baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.32 maximum_weight baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.32 strength baseline
sbatch run_compute_annot_sd.sh CNS_consensus_0.32 page_rank baseline

sbatch run_compute_annot_sd.sh CNS_GTEx_0.22 betweeness baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.22 closeness baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.22 degree baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.22 eigen_centrality baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.22 maximum_weight baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.22 strength baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.22 page_rank baseline

sbatch run_compute_annot_sd.sh CNS_GTEx_0.24 betweeness baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.24 closeness baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.24 degree baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.24 eigen_centrality baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.24 maximum_weight baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.24 strength baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.24 page_rank baseline

sbatch run_compute_annot_sd.sh CNS_GTEx_0.26 betweeness baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.26 closeness baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.26 degree baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.26 eigen_centrality baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.26 maximum_weight baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.26 strength baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.26 page_rank baseline

sbatch run_compute_annot_sd.sh CNS_GTEx_0.28 betweeness baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.28 closeness baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.28 degree baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.28 eigen_centrality baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.28 maximum_weight baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.28 strength baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.28 page_rank baseline

sbatch run_compute_annot_sd.sh CNS_GTEx_0.30 betweeness baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.30 closeness baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.30 degree baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.30 eigen_centrality baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.30 maximum_weight baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.30 strength baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.30 page_rank baseline

sbatch run_compute_annot_sd.sh CNS_GTEx_0.32 betweeness baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.32 closeness baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.32 degree baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.32 eigen_centrality baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.32 maximum_weight baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.32 strength baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.32 page_rank baseline

sbatch run_compute_annot_sd.sh CNS_GTEx_0.34 betweeness baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.34 closeness baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.34 degree baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.34 eigen_centrality baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.34 maximum_weight baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.34 strength baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.34 page_rank baseline

sbatch run_compute_annot_sd.sh CNS_GTEx_0.36 betweeness baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.36 closeness baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.36 degree baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.36 eigen_centrality baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.36 maximum_weight baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.36 strength baseline
sbatch run_compute_annot_sd.sh CNS_GTEx_0.36 page_rank baseline