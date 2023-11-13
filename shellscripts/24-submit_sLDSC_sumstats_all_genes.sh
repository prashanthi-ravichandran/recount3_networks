#!/bin/bash -l
#SBATCH
#SBATCH --job-name=submit_make_annot
#SBATCH --mail-type=END
#SBATCH --mail-user=pravich2@jhu.edu
#SBATCH --time=0:20:0
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

for trait in $(cat /data/abattle4/prashanthi/recount3_paper/data/s_LDSC/all_sumstats/traits.txt);
do
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.14 betweeness $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.14 closeness $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.14 degree $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.14 eigen_centrality $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.14 maximum_weight $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.14 strength $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.14 page_rank $trait all_sumstats
  
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.16 betweeness $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.16 closeness $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.16 degree $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.16 eigen_centrality $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.16 maximum_weight $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.16 strength $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.16 page_rank $trait all_sumstats
  
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.18 betweeness $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.18 closeness $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.18 degree $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.18 eigen_centrality $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.18 maximum_weight $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.18 strength $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.18 page_rank $trait all_sumstats
  
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.20 betweeness $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.20 closeness $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.20 degree $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.20 eigen_centrality $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.20 maximum_weight $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.20 strength $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.20 page_rank $trait all_sumstats
  
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.14 betweeness $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.14 closeness $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.14 degree $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.14 eigen_centrality $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.14 maximum_weight $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.14 strength $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.14 page_rank $trait all_sumstats
  
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.16 betweeness $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.16 closeness $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.16 degree $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.16 eigen_centrality $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.16 maximum_weight $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.16 strength $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.16 page_rank $trait all_sumstats
  
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.18 betweeness $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.18 closeness $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.18 degree $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.18 eigen_centrality $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.18 maximum_weight $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.18 strength $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.18 page_rank $trait all_sumstats
  
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.20 betweeness $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.20 closeness $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.20 degree $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.20 eigen_centrality $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.20 maximum_weight $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.20 strength $trait all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.20 page_rank $trait all_sumstats
done

for trait in $(cat /data/abattle4/prashanthi/recount3_paper/data/s_LDSC/independent_sumstats/traits.txt);
do
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.14 betweeness $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.14 closeness $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.14 degree $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.14 eigen_centrality $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.14 maximum_weight $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.14 strength $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.14 page_rank $trait independent_sumstats

  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.16 betweeness $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.16 closeness $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.16 degree $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.16 eigen_centrality $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.16 maximum_weight $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.16 strength $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.16 page_rank $trait independent_sumstats

  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.18 betweeness $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.18 closeness $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.18 degree $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.18 eigen_centrality $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.18 maximum_weight $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.18 strength $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.18 page_rank $trait independent_sumstats

  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.20 betweeness $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.20 closeness $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.20 degree $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.20 eigen_centrality $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.20 maximum_weight $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.20 strength $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh all_consensus_0.20 page_rank $trait independent_sumstats

  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.14 betweeness $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.14 closeness $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.14 degree $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.14 eigen_centrality $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.14 maximum_weight $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.14 strength $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.14 page_rank $trait independent_sumstats

  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.16 betweeness $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.16 closeness $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.16 degree $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.16 eigen_centrality $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.16 maximum_weight $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.16 strength $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.16 page_rank $trait independent_sumstats

  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.18 betweeness $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.18 closeness $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.18 degree $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.18 eigen_centrality $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.18 maximum_weight $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.18 strength $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.18 page_rank $trait independent_sumstats

  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.20 betweeness $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.20 closeness $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.20 degree $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.20 eigen_centrality $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.20 maximum_weight $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.20 strength $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh normal_consensus_0.20 page_rank $trait independent_sumstats
done

for trait in $(cat /data/abattle4/prashanthi/recount3_paper/data/s_LDSC/independent_sumstats/blood_traits.txt);
do
  sbatch run_sLDSC_sumstats_all_genes.sh blood_consensus_0.18 betweeness $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_consensus_0.18 closeness $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_consensus_0.18 degree $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_consensus_0.18 eigen_centrality $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_consensus_0.18 maximum_weight $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_consensus_0.18 strength $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_consensus_0.18 page_rank $trait independent_sumstats

  sbatch run_sLDSC_sumstats_all_genes.sh blood_consensus_0.20 betweeness $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_consensus_0.20 closeness $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_consensus_0.20 degree $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_consensus_0.20 eigen_centrality $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_consensus_0.20 maximum_weight $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_consensus_0.20 strength $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_consensus_0.20 page_rank $trait independent_sumstats

  sbatch run_sLDSC_sumstats_all_genes.sh blood_consensus_0.22 betweeness $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_consensus_0.22 closeness $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_consensus_0.22 degree $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_consensus_0.22 eigen_centrality $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_consensus_0.22 maximum_weight $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_consensus_0.22 strength $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_consensus_0.22 page_rank $trait independent_sumstats

  sbatch run_sLDSC_sumstats_all_genes.sh blood_consensus_0.24 betweeness $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_consensus_0.24 closeness $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_consensus_0.24 degree $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_consensus_0.24 eigen_centrality $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_consensus_0.24 maximum_weight $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_consensus_0.24 strength $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_consensus_0.24 page_rank $trait independent_sumstats

  sbatch run_sLDSC_sumstats_all_genes.sh blood_consensus_0.26 betweeness $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_consensus_0.26 closeness $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_consensus_0.26 degree $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_consensus_0.26 eigen_centrality $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_consensus_0.26 maximum_weight $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_consensus_0.26 strength $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_consensus_0.26 page_rank $trait independent_sumstats

  sbatch run_sLDSC_sumstats_all_genes.sh blood_GTEx_0.24 betweeness $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_GTEx_0.24 closeness $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_GTEx_0.24 degree $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_GTEx_0.24 eigen_centrality $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_GTEx_0.24 maximum_weight $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_GTEx_0.24 strength $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_GTEx_0.24 page_rank $trait independent_sumstats

  sbatch run_sLDSC_sumstats_all_genes.sh blood_GTEx_0.26 betweeness $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_GTEx_0.26 closeness $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_GTEx_0.26 degree $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_GTEx_0.26 eigen_centrality $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_GTEx_0.26 maximum_weight $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_GTEx_0.26 strength $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_GTEx_0.26 page_rank $trait independent_sumstats

  sbatch run_sLDSC_sumstats_all_genes.sh blood_GTEx_0.28 betweeness $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_GTEx_0.28 closeness $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_GTEx_0.28 degree $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_GTEx_0.28 eigen_centrality $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_GTEx_0.28 maximum_weight $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_GTEx_0.28 strength $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_GTEx_0.28 page_rank $trait independent_sumstats

  sbatch run_sLDSC_sumstats_all_genes.sh blood_GTEx_0.30 betweeness $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_GTEx_0.30 closeness $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_GTEx_0.30 degree $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_GTEx_0.30 eigen_centrality $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_GTEx_0.30 maximum_weight $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_GTEx_0.30 strength $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_GTEx_0.30 page_rank $trait independent_sumstats

  sbatch run_sLDSC_sumstats_all_genes.sh blood_GTEx_0.32 betweeness $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_GTEx_0.32 closeness $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_GTEx_0.32 degree $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_GTEx_0.32 eigen_centrality $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_GTEx_0.32 maximum_weight $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_GTEx_0.32 strength $trait independent_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh blood_GTEx_0.32 page_rank $trait independent_sumstats
done

for trait in $(cat /data/abattle4/prashanthi/recount3_paper/data/s_LDSC/all_sumstats/CNS_traits.txt);
do
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.20 betweeness $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.20 closeness $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.20 degree $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.20 eigen_centrality $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.20 maximum_weight $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.20 strength $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.20 page_rank $chr all_sumstats

  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.22 betweeness $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.22 closeness $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.22 degree $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.22 eigen_centrality $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.22 maximum_weight $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.22 strength $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.22 page_rank $chr all_sumstats

  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.24 betweeness $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.24 closeness $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.24 degree $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.24 eigen_centrality $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.24 maximum_weight $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.24 strength $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.24 page_rank $chr all_sumstats

  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.26 betweeness $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.26 closeness $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.26 degree $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.26 eigen_centrality $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.26 maximum_weight $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.26 strength $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.26 page_rank $chr all_sumstats

  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.28 betweeness $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.28 closeness $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.28 degree $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.28 eigen_centrality $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.28 maximum_weight $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.28 strength $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.28 page_rank $chr all_sumstats

  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.30 betweeness $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.30 closeness $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.30 degree $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.30 eigen_centrality $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.30 maximum_weight $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.30 strength $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.30 page_rank $chr all_sumstats

  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.32 betweeness $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.32 closeness $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.32 degree $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.32 eigen_centrality $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.32 maximum_weight $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.32 strength $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_consensus_0.32 page_rank $chr all_sumstats

  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.22 betweeness $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.22 closeness $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.22 degree $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.22 eigen_centrality $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.22 maximum_weight $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.22 strength $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.22 page_rank $chr all_sumstats

  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.24 betweeness $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.24 closeness $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.24 degree $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.24 eigen_centrality $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.24 maximum_weight $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.24 strength $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.24 page_rank $chr all_sumstats

  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.26 betweeness $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.26 closeness $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.26 degree $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.26 eigen_centrality $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.26 maximum_weight $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.26 strength $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.26 page_rank $chr all_sumstats

  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.28 betweeness $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.28 closeness $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.28 degree $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.28 eigen_centrality $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.28 maximum_weight $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.28 strength $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.28 page_rank $chr all_sumstats

  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.30 betweeness $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.30 closeness $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.30 degree $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.30 eigen_centrality $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.30 maximum_weight $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.30 strength $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.30 page_rank $chr all_sumstats

  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.32 betweeness $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.32 closeness $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.32 degree $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.32 eigen_centrality $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.32 maximum_weight $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.32 strength $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.32 page_rank $chr all_sumstats

  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.34 betweeness $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.34 closeness $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.34 degree $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.34 eigen_centrality $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.34 maximum_weight $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.34 strength $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.34 page_rank $chr all_sumstats

  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.36 betweeness $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.36 closeness $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.36 degree $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.36 eigen_centrality $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.36 maximum_weight $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.36 strength $chr all_sumstats
  sbatch run_sLDSC_sumstats_all_genes.sh CNS_GTEx_0.36 page_rank $chr all_sumstats
done


