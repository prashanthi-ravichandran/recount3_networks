#!/bin/bash -l
#SBATCH
#SBATCH --job-name=submit_make_annot
#SBATCH --mail-type=END
#SBATCH --mail-user=pravich2@jhu.edu
#SBATCH --time=0:20:0
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

for trait in $(cat /data/abattle4/prashanthi/recount3_paper/data/s_LDSC/UKBB/bgz_traits.txt);
do
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh all_consensus_0.14 betweeness $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh all_consensus_0.14 closeness $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh all_consensus_0.14 degree $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh all_consensus_0.14 eigen_centrality $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh all_consensus_0.14 maximum_weight $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh all_consensus_0.14 strength $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh all_consensus_0.14 page_rank $trait 
  
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh all_consensus_0.16 betweeness $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh all_consensus_0.16 closeness $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh all_consensus_0.16 degree $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh all_consensus_0.16 eigen_centrality $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh all_consensus_0.16 maximum_weight $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh all_consensus_0.16 strength $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh all_consensus_0.16 page_rank $trait 
  
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh all_consensus_0.18 betweeness $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh all_consensus_0.18 closeness $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh all_consensus_0.18 degree $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh all_consensus_0.18 eigen_centrality $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh all_consensus_0.18 maximum_weight $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh all_consensus_0.18 strength $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh all_consensus_0.18 page_rank $trait 
  
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh all_consensus_0.20 betweeness $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh all_consensus_0.20 closeness $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh all_consensus_0.20 degree $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh all_consensus_0.20 eigen_centrality $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh all_consensus_0.20 maximum_weight $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh all_consensus_0.20 strength $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh all_consensus_0.20 page_rank $trait 
  
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh normal_consensus_0.14 betweeness $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh normal_consensus_0.14 closeness $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh normal_consensus_0.14 degree $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh normal_consensus_0.14 eigen_centrality $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh normal_consensus_0.14 maximum_weight $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh normal_consensus_0.14 strength $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh normal_consensus_0.14 page_rank $trait 
  
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh normal_consensus_0.16 betweeness $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh normal_consensus_0.16 closeness $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh normal_consensus_0.16 degree $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh normal_consensus_0.16 eigen_centrality $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh normal_consensus_0.16 maximum_weight $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh normal_consensus_0.16 strength $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh normal_consensus_0.16 page_rank $trait 
  
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh normal_consensus_0.18 betweeness $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh normal_consensus_0.18 closeness $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh normal_consensus_0.18 degree $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh normal_consensus_0.18 eigen_centrality $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh normal_consensus_0.18 maximum_weight $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh normal_consensus_0.18 strength $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh normal_consensus_0.18 page_rank $trait 
  
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh normal_consensus_0.20 betweeness $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh normal_consensus_0.20 closeness $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh normal_consensus_0.20 degree $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh normal_consensus_0.20 eigen_centrality $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh normal_consensus_0.20 maximum_weight $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh normal_consensus_0.20 strength $trait 
  sbatch run_sLDSC_UKBB_bgz_all_genes.sh normal_consensus_0.20 page_rank $trait 
done

for trait in $(cat /data/abattle4/prashanthi/recount3_paper/data/s_LDSC/UKBB/bgz_traits.txt);
do
  sbatch run_sLDSC_UKBB_bgz_baseline.sh all_consensus_0.14 betweeness $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh all_consensus_0.14 closeness $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh all_consensus_0.14 degree $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh all_consensus_0.14 eigen_centrality $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh all_consensus_0.14 maximum_weight $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh all_consensus_0.14 strength $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh all_consensus_0.14 page_rank $trait 

  sbatch run_sLDSC_UKBB_bgz_baseline.sh all_consensus_0.16 betweeness $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh all_consensus_0.16 closeness $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh all_consensus_0.16 degree $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh all_consensus_0.16 eigen_centrality $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh all_consensus_0.16 maximum_weight $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh all_consensus_0.16 strength $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh all_consensus_0.16 page_rank $trait 

  sbatch run_sLDSC_UKBB_bgz_baseline.sh all_consensus_0.18 betweeness $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh all_consensus_0.18 closeness $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh all_consensus_0.18 degree $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh all_consensus_0.18 eigen_centrality $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh all_consensus_0.18 maximum_weight $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh all_consensus_0.18 strength $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh all_consensus_0.18 page_rank $trait 

  sbatch run_sLDSC_UKBB_bgz_baseline.sh all_consensus_0.20 betweeness $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh all_consensus_0.20 closeness $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh all_consensus_0.20 degree $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh all_consensus_0.20 eigen_centrality $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh all_consensus_0.20 maximum_weight $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh all_consensus_0.20 strength $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh all_consensus_0.20 page_rank $trait 

  sbatch run_sLDSC_UKBB_bgz_baseline.sh normal_consensus_0.14 betweeness $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh normal_consensus_0.14 closeness $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh normal_consensus_0.14 degree $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh normal_consensus_0.14 eigen_centrality $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh normal_consensus_0.14 maximum_weight $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh normal_consensus_0.14 strength $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh normal_consensus_0.14 page_rank $trait 

  sbatch run_sLDSC_UKBB_bgz_baseline.sh normal_consensus_0.16 betweeness $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh normal_consensus_0.16 closeness $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh normal_consensus_0.16 degree $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh normal_consensus_0.16 eigen_centrality $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh normal_consensus_0.16 maximum_weight $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh normal_consensus_0.16 strength $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh normal_consensus_0.16 page_rank $trait 

  sbatch run_sLDSC_UKBB_bgz_baseline.sh normal_consensus_0.18 betweeness $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh normal_consensus_0.18 closeness $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh normal_consensus_0.18 degree $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh normal_consensus_0.18 eigen_centrality $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh normal_consensus_0.18 maximum_weight $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh normal_consensus_0.18 strength $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh normal_consensus_0.18 page_rank $trait 

  sbatch run_sLDSC_UKBB_bgz_baseline.sh normal_consensus_0.20 betweeness $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh normal_consensus_0.20 closeness $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh normal_consensus_0.20 degree $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh normal_consensus_0.20 eigen_centrality $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh normal_consensus_0.20 maximum_weight $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh normal_consensus_0.20 strength $trait 
  sbatch run_sLDSC_UKBB_bgz_baseline.sh normal_consensus_0.20 page_rank $trait 
done

for trait in $(cat /data/abattle4/prashanthi/recount3_paper/data/s_LDSC/UKBB/gz_traits.txt);
do
  sbatch run_sLDSC_UKBB_gz_all_genes.sh all_consensus_0.14 betweeness $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh all_consensus_0.14 closeness $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh all_consensus_0.14 degree $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh all_consensus_0.14 eigen_centrality $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh all_consensus_0.14 maximum_weight $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh all_consensus_0.14 strength $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh all_consensus_0.14 page_rank $trait 

  sbatch run_sLDSC_UKBB_gz_all_genes.sh all_consensus_0.16 betweeness $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh all_consensus_0.16 closeness $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh all_consensus_0.16 degree $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh all_consensus_0.16 eigen_centrality $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh all_consensus_0.16 maximum_weight $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh all_consensus_0.16 strength $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh all_consensus_0.16 page_rank $trait 

  sbatch run_sLDSC_UKBB_gz_all_genes.sh all_consensus_0.18 betweeness $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh all_consensus_0.18 closeness $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh all_consensus_0.18 degree $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh all_consensus_0.18 eigen_centrality $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh all_consensus_0.18 maximum_weight $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh all_consensus_0.18 strength $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh all_consensus_0.18 page_rank $trait 

  sbatch run_sLDSC_UKBB_gz_all_genes.sh all_consensus_0.20 betweeness $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh all_consensus_0.20 closeness $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh all_consensus_0.20 degree $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh all_consensus_0.20 eigen_centrality $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh all_consensus_0.20 maximum_weight $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh all_consensus_0.20 strength $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh all_consensus_0.20 page_rank $trait 

  sbatch run_sLDSC_UKBB_gz_all_genes.sh normal_consensus_0.14 betweeness $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh normal_consensus_0.14 closeness $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh normal_consensus_0.14 degree $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh normal_consensus_0.14 eigen_centrality $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh normal_consensus_0.14 maximum_weight $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh normal_consensus_0.14 strength $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh normal_consensus_0.14 page_rank $trait 

  sbatch run_sLDSC_UKBB_gz_all_genes.sh normal_consensus_0.16 betweeness $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh normal_consensus_0.16 closeness $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh normal_consensus_0.16 degree $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh normal_consensus_0.16 eigen_centrality $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh normal_consensus_0.16 maximum_weight $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh normal_consensus_0.16 strength $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh normal_consensus_0.16 page_rank $trait 

  sbatch run_sLDSC_UKBB_gz_all_genes.sh normal_consensus_0.18 betweeness $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh normal_consensus_0.18 closeness $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh normal_consensus_0.18 degree $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh normal_consensus_0.18 eigen_centrality $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh normal_consensus_0.18 maximum_weight $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh normal_consensus_0.18 strength $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh normal_consensus_0.18 page_rank $trait 

  sbatch run_sLDSC_UKBB_gz_all_genes.sh normal_consensus_0.20 betweeness $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh normal_consensus_0.20 closeness $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh normal_consensus_0.20 degree $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh normal_consensus_0.20 eigen_centrality $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh normal_consensus_0.20 maximum_weight $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh normal_consensus_0.20 strength $trait 
  sbatch run_sLDSC_UKBB_gz_all_genes.sh normal_consensus_0.20 page_rank $trait 
done

for trait in $(cat /data/abattle4/prashanthi/recount3_paper/data/s_LDSC/UKBB/gz_traits.txt);
do
sbatch run_sLDSC_UKBB_gz_baseline.sh all_consensus_0.14 betweeness $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh all_consensus_0.14 closeness $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh all_consensus_0.14 degree $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh all_consensus_0.14 eigen_centrality $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh all_consensus_0.14 maximum_weight $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh all_consensus_0.14 strength $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh all_consensus_0.14 page_rank $trait 

sbatch run_sLDSC_UKBB_gz_baseline.sh all_consensus_0.16 betweeness $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh all_consensus_0.16 closeness $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh all_consensus_0.16 degree $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh all_consensus_0.16 eigen_centrality $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh all_consensus_0.16 maximum_weight $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh all_consensus_0.16 strength $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh all_consensus_0.16 page_rank $trait 

sbatch run_sLDSC_UKBB_gz_baseline.sh all_consensus_0.18 betweeness $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh all_consensus_0.18 closeness $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh all_consensus_0.18 degree $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh all_consensus_0.18 eigen_centrality $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh all_consensus_0.18 maximum_weight $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh all_consensus_0.18 strength $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh all_consensus_0.18 page_rank $trait 

sbatch run_sLDSC_UKBB_gz_baseline.sh all_consensus_0.20 betweeness $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh all_consensus_0.20 closeness $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh all_consensus_0.20 degree $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh all_consensus_0.20 eigen_centrality $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh all_consensus_0.20 maximum_weight $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh all_consensus_0.20 strength $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh all_consensus_0.20 page_rank $trait 

sbatch run_sLDSC_UKBB_gz_baseline.sh normal_consensus_0.14 betweeness $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh normal_consensus_0.14 closeness $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh normal_consensus_0.14 degree $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh normal_consensus_0.14 eigen_centrality $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh normal_consensus_0.14 maximum_weight $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh normal_consensus_0.14 strength $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh normal_consensus_0.14 page_rank $trait 

sbatch run_sLDSC_UKBB_gz_baseline.sh normal_consensus_0.16 betweeness $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh normal_consensus_0.16 closeness $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh normal_consensus_0.16 degree $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh normal_consensus_0.16 eigen_centrality $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh normal_consensus_0.16 maximum_weight $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh normal_consensus_0.16 strength $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh normal_consensus_0.16 page_rank $trait 

sbatch run_sLDSC_UKBB_gz_baseline.sh normal_consensus_0.18 betweeness $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh normal_consensus_0.18 closeness $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh normal_consensus_0.18 degree $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh normal_consensus_0.18 eigen_centrality $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh normal_consensus_0.18 maximum_weight $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh normal_consensus_0.18 strength $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh normal_consensus_0.18 page_rank $trait 

sbatch run_sLDSC_UKBB_gz_baseline.sh normal_consensus_0.20 betweeness $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh normal_consensus_0.20 closeness $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh normal_consensus_0.20 degree $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh normal_consensus_0.20 eigen_centrality $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh normal_consensus_0.20 maximum_weight $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh normal_consensus_0.20 strength $trait 
sbatch run_sLDSC_UKBB_gz_baseline.sh normal_consensus_0.20 page_rank $trait 
done
