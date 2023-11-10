#!/bin/bash -l
#SBATCH
#SBATCH --job-name=submit_SNAP_enrichment
#SBATCH --mail-type=END
#SBATCH --mail-user=pravich2@jhu.edu
#SBATCH --time=0:20:0
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

sbatch run_SNAP_enrichment.sh adipose all 11
sbatch run_SNAP_enrichment.sh airway all 8
sbatch run_SNAP_enrichment.sh blood all 65
sbatch run_SNAP_enrichment.sh B_cells all 17
sbatch run_SNAP_enrichment.sh breast all 12
sbatch run_SNAP_enrichment.sh cardiac all 13
sbatch run_SNAP_enrichment.sh central_nervous_system all 53
sbatch run_SNAP_enrichment.sh colon all 10
sbatch run_SNAP_enrichment.sh eye all 2
sbatch run_SNAP_enrichment.sh esophagus all 3
sbatch run_SNAP_enrichment.sh fibroblasts all 13
sbatch run_SNAP_enrichment.sh hescs all 21
sbatch run_SNAP_enrichment.sh intestine all 9
sbatch run_SNAP_enrichment.sh ipscs all 14
sbatch run_SNAP_enrichment.sh kidney all 26
sbatch run_SNAP_enrichment.sh liver all 28
sbatch run_SNAP_enrichment.sh lung all 10
sbatch run_SNAP_enrichment.sh multipotent_cells all 6
sbatch run_SNAP_enrichment.sh myeloid_cells all 29
sbatch run_SNAP_enrichment.sh pancreas all 3
sbatch run_SNAP_enrichment.sh pbmcs_t_cells all 58
sbatch run_SNAP_enrichment.sh prostate all 6
sbatch run_SNAP_enrichment.sh skeletal_muscle all 8
sbatch run_SNAP_enrichment.sh skin all 20
sbatch run_SNAP_enrichment.sh stomach all 4
sbatch run_SNAP_enrichment.sh vascular all 4

sbatch run_SNAP_enrichment.sh adipose GTEx 2
sbatch run_SNAP_enrichment.sh blood GTEx 1
sbatch run_SNAP_enrichment.sh B_cells GTEx 1
sbatch run_SNAP_enrichment.sh breast GTEx 1
sbatch run_SNAP_enrichment.sh cardiac GTEx 2
sbatch run_SNAP_enrichment.sh central_nervous_system GTEx 13
sbatch run_SNAP_enrichment.sh colon GTEx 2
sbatch run_SNAP_enrichment.sh esophagus GTEx 3
sbatch run_SNAP_enrichment.sh fibroblasts GTEx 1
sbatch run_SNAP_enrichment.sh intestine GTEx 1
sbatch run_SNAP_enrichment.sh kidney GTEx 1
sbatch run_SNAP_enrichment.sh liver GTEx 1
sbatch run_SNAP_enrichment.sh lung GTEx 1
sbatch run_SNAP_enrichment.sh pancreas GTEx 1
sbatch run_SNAP_enrichment.sh prostate GTEx 1
sbatch run_SNAP_enrichment.sh skeletal_muscle GTEx 1
sbatch run_SNAP_enrichment.sh skin GTEx 2
sbatch run_SNAP_enrichment.sh stomach GTEx 1
sbatch run_SNAP_enrichment.sh vascular GTEx 3


