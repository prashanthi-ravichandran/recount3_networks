#!/bin/bash -l
#SBATCH
#SBATCH --job-name=submit_ll_context_specific
#SBATCH --mail-type=END
#SBATCH --mail-user=pravich2@jhu.edu
#SBATCH --time=0:20:0
#SBATCH --partition=defq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

sbatch run_held_out_ll_context_specific.sh blood 1
sbatch run_held_out_ll_context_specific.sh blood 5
sbatch run_held_out_ll_context_specific.sh blood 10
sbatch run_held_out_ll_context_specific.sh blood 15
sbatch run_held_out_ll_context_specific.sh blood 20
sbatch run_held_out_ll_context_specific.sh blood 25
sbatch run_held_out_ll_context_specific.sh blood 30
sbatch run_held_out_ll_context_specific.sh blood 35
sbatch run_held_out_ll_context_specific.sh blood 40
sbatch run_held_out_ll_context_specific.sh blood 45
sbatch run_held_out_ll_context_specific.sh blood 50
sbatch run_held_out_ll_context_specific.sh blood 55
sbatch run_held_out_ll_context_specific.sh blood 60
sbatch run_held_out_ll_context_specific.sh blood 64

sbatch run_held_out_ll_context_specific.sh central_nervous_system 1
sbatch run_held_out_ll_context_specific.sh central_nervous_system 5
sbatch run_held_out_ll_context_specific.sh central_nervous_system 10
sbatch run_held_out_ll_context_specific.sh central_nervous_system 15
sbatch run_held_out_ll_context_specific.sh central_nervous_system 20
sbatch run_held_out_ll_context_specific.sh central_nervous_system 25
sbatch run_held_out_ll_context_specific.sh central_nervous_system 30
sbatch run_held_out_ll_context_specific.sh central_nervous_system 35
sbatch run_held_out_ll_context_specific.sh central_nervous_system 40


