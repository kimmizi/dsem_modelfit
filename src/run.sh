#!/bin/bash
#
#SBATCH --array=400-410
#SBATCH --cpus-per-task=4
#SBATCH --output=log/"sim_%j_%a.out"
#SBATCH --error=log/"sim_%j_%a.err"
#SBATCH --job-name=dsem_sim_big
#SBATCH --time=120:00:00
#SBATCH --partition=cpu-single
#SBATCH --mem-per-cpu=256000

# Get simulation number from command line
SIM_NUM=$1

module purge
module load math/R

Rscript run_par_new.R ${SIM_NUM}
