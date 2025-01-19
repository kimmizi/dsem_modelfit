#!/bin/bash
#
#SBATCH --array=3-256
#SBATCH --cpus-per-task=4
# Check if simulation number was provided
if [ $# -ne 1 ]; then
    echo "Usage: sbatch run_simulation.sh <simulation_number>"
    exit 1
fi

# Get simulation number from command line
SIM_NUM=$1

# Set job name (optional, as SLURM's #SBATCH job-name is static)
JOB_NAME="dsem_sim_${SIM_NUM}"
# Dynamically set log files using the SLURM_JOB_ID
OUTPUT_FILE="log/sim${SIM_NUM}_${SLURM_JOB_ID}_%a.out"
ERROR_FILE="log/sim${SIM_NUM}_${SLURM_JOB_ID}_%a.err"
# Set job name and log files with simulation number
#SBATCH --job-name=dsem_sim_${SIM_NUM}
#SBATCH --time=48:00:00
#SBATCH --partition=cpu-single
#SBATCH --mem-per-cpu=64000
module purge
module load math/R

Rscript run_sim_new.R