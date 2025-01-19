#!/bin/bash
#
#SBATCH --array=3-256
#SBATCH --cpus-per-task=4
# Check if simulation number was provided
# Ensure the log directory exists
mkdir -p log

# Check if simulation number was provided
if [ $# -ne 1 ]; then
    echo "Usage: sbatch run_simulation.sh <simulation_number>"
    exit 1
fi

# Get simulation number from command line
SIM_NUM=$1

# Set job name (optional, as SLURM's #SBATCH job-name is static)
JOB_NAME="dsem_sim_1"
#SBATCH --output=log/sim_w1_%a_%j.out  # Logs per array task (%a) and job ID (%j)
#SBATCH --error=log/sim__w1_%a_%j.err   # Logs errors to the same folder
# Set job name and log files with simulation number
#SBATCH --job-name=dsem_sim_${SIM_NUM}
#SBATCH --time=48:00:00
#SBATCH --partition=cpu-single
#SBATCH --mem-per-cpu=64000
module purge
module load math/R

Rscript run_par_new.R