#!/bin/bash

#SBATCH --job-name=pwcet-safety-array
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8g
#SBATCH --time=0-2:00:00
#SBATCH --mail-type=all                         # send email on job start, end and fault
#SBATCH --mail-user=sxunique@cs.unc.edu

echo "SLURM_ARRAY_JOB_ID: $SLURM_ARRAY_JOB_ID."
echo "Executing on the machine:" $(hostname)

module purge
module use $HOME/modulefiles
module add julia/1.10.0

echo "Start running Julia: " $(date)
SECONDS=0

(set -x; julia --project summarize_data.jl)

echo "End running Julia: " $(date)
echo "Elapsed time: $SECONDS"
