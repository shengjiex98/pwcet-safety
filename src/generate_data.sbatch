#!/bin/bash

# 4GB of memory, 16 core, and a one day time limit.
#SBATCH --partition=general
#SBATCH --ntasks=16
#SBATCH --time=1-00:00:00
#SBATCH --mem=16g

module use $HOME/modulefiles
module add julia/1.10.0

# for b in 100000 1000000
for b in 10000
do
    for q in 0.9 0.99 0.999
    do
        for n in 100
        do
            echo $b $q $n
            julia --project --threads=16 generate_data.jl batch $b $q $n
        done
    done
done
