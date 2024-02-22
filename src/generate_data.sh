#!/bin/bash

module use $HOME/modulefiles
module add julia/1.10.0

for th in 2 4 8 16
do
    # for b in 100000 1000000
    for b in 10000000
    do
        for q in 0.9 0.99 0.999
        do
            # 20GB of memory, 16 core, and a one day time limit.
            (set -x; sbatch -p general -c $th -t 01-00:00:00 --mem=20g --wrap="julia --project --threads=$th generate_data.jl $b $q")
        done
    done
done
