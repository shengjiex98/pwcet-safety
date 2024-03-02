#!/bin/bash

module use $HOME/modulefiles
module add julia/1.10.0

# for th in 2 4 8 16
for th in 16
do
    # for b in 100000 1000000
    # for b in 10000000
    for b in 10
    do
        for n in 10
        do
            for q in 0.1 0.2
            do
            for h in 0.005 0.01
                # 20GB of memory, 16 core, and a one day time limit.
                (set -x; sbatch -p general -c $th -t 01-00:00:00 --mem=20g --wrap="julia --project --threads=$th generate_data.jl normal $b $q")
                # (set -x; sbatch -p general -c $th -t 01-00:00:00 --mem=20g --wrap="julia --project --threads=$th generate_data.jl batch $b $q $n")
            done
        done
    done
done
