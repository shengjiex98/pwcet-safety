#!/bin/bash

module use $HOME/modulefiles
module add julia/1.10.0

# for th in 2 4 8 16
for th in 16
do
    # for b in 100000 1000000
    # for b in 10000000
    for b in 10000
    do
        for n in 100
        do
            for q in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99 0.999
            do
                for h in 0.005 0.01 0.015 0.02 0.025 0.03 0.035 0.04 0.045 0.05 0.1 0.2 0.5 1.0
                do
                    # 20GB of memory, 16 core, and a one day time limit.
                    (set -x; sbatch -p general -c $th -t 01-00:00:00 --mem=8g --wrap="julia --project --threads=$th generate_data.jl normal $b $q $h")
                    # (set -x; sbatch -p general -c $th -t 01-00:00:00 --mem=20g --wrap="julia --project --threads=$th generate_data.jl batch $b $q $h $n")
                done
            done
        done
    done
done
