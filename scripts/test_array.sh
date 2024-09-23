# (set -x; SLURM_ARRAY_JOB_ID=1 SLURM_ARRAY_TASK_ID=30 COMPARE_MODE=y julia --project --threads=1 generate_data_mpc.jl)
# (set -x; SLURM_ARRAY_JOB_ID=1 SLURM_ARRAY_TASK_ID=30 COMPARE_MODE=y julia --project --threads=1 generate_grid_mpc.jl)
(set -x; SLURM_ARRAY_JOB_ID=1 SLURM_ARRAY_TASK_ID=30 BATCHSIZE=100 julia --project --threads=4 generate_grid.jl)
