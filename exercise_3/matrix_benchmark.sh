#!/bin/sh
#SBATCH --job-name="matrix_benchmark"
#SBATCH --partition=gpu-a100-small
#SBATCH --account=Education-EEMCS-Courses-WI4049TU
#SBATCH --time=00:05:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-task=1
#SBATCH --mem-per-cpu=4GB
module load 2025
module load openmpi
module load slurm
module load cuda

for m_size in 50 500 2000 4000
do
    for block_size in 32 64 100
    do
        srun build/power_gpu_${block_size}.x -size $m_size > ./step_2/output_s${m_size}_b${block_size}.out
    done
done
