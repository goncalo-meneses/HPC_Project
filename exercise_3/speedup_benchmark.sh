#!/bin/sh
#SBATCH --job-name="speedup_benchmark"
#SBATCH --partition=gpu-a100-small
#SBATCH --account=Education-EEMCS-Courses-WI4049TU
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gpus-per-task=1
#SBATCH --mem-per-cpu=4GB
module load 2025
module load openmpi
module load slurm
module load cuda

mkdir -p step_3

for m_size in 500 1000 2000 4000 6000 8000 10000
do
    srun build/power_gpu_32.x -size $m_size > ./speedup_results/output_s${m_size}.out
done
