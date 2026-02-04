#!/bin/sh
#SBATCH --job-name="shared_v_global"
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

for m_size in 500 2000 4000 8000 10000
do
    srun build/power_gpu_32.x -size $m_size > ./step_1/shared_s${m_size}.out
    srun build/global_power_gpu.x -size $m_size > ./step_1/global_s${m_size}.out
done
