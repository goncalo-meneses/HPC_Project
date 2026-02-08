#!/bin/bash
#SBATCH --job-name=mm_scaling
#SBATCH --partition=compute-p1
#SBATCH --time=00:10:00
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --account=Education-EEMCS-Courses-WI4049TU

#SBATCH --output=mm_scaling_%j.out

module load 2025 openmpi slurm

for N in 200 400 800 1000; do
  echo "=== Serial N=$N ==="
  srun -n 1 serial_MM_product $N
  echo
  for p in 1 2 8 24 48 64; do
    if [ $p -le $SLURM_NTASKS ]; then
      echo "=== N=$N, p=$p ==="
      srun -n $p MM_product $N
      echo
    fi
  done
done
