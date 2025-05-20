#!/bin/bash
#SBATCH -N 1
#SBATCH -t 10:00:00
#SBATCH --mail-type ALL
#SBATCH --mail-user andrescb@kth.se
#SBATCH --output ./check.out
#SBATCH --error  ./check.error
#SBATCH --account=Berzelius-2024-98

cd ../
singularity run --nv ../../acbtensorflow_2.6.1-gpu.sif python check_data.py
