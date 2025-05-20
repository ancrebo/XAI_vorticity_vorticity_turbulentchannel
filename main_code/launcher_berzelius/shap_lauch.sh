#!/bin/bash
#SBATCH --gpus 1
#SBATCH -N 1
#SBATCH -t 60:00:00
#SBATCH --mail-type ALL
#SBATCH --mail-user andrescb@kth.se
#SBATCH --output ./shap.out
#SBATCH --error  ./shap.error
#SBATCH --account=Berzelius-2024-98

cd ../
singularity run --nv ../../acbtensorflow_2.6.1-gpu.sif python main_SHAP.py
