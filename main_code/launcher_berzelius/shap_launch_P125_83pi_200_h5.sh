#!/bin/bash
#SBATCH --gpus 1 -C "fat"
#SBATCH -J shap_0
#SBATCH -N 1
#SBATCH -t 60:00:00
#SBATCH --mail-type ALL
#SBATCH --mail-user andrescb@kth.se
#SBATCH --output ./logs/shap.%j.out
#SBATCH --error  ./logs/shap.%j.error
#SBATCH --account=Berzelius-2024-98

cd ../
singularity run --nv ../../../../tensorflow-2.11.0.sif python main_shap_200.py
