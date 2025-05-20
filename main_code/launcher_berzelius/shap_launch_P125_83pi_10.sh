#!/bin/bash
#SBATCH --gpus 1
#SBATCH -J shap_10
#SBATCH -N 1
#SBATCH -t 60:00:00
#SBATCH --mail-type ALL
#SBATCH --mail-user andrescb@kth.se
#SBATCH --output ./logs/shap.%j.out
#SBATCH --error  ./logs/shap.%j.error
#SBATCH --account=Berzelius-2024-98

cd ../
module load Mambaforge/23.3.1-1-hpc1-bdist
mamba activate tensorflow_2.17.0
python main_shap_10.py
