#!/usr/bin/env bash
#SBATCH -A NAISS2024-5-129 -p alvis
#SBATCH -N 2 --gpus-per-node=A40:4
#SBATCH -t 10:00:00
#SBATCH --mail-type ALL
#SBATCH --mail-user andrescb@kth.se
#SBATCH --output ./logs/train.out
#SBATCH --error  ./logs/train.error

cd ../
module load  TensorFlow/2.11.0-foss-2022a-CUDA-11.7.0
python main_CNN.py
