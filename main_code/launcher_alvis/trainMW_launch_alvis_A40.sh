#!/usr/bin/env bash
#SBATCH --job-name=MW_A40_2
#SBATCH -A NAISS2024-5-129 -p alvis
#SBATCH --nodes 2
#SBATCH --ntasks-per-node=1
#SBATCH --gpus-per-node=A40:4
#SBATCH --cpus-per-task=64
#SBATCH --hint=nomultithread
#SBATCH --distribution=block:block
#SBATCH --time=10:00:00
#SBATCH --mail-type ALL
#SBATCH --mail-user andrescb@kth.se
#SBATCH --output ./logs/trainA40.out
#SBATCH --error  ./logs/trainA40.error


unset http_proxy https_proxy HTTP_PROXY HTTPS_PROXY
set -x
cd ../


module purge
module load  TensorFlow/2.6.0-foss-2021a-CUDA-11.3.1

srun python main_CNN.py
