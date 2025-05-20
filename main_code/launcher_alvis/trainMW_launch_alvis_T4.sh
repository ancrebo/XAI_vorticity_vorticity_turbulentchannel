#!/usr/bin/env bash
#SBATCH --job-name=MW_T4_2
#SBATCH -A NAISS2024-5-129 -p alvis
#SBATCH --nodes 2
#SBATCH --ntasks-per-node=8
#SBATCH --gpus-per-node=T4:8
#SBATCH --cpus-per-task=4
#SBATCH --hint=nomultithread
#SBATCH --distribution=block:block
#SBATCH --time=10:00:00
#SBATCH --mail-type ALL
#SBATCH --mail-user andrescb@kth.se
#SBATCH --output ./logs/trainT4.out
#SBATCH --error  ./logs/trainT4.error


unset http_proxy https_proxy HTTP_PROXY HTTPS_PROXY
set -x
cd ../


module purge
module load  TensorFlow/2.6.0-foss-2021a-CUDA-11.3.1

srun python main_CNN.py
