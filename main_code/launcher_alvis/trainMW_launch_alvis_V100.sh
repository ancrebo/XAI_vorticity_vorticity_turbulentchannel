#!/usr/bin/env bash
#SBATCH --job-name=MW_V100_2
#SBATCH -A NAISS2024-5-129 -p alvis
#SBATCH --nodes 2
#SBATCH --ntasks-per-node=2
#SBATCH --gpus-per-node=V100:2
#SBATCH --cpus-per-task=8
#SBATCH --hint=nomultithread
#SBATCH --distribution=block:block
#SBATCH --time=10:00:00
#SBATCH --mail-type ALL
#SBATCH --mail-user andrescb@kth.se
#SBATCH --output ./logs/trainV100.out
#SBATCH --error  ./logs/trainV100.error


unset http_proxy https_proxy HTTP_PROXY HTTPS_PROXY
set -x
cd ../


module purge
module load  TensorFlow/2.6.0-foss-2021a-CUDA-11.3.1

srun python main_CNN.py
