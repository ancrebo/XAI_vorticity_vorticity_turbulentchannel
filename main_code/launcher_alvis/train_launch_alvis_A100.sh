#!/usr/bin/env bash
#SBATCH --job-name=MW_A100_1
#SBATCH -A NAISS2024-5-129 -p alvis
#SBATCH --nodes 1
#SBATCH --exclude=alvis[3-4]-[01-09]
#SBATCH --ntasks-per-node=1
#SBATCH --gpus-per-node=A100:4
#SBATCH --cpus-per-task=64
#SBATCH --hint=nomultithread
#SBATCH --distribution=block:block
#SBATCH --time=60:00:00
#SBATCH --mail-type ALL
#SBATCH --mail-user andrescb@kth.se
#SBATCH --output ./logs/trainA100.out
#SBATCH --error  ./logs/trainA100.error


unset http_proxy https_proxy HTTP_PROXY HTTPS_PROXY
set -x
cd ../


module purge
module load  TensorFlow/2.7.1-foss-2021b-CUDA-11.4.1

srun python main_CNN.py