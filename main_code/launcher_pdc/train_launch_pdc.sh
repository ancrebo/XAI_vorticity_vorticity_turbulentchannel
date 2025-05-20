#!/usr/bin/env bash
#SBATCH --job-name=MW_A100_1
#SBATCH -A naiss2023-3-13 -p gpu
#SBATCH --nodes 1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1:00:00
#SBATCH --mail-type ALL
#SBATCH --mail-user andrescb@kth.se
#SBATCH --output ./logs/train.out
#SBATCH --error  ./logs/train.error


unset http_proxy https_proxy HTTP_PROXY HTTPS_PROXY
set -x
cd ../


module purge
ml add PDC/22.06
ml add singularity/3.10.4-cpeGNU-22.06

singularity exec --rocm -B /cfs/klemming /pdc/software/resources/sing_hub/rocm5.4-tf2.10 python3 main_CNN.py