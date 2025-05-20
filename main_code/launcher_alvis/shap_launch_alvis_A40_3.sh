#!/usr/bin/env bash
#SBATCH --job-name=A40_3
#SBATCH -A NAISS2024-5-129 -p alvis
#SBATCH --nodes 1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus-per-node=A40:2
#SBATCH --hint=nomultithread
#SBATCH --distribution=block:block
#SBATCH --time=60:00:00
#SBATCH --mail-type ALL
#SBATCH --mail-user andrescb@kth.se
#SBATCH --output ./logs/shapA40_3.out
#SBATCH --error  ./logs/shapA40_3.error


unset http_proxy https_proxy HTTP_PROXY HTTPS_PROXY
set -x
cd ../


module purge
module load  TensorFlow/2.7.1-foss-2021b-CUDA-11.4.1
module load  scikit-learn/1.0.1-foss-2021b
module load  numba/0.54.1-foss-2021b-CUDA-11.4.1
module load  tqdm/4.62.3-GCCcore-11.2.0
python main_SHAP_3.py