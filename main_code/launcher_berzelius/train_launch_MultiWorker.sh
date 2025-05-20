#!/bin/bash
#SBATCH --gpus 16
#SBATCH -N 2
#SBATCH -t 60:00:00
#SBATCH --mail-type ALL
#SBATCH --mail-user andrescb@kth.se
#SBATCH --output ./train.out
#SBATCH --error  ./train.error
#SBATCH --account=Berzelius-2024-98

cd ../
module load buildenv-gcccuda/11.4-8.3.1-bare
srun apptainer run --nv ../../tensorflow-24.03-tf2-py3.sif python main_CNN.py