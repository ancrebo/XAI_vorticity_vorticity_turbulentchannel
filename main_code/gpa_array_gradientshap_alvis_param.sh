#!/usr/bin/env bash
#SBATCH --job-name=shapvor
#SBATCH -A NAISS2025-5-144 -p alvis
#SBATCH --nodes 1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus-per-node=A100:1
#SBATCH --hint=nomultithread
#SBATCH --distribution=block:block
#SBATCH --time=4:00:00
#SBATCH --mail-type ALL
#SBATCH --mail-user andrescb@kth.se
#SBATCH --output ./logsalvis/shapvor_%j.out
#SBATCH --error  ./logsalvis/shapvor_%j.error


module purge
module load  TensorFlow/2.7.1-foss-2021b-CUDA-11.4.1
module load  scikit-learn/1.0.1-foss-2021b
module load  numba/0.54.1-foss-2021b-CUDA-11.4.1
module load  tqdm/4.62.3-GCCcore-11.2.0



# Compute field range using SLURM_ARRAY_TASK_ID
initial_field=$1
delta_fields=$2

fini=$(( initial_field + (SLURM_ARRAY_TASK_ID - 1) * delta_fields ))
ffin=$(( fini + delta_fields ))

sleep $timewait

echo "$fini $ffin"

echo "hola"

nvidia-smi


cd /mimer/NOBACKUP/groups/deepmechalvis/andres/SHAP/XAI_vorticity_vorticity_turbulentchannel/main_code

pathname="P125_83pi_250507_v0_definitions"
foldername="folders_${fini}_${ffin}"
shapname="shap_data_${fini}_${ffin}"
mainname="main_SHAP_params_${fini}_${ffin}.py"

cp  "$pathname"/folders_base.py  "$pathname"/"$foldername".py
cp  "$pathname"/shap_data_base.py  "$pathname"/"$shapname".py
cp  main_SHAP_params.py  "$mainname"

sed -i "s/%INDW%/0001/" "$pathname"/"$foldername".py
sed -i "s/%INDR%/0001/" "$pathname"/"$foldername".py
sed -i "s/%FINI%/$fini/g" "$pathname"/"$shapname".py
sed -i "s/%FFIN%/$ffin/g" "$pathname"/"$shapname".py

sed -i "s/%FOLDER%/$foldername/g"  "$mainname"
sed -i "s/%SHAP%/$shapname/g" "$mainname"

srun python3 "$mainname"

rm "$pathname"/"$foldername".py
rm "$pathname"/"$shapname".py
rm "$mainname"

