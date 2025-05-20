#!/bin/bash
#SBATCH --gpus 1 -C "fat"
#SBATCH -J shap_0
#SBATCH -N 1
#SBATCH -t 60:00:00
#SBATCH --mail-type ALL
#SBATCH --mail-user andrescb@kth.se
#SBATCH --output ./logs/shap.%j.out
#SBATCH --error  ./logs/shap.%j.error
#SBATCH --account=Berzelius-2024-98

cd ../
fini=$1
ffin=$2
timewait=$3
sleep $timewait
echo "$fini $ffin"

cp  P125_83pi_240603_v0_definitions/folders_base.py  P125_83pi_240603_v0_definitions/folders.py
cp  P125_83pi_240603_v0_definitions/shap_data_base.py  P125_83pi_240603_v0_definitions/shap_data.py

sed -i "s/%INDW%/5285/" P125_83pi_240603_v0_definitions/folders.py
sed -i "s/%INDR%/5285/" P125_83pi_240603_v0_definitions/folders.py
sed -i "s/%FINI%/$fini/g" P125_83pi_240603_v0_definitions/shap_data.py
sed -i "s/%FFIN%/$ffin/g" P125_83pi_240603_v0_definitions/shap_data.py

singularity run --nv ../../../../tensorflow-2.11.0.sif python main_shap.py
