#!/bin/bash
#SBATCH --gpus 1
#SBATCH -N 1
#SBATCH -t 1:00:00
#SBATCH --mail-type ALL
#SBATCH --mail-user andrescb@kth.se
#SBATCH --output ./launcher_berzelius/logs/shap550_%j.out
#SBATCH --error  ./launcher_berzelius/logs/shap550_%j.error
#SBATCH --account=Berzelius-2024-450

fini=$1
ffin=$2
timewait=$3

sleep $timewait

echo "$fini $ffin"

echo "hola"

nvidia-smi


cd /proj/wasp_wise/users/x_andcr/SHAP/XAI_TurbulentChannel_optimized-main/main_code

foldername="folders_${fini}_${ffin}"
shapname="shap_data_${fini}_${ffin}"
mainname="main_SHAP_params_${fini}_${ffin}.py"

cp  P550_21pi_250225_v2_definitions/folders_base.py  P550_21pi_250225_v2_definitions/"$foldername".py
cp  P550_21pi_250225_v2_definitions/shap_data_base.py  P550_21pi_250225_v2_definitions/"$shapname".py
cp  main_SHAP_params.py  "$mainname"

sed -i "s/%INDW%/0001/" P550_21pi_250225_v2_definitions/"$foldername".py
sed -i "s/%INDR%/0001/" P550_21pi_250225_v2_definitions/"$foldername".py
sed -i "s/%FINI%/$fini/g" P550_21pi_250225_v2_definitions/"$shapname".py
sed -i "s/%FFIN%/$ffin/g" P550_21pi_250225_v2_definitions/"$shapname".py

sed -i "s/%FOLDER%/$foldername/g"  "$mainname"
sed -i "s/%SHAP%/$shapname/g" "$mainname"

singularity run --nv ../../tensorflow-2.9.1.sif python3 "$mainname"

rm P550_21pi_250225_v2_definitions/"$foldername".py
rm P550_21pi_250225_v2_definitions/"$shapname".py
rm "$mainname"
