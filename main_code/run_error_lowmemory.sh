INIFIELD=20000
FINFIELD=27000
source ../../SHAP_enviroment_tf2_17/bin/activate

for ((ii=$INIFIELD; ii<$FINFIELD; ii++))
do
	echo "Field $ii ...$
    cp  cal_pred_error_lowmemory_restart.py  cal_pred_error_lowmemory_restart_bis.py
    ini=$ii
	fin=$ii+1
    sed -i "s/%NAMEFILE%/restart1/" cal_pred_error_lowmemory_restart_bis.py
    sed -i "s/%INI%/$ini/" cal_pred_error_lowmemory_restart_bis.py
    sed -i "s/%FIN%/$fin/" cal_pred_error_lowmemory_restart_bis.py
	python cal_pred_error_lowmemory_restart_bis.py
done