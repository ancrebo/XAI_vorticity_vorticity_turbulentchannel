# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------------------------------------------------------------------
plot_Q_streak_chong_hunt_coinc_y_matrix_cont_shap.py
-------------------------------------------------------------------------------------------------------------------------
Created on Wed Oct 23 11:53:24 2024

@author: Andres Cremades Botella

Function to save the percolation from a h5 file:
    - folder_def  : (str) name of the folder containing the files for configuring the case of analysis.
    - folders_str : (str) name of the file containing the folders and files used in the problem.
"""

import h5py
import numpy as np
from py_bin.py_functions.percolation import percolation,save_percolation

folder_def  = "d20240703_definitions"
folders_str = "folders"
exec("from "+folder_def+" import "+folders_str+" as folders")

perc_SHAP_file = folders.perc_SHAP_file
data_folder    = folders.data_folder

fperc = h5py.File(data_folder+'/percolation_125_8pi3pi.h5','r')
print(fperc.keys())
Hvec = np.array(fperc["H_sh"])
Nvec = np.array(fperc["N_st_coef"])/np.max(np.array(fperc["N_st_coef"]))
Vvec = np.array(fperc["V_st"])/np.array(fperc["V_st_sum"])




save_percolation(data_in={"nstruc":Nvec,"Vstruc":Vvec,"H_perc":Hvec,"perc_file":perc_SHAP_file,
                          "folder":data_folder})