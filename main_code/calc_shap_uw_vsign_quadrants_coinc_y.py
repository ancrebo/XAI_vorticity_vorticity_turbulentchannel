# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------------------------------------------------------------------
calc_shap_uw_vsign_quadrants_coinc_y.py
-------------------------------------------------------------------------------------------------------------------------
Created on Tue Jun 18 11:30:49 2024

@author: Andres Cremades Botella

Function to calculate the volume of coincidence of the different shap streaks based on shap_v structures
as a function of y:
    - folder_def  : (str) name of the folder containing the files for configuring the case of analysis.
    - chd_str     : (str) name of the file containing the data of the channel.
    - folders_str : (str) name of the file containing the folders and files used in the problem.
    - st_data_str : (str) name of the file containing the information required for the statistics.
"""
# -----------------------------------------------------------------------------------------------------------------------
# Define the names of the files containing the definitios of the parameters
# - folder_def : folder containing the files with the definitions required in the problem
# - chd_str    : file containing the data of the channel
# - folders    : file containing the folder and file structures
# - st_data    : file containing the data of the statistics
# -----------------------------------------------------------------------------------------------------------------------
folder_def  = "d20240603_definitions"
chd_str     = "channel_data"
folders_str = "folders"
st_data_str = "stats_data_shap"
sh_data_str = "shap_data"
tr_data_str = "training_data"

# -----------------------------------------------------------------------------------------------------------------------
# Import packages
# -----------------------------------------------------------------------------------------------------------------------
import os
from py_bin.py_functions.calc_coinc_shap_uw_vsign import calc_coinc,save_coinc
import numpy as np
from py_bin.py_class.shap_uw_vsign_structure import shap_structure

# -----------------------------------------------------------------------------------------------------------------------
# Unlock the h5 files for avoiding problems in some clusters
# -----------------------------------------------------------------------------------------------------------------------
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

# -----------------------------------------------------------------------------------------------------------------------
# Import information files
# -----------------------------------------------------------------------------------------------------------------------
exec("from "+folder_def+" import "+chd_str+" as chd")
exec("from "+folder_def+" import "+folders_str+" as folders")
exec("from "+folder_def+" import "+st_data_str+" as st_data")
exec("from "+folder_def+" import "+sh_data_str+" as sh_data")
exec("from "+folder_def+" import "+tr_data_str+" as tr_data")


# -----------------------------------------------------------------------------------------------------------------------
# Data for the statistics:
#     - index         : index of the field
#     - Hperc         : percolation index
#     - uvw_folder    : folder of the flow field data
#     - uvw_file      : file of the flow field data
#     - umean_file    : file to save the mean velocity
#     - data_folder   : folder to store the calculated data
#     - dx            : downsampling in x
#     - dy            : downsampling in y
#     - dz            : downsampling in z
#     - L_x           : length of the channel in the streamwise direction
#     - L_y           : half-width of the channel in the wall-normal direction
#     - L_z           : length of the channel in the spanwise direction
#     - urms_file     : file to save the rms of the velocity
#     - rey           : Friction Reynolds number
#     - utau          : Friction velocity
#     - padding       : padding of the flow field
#     - sym_quad      : flag for using the symmetry in the direction 2 of the field for the quadrant selection
#     - filvol        : volume for filtering the structures+
#     - shap_folder   : folder of the shap values
#     - shap_folder   : file of the shap values
#     - hunt_folder   : folder of the hunt structures
#     - hunt_file     : file of the hunt structures
#     - padding       : padding of the field
#     - data_type     : type of data used by the model
#     - SHAPq_folder  : folder of the shap structures
#     - SHAPq_file    : file of the uv structures
#     - nsamples      : number of samples of the shap calculation
#     - SHAPrms_file  : file of the rms of the shap
# -----------------------------------------------------------------------------------------------------------------------
index                   = 7000
Hperc                   = 0.83
uvw_folder              = folders.uvw_folder
uvw_file                = folders.uvw_file
umean_file              = folders.umean_file
data_folder             = folders.data_folder
dx                      = chd.dx
dy                      = chd.dy
dz                      = chd.dz
L_x                     = chd.L_x
L_y                     = chd.L_y
L_z                     = chd.L_z
urms_file               = folders.urms_file
rey                     = chd.rey
utau                    = chd.utau
padding                 = chd.padding
sym_quad                = True
filvol                  = chd.filvol
padding                 = chd.padding
data_type               = tr_data.data_type
plot_folder             = folders.plot_folder
shap_folder             = folders.shap_folder
shap_file               = folders.shap_file
SHAPq_uw_vsign_folder   = folders.SHAPq_uw_vsign_folder
SHAPq_uw_vsign_file     = folders.SHAPq_uw_vsign_file
nsamples                = sh_data.nsamples
SHAPrms_file            = folders.SHAPrms_file
SHAPmean_file           = folders.SHAPmean_file
hunt_shap_file          = folders.hunt_shap_file
umax_file               = folders.umax_file
calc_coin_shap_uw_vsign = folders.calc_coin_shap_uw_vsign

# -----------------------------------------------------------------------------------------------------------------------
# Create the data of the shap structures
# -----------------------------------------------------------------------------------------------------------------------
shap_data  = {"uvw_folder":uvw_folder,"uvw_file":uvw_file,"Hperc":Hperc,"index":index,"dx":dx,
              "dy":dy,"dz":dz,"L_x":L_x,"L_y":L_y,"L_z":L_z,"rey":rey,"utau":utau,
              "padding":padding,"data_folder":data_folder,"umean_file":umean_file,
              "urms_file":urms_file,"sym_quad":True,"filvol":filvol,"shap_folder":shap_folder,
              "shap_file":shap_file,"folder":SHAPq_uw_vsign_folder,"file":SHAPq_uw_vsign_file,
              "padding":padding,"data_type":data_type,"nsamples":nsamples,
              "SHAPrms_file":SHAPrms_file,"SHAPmean_file":SHAPmean_file}

# -----------------------------------------------------------------------------------------------------------------------
# calculate the coincidence between all the structures as a function of y
# -----------------------------------------------------------------------------------------------------------------------


shap_struc            = shap_structure(data_in=shap_data)
shap_struc.read_struc()
data_out = calc_coinc(data_in={"data_struc":shap_struc,"save_data":True,
                               "calc_coin_file":calc_coin_shap_uw_vsign,"folder":data_folder,
                               "dy":dy,"dx":dx,"dz":dz,"uvw_folder":uvw_folder,
                               "uvw_file":uvw_file,"L_x":L_x,"L_y":L_y,"L_z":L_z,
                               "rey":rey,"utau":utau})    
 