# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------------------------------------------------------------------
calc_streak_chong_coinc_y.py
-------------------------------------------------------------------------------------------------------------------------
Created on Tue Jun 18 11:30:49 2024

@author: Andres Cremades Botella

Function to calculate the volume occupied by the streaks the shap and the coincidence as a function of y:
    - folder_def  : (str) name of the folder containing the files for configuring the case of analysis.
    - chd_str     : (str) name of the file containing the data of the channel.
    - folders_str : (str) name of the file containing the folders and files used in the problem.
    - st_data_str : (str) name of the file containing the information required for the statistics.
For more information about the tangential Reynolds stress structures:
    - Lozano-Durán, A., Flores, O., & Jiménez, J. (2012). The three-dimensional structure of momentum transfer in
      turbulent channels. Journal of Fluid Mechanics, 694, 100-130.
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
from py_bin.py_class.streak_structure import streak_structure
from py_bin.py_class.chong_structure import chong_structure
from py_bin.py_class.uv_structure import uv_structure
from py_bin.py_class.shap_structure import shap_structure
import os
from py_bin.py_functions.calc_coinc import calc_coinc_4struc,save_coinc_4struc

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
#     - index            : index of the field
#     - Hperc            : percolation index
#     - uvw_folder       : folder of the flow field data
#     - uvw_file         : file of the flow field data
#     - umean_file       : file to save the mean velocity
#     - data_folder      : folder to store the calculated data
#     - dx               : downsampling in x
#     - dy               : downsampling in y
#     - dz               : downsampling in z
#     - L_x              : length of the channel in the streamwise direction
#     - L_y              : half-width of the channel in the wall-normal direction
#     - L_z              : length of the channel in the spanwise direction
#     - urms_file        : file to save the rms of the velocity
#     - rey              : Friction Reynolds number
#     - utau             : Friction velocity
#     - padding          : padding of the flow field
#     - sym_quad         : flag for using the symmetry in the direction 2 of the field for the quadrant selection
#     - filvol           : volume for filtering the structures+
#     - shap_folder      : folder of the shap values
#     - shap_folder      : file of the shap values
#     - uv_folder        : folder of the streak structures
#     - uv_file          : file of the streak structures
#     - streak_folder    : folder of the streak structures
#     - streak_file      : file of the streak structures
#     - chong_folder     : folder of the chong structures
#     - chong_file       : file of the chong structures
#     - padding          : padding of the field
#     - data_type        : type of data used by the model
#     - SHAPq_folder     : folder of the shap structures
#     - SHAPq_file       : file of the shap structures
#     - nsamples         : number of samples of the shap calculation
#     - SHAPrms_file     : file of the rms of the shap
#     - streak_shap_file : file for saving the coincidence between streak vortices and shap structures
# -----------------------------------------------------------------------------------------------------------------------
index_ini                 = st_data.field_ini
index_fin                 = st_data.field_fin
index_delta               = st_data.field_delta*100
Hperc                     = 1.41
uvw_folder                = folders.uvw_folder
uvw_file                  = folders.uvw_file
umean_file                = folders.umean_file
data_folder               = folders.data_folder
dx                        = chd.dx
dy                        = chd.dy
dz                        = chd.dz
L_x                       = chd.L_x
L_y                       = chd.L_y
L_z                       = chd.L_z
urms_file                 = folders.urms_file
rey                       = chd.rey
utau                      = chd.utau
padding                   = chd.padding
sym_quad                  = True
filvol                    = chd.filvol
shap_folder               = folders.shap_folder
shap_file                 = folders.shap_file
streak_folder             = folders.streak_folder
streak_file               = folders.streak_file
padding                   = chd.padding
data_type                 = tr_data.data_type
plot_folder               = folders.plot_folder
uv_folder                 = folders.uv_folder
uv_file                   = folders.uv_file
chong_folder              = folders.chong_folder
chong_file                = folders.chong_file
SHAPq_folder              = folders.SHAPq_folder
SHAPq_file                = folders.SHAPq_file
nsamples                  = sh_data.nsamples
SHAPrms_file              = folders.SHAPrms_file
uv_chong_streak_shap_file = folders.uv_chong_streak_shap_file

# -----------------------------------------------------------------------------------------------------------------------
# Create the data of the streak structure
# -----------------------------------------------------------------------------------------------------------------------
data_streak  = {"uvw_folder":uvw_folder,"uvw_file":uvw_file,"Hperc":Hperc,"index":0,"dx":dx,
                "dy":dy,"dz":dz,"L_x":L_x,"L_y":L_y,"L_z":L_z,"rey":rey,"utau":utau,
                "padding":padding,"data_folder":data_folder,"umean_file":umean_file,
                "urms_file":urms_file,"sym_quad":True,"filvol":filvol,"shap_folder":shap_folder,
                "shap_file":shap_file,"folder":streak_folder,"file":streak_file,"padding":padding,
                "data_type":data_type}

# -----------------------------------------------------------------------------------------------------------------------
# Create the data of the chong structure
# -----------------------------------------------------------------------------------------------------------------------
data_chong = {"uvw_folder":uvw_folder,"uvw_file":uvw_file,"Hperc":Hperc,"index":0,"dx":dx,
              "dy":dy,"dz":dz,"L_x":L_x,"L_y":L_y,"L_z":L_z,"rey":rey,"utau":utau,
              "padding":padding,"data_folder":data_folder,"umean_file":umean_file,
              "urms_file":urms_file,"sym_quad":True,"filvol":filvol,"shap_folder":shap_folder,
              "shap_file":shap_file,"folder":chong_folder,"file":chong_file,"padding":padding,
              "data_type":data_type}

# -----------------------------------------------------------------------------------------------------------------------
# Create the data of the uv structure
# -----------------------------------------------------------------------------------------------------------------------
data_uv = {"uvw_folder":uvw_folder,"uvw_file":uvw_file,"Hperc":Hperc,"index":0,"dx":dx,
           "dy":dy,"dz":dz,"L_x":L_x,"L_y":L_y,"L_z":L_z,"rey":rey,"utau":utau,
           "padding":padding,"data_folder":data_folder,"umean_file":umean_file,
           "urms_file":urms_file,"sym_quad":True,"filvol":filvol,"shap_folder":shap_folder,
           "shap_file":shap_file,"folder":uv_folder,"file":uv_file,"padding":padding,
           "data_type":data_type}
           
# -----------------------------------------------------------------------------------------------------------------------
# Create the data of the shap structure
# -----------------------------------------------------------------------------------------------------------------------
data_shap = {"uvw_folder":uvw_folder,"uvw_file":uvw_file,"Hperc":Hperc,"index":0,"dx":dx,
             "dy":dy,"dz":dz,"L_x":L_x,"L_y":L_y,"L_z":L_z,"rey":rey,"utau":utau,
             "padding":padding,"data_folder":data_folder,"umean_file":umean_file,
             "urms_file":urms_file,"sym_quad":True,"filvol":filvol,"shap_folder":shap_folder,
             "shap_file":shap_file,"folder":SHAPq_folder,"file":SHAPq_file,"padding":padding,
             "data_type":data_type,"nsamples":nsamples,"SHAPrms_file":SHAPrms_file}

# -----------------------------------------------------------------------------------------------------------------------
# calculate the coincidence between the chong and the streaks structures as a function of y
# -----------------------------------------------------------------------------------------------------------------------
index_range = range(index_ini,index_fin,index_delta)
for ii in index_range:
    print(ii,flush=True)
    data_streak["index"] = ii
    data_chong["index"]  = ii
    data_uv["index"]     = ii
    data_shap["index"]   = ii
    streak_struc         = streak_structure(data_in=data_streak)
    streak_struc.read_struc()
    chong_struc          = chong_structure(data_in=data_chong)
    chong_struc.read_struc()
    uv_struc             = uv_structure(data_in=data_uv)
    uv_struc.read_struc()
    shap_struc           = shap_structure(data_in=data_shap)
    shap_struc.read_struc()
    data_out             = calc_coinc_4struc(data_in={"data_struc1":uv_struc,"data_struc2":streak_struc,
                                                      "data_struc3":chong_struc,"data_struc4":shap_struc,
                                                      "save_data":False,"calc_coin_file":uv_chong_streak_shap_file,
                                                      "folder":data_folder,"dy":dy,"dx":dx,"dz":dz,
                                                      "uvw_folder":uvw_folder,"uvw_file":uvw_file,"L_x":L_x,
                                                      "L_y":L_y,"L_z":L_z,"rey":rey,"utau":utau})
    if ii == index_range[0]:
        frac_struc1  = data_out["frac_struc1"]
        frac_struc2  = data_out["frac_struc2"]
        frac_struc3  = data_out["frac_struc3"]
        frac_struc4  = data_out["frac_struc4"]
        frac_coinc   = data_out["frac_coinc"]
        yplus        = data_out["yplus"]
    else:
        frac_struc1 += data_out["frac_struc1"]
        frac_struc2 += data_out["frac_struc2"]
        frac_struc3 += data_out["frac_struc3"]
        frac_struc4 += data_out["frac_struc4"]
        frac_coinc  += data_out["frac_coinc"]
frac_struc1 /= len(index_range)
frac_struc2 /= len(index_range)
frac_coinc  /= len(index_range)
save_coinc_4struc(data_in={"frac_struc1":frac_struc1,"frac_struc2":frac_struc2,"frac_struc3":frac_struc3,
                           "frac_struc4":frac_struc4,"frac_coinc":frac_coinc,
                           "yplus":yplus,"calc_coin_file":uv_chong_streak_shap_file,"folder":data_folder})