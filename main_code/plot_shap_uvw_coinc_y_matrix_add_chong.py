# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------------------------------------------------------------------
plot_shap_3d.py
-------------------------------------------------------------------------------------------------------------------------
Created on Tue Jun 18 11:30:49 2024

@author: Andres Cremades Botella

Function to plot the 3d Reynolds stress structures:
    - folder_def  : (str) name of the folder containing the files for configuring the case of analysis.
    - chd_str     : (str) name of the file containing the data of the channel.
    - folders_str : (str) name of the file containing the folders and files used in the problem.
    - st_data_str : (str) name of the file containing the information required for the statistics.
"""
# -----------------------------------------------------------------------------------------------------------------------
# Define the variables required for the plot
#     - xlabel      : label of the x axis
#     - ylabel      : label of the y axis
#     - fontsize    : size of the text in the figure
#     - figsize_x   : size of the figure in axis x
#     - figsize_y   : size of the figure in axis y
#     - colormap    : colormap used in the plot
#     - colornum    : number of levels required in the colormap
#     - fig_name    : name of the figure after saving
#     - dpi         : dots per inch of the figure
# -----------------------------------------------------------------------------------------------------------------------
xlabel           = "$x^+$"
ylabel           = "$z^+$"
fontsize         = 18
figsize_x        = 7
figsize_y        = 5
colormap         = "viridis"
fig_name         = "shap_uvw_quadrants_coinc_mat_add_chong"
plot_coin_folder = "shap_uvw_quadrants_coinc_mat_add_chong"
dpi              = 200
colors_11        = ['#1f77b4',  # Blue
                    '#ff7f0e',  # Orange
                    '#2ca02c',  # Green
                    '#d62728',  # Red
                    '#9467bd',  # Purple
                    '#8c564b',  # Brown
                    '#e377c2',  # Pink
                    '#7f7f7f',  # Gray
                    '#bcbd22',  # Olive
                    '#17becf',  # Cyan
                    '#aec7e8'   # Light blue
                    ]

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
st_data_str = "stats_data"
sh_data_str = "shap_data"
tr_data_str = "training_data"

# -----------------------------------------------------------------------------------------------------------------------
# Import packages
# -----------------------------------------------------------------------------------------------------------------------
from py_bin.py_class.chong_structure import chong_structure
from py_bin.py_class.shap_uvw_structure import shap_structure
import os
from py_bin.py_plots.plot_coinc_mat_quadrants_uvw import plot_coinc_mat_11color_withcontour

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
#     - index        : index of the field
#     - Hperc        : percolation index
#     - uvw_folder   : folder of the flow field data
#     - uvw_file     : file of the flow field data
#     - umean_file   : file to save the mean velocity
#     - data_folder  : folder to store the calculated data
#     - dx           : downsampling in x
#     - dy           : downsampling in y
#     - dz           : downsampling in z
#     - L_x          : length of the channel in the streamwise direction
#     - L_y          : half-width of the channel in the wall-normal direction
#     - L_z          : length of the channel in the spanwise direction
#     - urms_file    : file to save the rms of the velocity
#     - rey          : Friction Reynolds number
#     - utau         : Friction velocity
#     - padding      : padding of the flow field
#     - sym_quad     : flag for using the symmetry in the direction 2 of the field for the quadrant selection
#     - filvol       : volume for filtering the structures+
#     - shap_folder  : folder of the shap values
#     - shap_folder  : file of the shap values
#     - SHAPq_folder : folder of the shap structures
#     - SHAPq_file   : file of the uv structures
#     - padding      : padding of the field
#     - data_type    : type of data used by the model
#     - nsamples     : number of samples of the shap calculation
#     - SHAPrms_file : file of the rms of the shap
# -----------------------------------------------------------------------------------------------------------------------
index              = 7000
Hperc              = 1.18
uvw_folder         = folders.uvw_folder
uvw_file           = folders.uvw_file
umean_file         = folders.umean_file
data_folder        = folders.data_folder
dx                 = chd.dx
dy                 = chd.dy
dz                 = chd.dz
L_x                = chd.L_x
L_y                = chd.L_y
L_z                = chd.L_z
urms_file          = folders.urms_file
rey                = chd.rey
utau               = chd.utau
padding            = chd.padding
sym_quad           = True
filvol             = chd.filvol
shap_folder        = folders.shap_folder
shap_file          = folders.shap_file
SHAPq_uvw_folder   = folders.SHAPq_uvw_folder
SHAPq_uvw_file     = folders.SHAPq_uvw_file
chong_folder       = folders.chong_folder
chong_file         = folders.chong_file
padding            = chd.padding
data_type          = tr_data.data_type
plot_folder        = folders.plot_folder
nsamples           = sh_data.nsamples
SHAPrms_file       = folders.SHAPrms_file
SHAPmean_file      = folders.SHAPmean_file
calc_coin_shap_uvw = folders.calc_coin_shap_uvw
data_folder        = folders.data_folder

# -----------------------------------------------------------------------------------------------------------------------
# Create the data of the uv structure
# -----------------------------------------------------------------------------------------------------------------------
chong_struc = chong_structure(data_in={"uvw_folder":uvw_folder,"uvw_file":uvw_file,"Hperc":Hperc,"index":index,
                                         "dx":dx,"dy":dy,"dz":dz,"L_x":L_x,"L_y":L_y,"L_z":L_z,"rey":rey,"utau":utau,
                                         "padding":padding,"data_folder":data_folder,"umean_file":umean_file,
                                         "urms_file":urms_file,"sym_quad":True,"filvol":filvol,"shap_folder":shap_folder,
                                         "shap_file":shap_file,"folder":chong_folder,"file":chong_file,
                                         "padding":padding,"data_type":data_type})
chong_struc.read_struc()

# -----------------------------------------------------------------------------------------------------------------------
# Create the data of the structure
# -----------------------------------------------------------------------------------------------------------------------
shap_struc = shap_structure(data_in={"uvw_folder":uvw_folder,"uvw_file":uvw_file,"Hperc":Hperc,"index":index,"dx":dx,
                                     "dy":dy,"dz":dz,"L_x":L_x,"L_y":L_y,"L_z":L_z,"rey":rey,"utau":utau,
                                     "padding":padding,"data_folder":data_folder,"umean_file":umean_file,
                                     "urms_file":urms_file,"sym_quad":True,"filvol":filvol,"shap_folder":shap_folder,
                                     "shap_file":shap_file,"folder":SHAPq_uvw_folder,"file":SHAPq_uvw_file,
                                     "padding":padding,"data_type":data_type,"nsamples":nsamples,
                                     "SHAPrms_file":SHAPrms_file,"SHAPmean_file":SHAPmean_file})
shap_struc.read_struc()
shap_struc.read_struc_Q()

# -----------------------------------------------------------------------------------------------------------------------
# Plot the 3D field
# -----------------------------------------------------------------------------------------------------------------------
plot_coinc_mat_11color_withcontour(data_in={"struc":shap_struc,"plot_folder":plot_folder,"xlabel":xlabel,"ylabel":ylabel,
                                            "fontsize":fontsize,"figsize_x":figsize_x,"figsize_y":figsize_y,
                                            "colors_11":colors_11,"fig_name":fig_name,"dpi":dpi,"dy":dy,"dx":dx,"dz":dz,
                                            "uvw_folder":uvw_folder,"uvw_file":uvw_file,"L_x":L_x,"L_y":L_y,"L_z":L_z,
                                            "rey":rey,"utau":utau,"cmap_flag":False,"index_ii":index,
                                            "plot_coin_folder":plot_coin_folder,"calc_coin_shap_uvw":calc_coin_shap_uvw,
                                            "data_folder":data_folder,"cont_struc":chong_struc})