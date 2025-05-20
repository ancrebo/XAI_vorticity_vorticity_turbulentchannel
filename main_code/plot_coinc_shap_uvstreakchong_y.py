# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------------------------------------------------------------------
plot_uv_shap_coinc_y.py
-------------------------------------------------------------------------------------------------------------------------
Created on Thu Mar 28 13:38:39 2024

@author: Andres Cremades Botella

Create the plot of the mean velocities. The file requires to set the following paths:
    - folder_def  : (str) name of the folder containing the files for configuring the case of analysis.
    - chd_str     : (str) name of the file containing the data of the channel.
    - folders_str : (str) name of the file containing the folders and files used in the problem.
    - st_data_str : (str) name of the file containing the information required for the statistics.
In addition the following variables need to be set:
    - xlabel      : label of the x axis
    - ylabel      : label of the y axis
    - fontsize    : size of the text in the figure
    - figsize_x   : size of the figure in axis x
    - figsize_y   : size of the figure in axis y
    - colormap    : colormap used in the plot
    - colornum    : number of levels required in the colormap
    - fig_name    : name of the figure after saving
    - dpi         : dots per inch of the figure
"""
# ----------------------------------------------------------------------------------------------------------------------
# Define the names of the files containing the definitios of the parameters
# - folder_def : folder containing the files with the definitions required in the problem
# - folders    : file containing the folder and file structures
# ----------------------------------------------------------------------------------------------------------------------
folder_def  = "d20240603_definitions"
chd_str     = "channel_data"
folders_str = "folders"
st_data_str = "stats_data"

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
xlabel      = "$y^+$"
ylabel      = "$V/V_{tot}$"
ylabel2     = "$Coincidence (\%)$"
fontsize    = 18
figsize_x   = 10
figsize_y   = 5
colormap    = "viridis"
colornum    = 4
fig_name    = "Coincidence_shap_uvstreakchong"
dpi         = 200
struc1a_lab = "SHAP"
struc1b_lab = "uv"
struc2a_lab = "SHAP"
struc2b_lab = "streak"
struc3a_lab = "SHAP"
struc3b_lab = "chong"

# -----------------------------------------------------------------------------------------------------------------------
# Import packages
# -----------------------------------------------------------------------------------------------------------------------
from py_bin.py_plots.plot_coinc import plot_coinc_3coinc
import os

# ----------------------------------------------------------------------------------------------------------------------
# Unlock the h5 files for avoiding problems in some clusters
# ----------------------------------------------------------------------------------------------------------------------
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

# ----------------------------------------------------------------------------------------------------------------------
# Import information files
# ----------------------------------------------------------------------------------------------------------------------
exec("from "+folder_def+" import "+chd_str+" as chd")
exec("from "+folder_def+" import "+folders_str+" as folders")
exec("from "+folder_def+" import "+st_data_str+" as st_data")

# -----------------------------------------------------------------------------------------------------------------------
# Define the information of the data
#     - data_folder : folder containing the training data
#     - plot_folder : folder to save the figures
#     - file        : name of the file containing the mean velocity information
#     - dy          : downsampling of the wall-normal direction
#     - dx          : downsampling of the streamwise direction
#     - dz          : downsampling of the spanwise direction
#     - uvw_folder  : folder of the flow fields
#     - uvw_file    : file of the flow fields
#     - L_x         : dimension of the channel in the streamwise direction
#     - L_y         : dimension of the channel in the wall-normal direction
#     - L_z         : dimension of the channel in the spanwise direction
#     - file_trj    : file containing the data of Torroja
# -----------------------------------------------------------------------------------------------------------------------
data_folder = folders.data_folder
plot_folder = folders.plot_folder
file_uv     = folders.uv_shap_file
file_streak = folders.streak_shap_file
file_chong  = folders.chong_shap_file

# -----------------------------------------------------------------------------------------------------------------------
# Create the plot
# -----------------------------------------------------------------------------------------------------------------------
plot_format_data = {"file_1":file_uv,"file_2":file_streak,"file_3":file_chong,"folder":data_folder,
                    "plot_folder":plot_folder,"xlabel":xlabel,"ylabel":ylabel,"ylabel2":ylabel2,
                    "fontsize":fontsize,"figsize_x":figsize_x,"figsize_y":figsize_y,"colormap":colormap,
                    "colornum":colornum,"fig_name":fig_name,"dpi":dpi,"struc1a_lab":struc1a_lab,
                    "struc1b_lab":struc1b_lab,"struc2a_lab":struc2a_lab,"struc2b_lab":struc2b_lab,
                    "struc3a_lab":struc3a_lab,"struc3b_lab":struc3b_lab}
plot_coinc_3coinc(data_in=plot_format_data)