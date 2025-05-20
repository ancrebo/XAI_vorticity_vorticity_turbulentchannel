# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------------------------------------------------------------------
plot_shapintensity_uv_coinc_y_matrix.py
-------------------------------------------------------------------------------------------------------------------------
Created on Tue Jun 18 11:30:49 2024

@author: Andres Cremades Botella

Function to plot the coincidence between the shap structures and the uv structures based in the intensity
in the whole domain:
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
st_data_str = "stats_data"
sh_data_str = "shap_data"
tr_data_str = "training_data"

# -----------------------------------------------------------------------------------------------------------------------
# Import packages
# -----------------------------------------------------------------------------------------------------------------------
from py_bin.py_class.shap_intensity_structure import shap_structure
from py_bin.py_class.uv_structure import uv_structure
import os
from py_bin.py_functions.calc_coinc import calc_coinc_mat
from py_bin.py_plots.plot_coinc_mat import plot_coinc_mat
from py_bin.py_class.flow_field import flow_field

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
# Define the variables required for the plot
#     - xlabel                   : label of the x axis
#     - ylabel                   : label of the y axis
#     - fontsize                 : size of the text in the figure
#     - figsize_x                : size of the figure in axis x
#     - figsize_y                : size of the figure in axis y
#     - colormap                 : colormap used in the plot
#     - colornum                 : number of levels required in the colormap
#     - fig_name                 : name of the figure after saving
#     - dpi                      : dots per inch of the figure
#     - struc1_lab               : label of the structure 1
#     - struc2_lab               : label of the structure 2
#     - plot_coin_folder_uv_shap : folder to save the coincidence plots between uv and shap structures
#     - plot_coin_file_uv_shap   : folder to save the coincidence plots between uv and shap structures
# -----------------------------------------------------------------------------------------------------------------------
xlabel           = "$x^+$"
ylabel           = "$z^+$"
fontsize         = 18
figsize_x        = 7
figsize_y        = 3
colormap         = "viridis"
colornum         = 4
dpi              = 200
struc1_lab       = "uv"
struc2_lab       = "SHAP intensity"
plot_coin_folder = "shapintensity_uv_coin_y"
plot_coin_file   = "shapintensity_uv_coin_y"

# -----------------------------------------------------------------------------------------------------------------------
# Data for the statistics:
#     - index                   : index of the field
#     - Hperc                   : percolation index
#     - uvw_folder              : folder of the flow field data
#     - uvw_file                : file of the flow field data
#     - umean_file              : file to save the mean velocity
#     - data_folder             : folder to store the calculated data
#     - dx                      : downsampling in x
#     - dy                      : downsampling in y
#     - dz                      : downsampling in z
#     - L_x                     : length of the channel in the streamwise direction
#     - L_y                     : half-width of the channel in the wall-normal direction
#     - L_z                     : length of the channel in the spanwise direction
#     - urms_file               : file to save the rms of the velocity
#     - rey                     : Friction Reynolds number
#     - utau                    : Friction velocity
#     - padding                 : padding of the flow field
#     - sym_quad                : flag for using the symmetry in the direction 2 of the field for the quadrant selection
#     - filvol                  : volume for filtering the structures+
#     - shap_folder             : folder of the shap values
#     - shap_file               : file of the shap values
#     - SHAPq_intensity_folder  : folder of the SHAP values based on the intensity
#     - SHAPq_intensity_file    : file of the SHAP values based on the intensity
#     - padding                 : padding of the field
#     - data_type               : type of data used by the model
#     - plot_folder             : folder to plot the data
#     - uv_folder               : folder of the uv structures
#     - uv_file                 : file of the uv structures
#     - nsamples                : number of samples of the shap calculation
#     - SHAPrms_file            : file of the rms of the shap
#     - SHAPmean_file           : file of the mean of the shap
#     - shapintensity_shap_file : file to store the coincidence between the shap values based in the intensity and
#                                 based on the absolute value
#     - umax_file               : file with the maximum velocity
# -----------------------------------------------------------------------------------------------------------------------
index                   = 7000
Hperc                   = 1.41
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
shap_folder             = folders.shap_folder
shap_file               = folders.shap_file
SHAPq_intensity_folder  = folders.SHAPq_intensity_folder
SHAPq_intensity_file    = folders.SHAPq_intensity_file
padding                 = chd.padding
data_type               = tr_data.data_type
plot_folder             = folders.plot_folder
uv_folder               = folders.uv_folder
uv_file                 = folders.uv_file
nsamples                = sh_data.nsamples
SHAPrms_file            = folders.SHAPrms_file
SHAPmean_file           = folders.SHAPmean_file
umax_file               = folders.umax_file

# -----------------------------------------------------------------------------------------------------------------------
# Create the data of the SHAP intensity structure
# -----------------------------------------------------------------------------------------------------------------------
shap_struc = shap_structure(data_in={"uvw_folder":uvw_folder,"uvw_file":uvw_file,"Hperc":Hperc,
                                     "index":index,"dx":dx,"dy":dy,"dz":dz,"L_x":L_x,"L_y":L_y,
                                     "L_z":L_z,"rey":rey,"utau":utau,"padding":padding,
                                     "data_folder":data_folder,"umean_file":umean_file,
                                     "urms_file":urms_file,"sym_quad":True,"filvol":filvol,
                                     "shap_folder":shap_folder,"shap_file":shap_file,
                                     "folder":SHAPq_intensity_folder,"file":SHAPq_intensity_file,
                                     "padding":padding,"data_type":data_type,"nsamples":nsamples,
                                     "SHAPrms_file":SHAPrms_file,"SHAPmean_file":SHAPmean_file})
shap_struc.read_struc()

# -----------------------------------------------------------------------------------------------------------------------
# Create the data of the uv structure
# -----------------------------------------------------------------------------------------------------------------------
uv_struc = uv_structure(data_in={"uvw_folder":uvw_folder,"uvw_file":uvw_file,"Hperc":Hperc,"index":index,"dx":dx,
                                 "dy":dy,"dz":dz,"L_x":L_x,"L_y":L_y,"L_z":L_z,"rey":rey,"utau":utau,
                                 "padding":padding,"data_folder":data_folder,"umean_file":umean_file,
                                 "urms_file":urms_file,"sym_quad":True,"filvol":filvol,"shap_folder":shap_folder,
                                 "shap_file":shap_file,"folder":uv_folder,"file":uv_file,"padding":padding,
                                 "data_type":data_type})
uv_struc.read_struc()


# -----------------------------------------------------------------------------------------------------------------------
# calculate the coincidence between the uv and the shap structures as a function of y
# -----------------------------------------------------------------------------------------------------------------------
mat_comb = calc_coinc_mat(data_in={"mat_struc1":uv_struc.mat_struc,
                                   "mat_struc2":shap_struc.mat_struc})["mat_comb"]

    
# -----------------------------------------------------------------------------------------------------------------------
# Read the channel characteristics
# -----------------------------------------------------------------------------------------------------------------------
Data_flow={"folder":uvw_folder,"file":uvw_file,"down_x":dx,"down_y":dy,
           "down_z":dz,"L_x":L_x,"L_y":L_y,"L_z":L_z,"rey":rey,"utau":utau,"umax_file":umax_file}
flowfield = flow_field(data_in=Data_flow)
flowfield.shape_tensor()
flowfield.flow_grid()

# -----------------------------------------------------------------------------------------------------------------------
# Plot the data
# -----------------------------------------------------------------------------------------------------------------------
plot_format_data = {"plot_folder":plot_folder,"plot_coin_folder":plot_coin_folder,"plot_coin_file":plot_coin_file,
                    "xlabel":xlabel,"ylabel":ylabel,"fontsize":fontsize,"figsize_x":figsize_x,"figsize_y":figsize_y,
                    "colormap":colormap,"colornum":colornum,"dpi":dpi,"struc1_lab":struc1_lab,
                    "struc2_lab":struc2_lab,"flowfield":flowfield,"mat_comb":mat_comb,"index_ii":index}
plot_coinc_mat(data_in=plot_format_data)