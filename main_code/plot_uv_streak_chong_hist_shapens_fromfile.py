# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------------------------------------------------------------------
plot_uv_hist_shapens.py
-------------------------------------------------------------------------------------------------------------------------
Created on Tue Jun 18 11:30:49 2024

@author: Andres Cremades Botella

Function to plot the 3d Reynolds stress structures:
    - folder_def  : (str) name of the folder containing the files for configuring the case of analysis.
    - chd_str     : (str) name of the file containing the data of the channel.
    - folders_str : (str) name of the file containing the folders and files used in the problem.
    - st_data_str : (str) name of the file containing the information required for the statistics.
For more information about the tangential Reynolds stress structures:
    - Lozano-Durán, A., Flores, O., & Jiménez, J. (2012). The three-dimensional structure of momentum transfer in
      turbulent channels. Journal of Fluid Mechanics, 694, 100-130.
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
xlabel            = "$\Omega/\Omega_{tot}\cdot 10^{-2}$"
ylabel            = r"$|\phi^{\mathbf{u}}|\cdot 10^{-1}$"
xlabelvol         = r"$\Omega/(\Omega_{tot}V^+)\cdot 10^{-8}$"
ylabelvol         = r"$|\phi^{\mathbf{u}}/V^+|\cdot 10^{-7}$"
zlabel            = "$ $"
fontsize          = 24
figsize_x         = 9
figsize_y         = 7
colormap          = "viridis"
colornum          = 4
fig_name          = "uvstreakchong_hist_shap_ens_b_0001"
fig_namevol       = "uvstreakchong_hist_shapvol_ens_b_0001"
file_frequency    = "uvstreakchong_hist_shap_ens_b_0001.h5"
file_frequencyvol = "uvstreakchong_hist_shapvol_ens_b_0001.h5"
dpi               = 200
padtext_x         = 50
padtext_y         = 10
padtext_z         = 7
linewidth         = 2
lev_min           = 1e-4
bins              = 200
nlev              = 1
labels_pdf        = ["Ejections","Sweeps","Low-velocity\nstreaks","High-velocity\nstreaks","Vortices"]
colors_pdf        = ['#440154','#0072B2','#009E73','#E0115F','#F0E442'] #
foldersave2       = "/data2/andres/P125_83pi/sta"
# -----------------------------------------------------------------------------------------------------------------------
# Define the names of the files containing the definitios of the parameters
# - folder_def : folder containing the files with the definitions required in the problem
# - chd_str    : file containing the data of the channel
# - folders    : file containing the folder and file structures
# - st_data    : file containing the data of the statistics
# -----------------------------------------------------------------------------------------------------------------------
folder_def  = "P125_83pi_240603_v0_definitions"
chd_str     = "channel_data"
folders_str = "folders_msi" #"folders_local"
st_data_str = "stats_data"
sh_data_str = "shap_data"
tr_data_str = "training_data"

# -----------------------------------------------------------------------------------------------------------------------
# Import packages
# -----------------------------------------------------------------------------------------------------------------------
from py_bin.py_class.uv_structure import uv_structure
from py_bin.py_class.streak_structure import streak_structure
from py_bin.py_class.chong_structure import chong_structure
import os
from py_bin.py_plots.plot_hist_SHAP_magnitude import plot_hist_SHAP_vol_type
import numpy as np
from py_bin.py_class.shap_config import shap_config
import h5py

    

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
#     - index       : index of the field
#     - Hperc       : percolation index
#     - uvw_folder  : folder of the flow field data
#     - uvw_file    : file of the flow field data
#     - umean_file  : file to save the mean velocity
#     - data_folder : folder to store the calculated data
#     - dx          : downsampling in x
#     - dy          : downsampling in y
#     - dz          : downsampling in z
#     - L_x         : length of the channel in the streamwise direction
#     - L_y         : half-width of the channel in the wall-normal direction
#     - L_z         : length of the channel in the spanwise direction
#     - urms_file   : file to save the rms of the velocity
#     - rey         : Friction Reynolds number
#     - utau        : Friction velocity
#     - padding     : padding of the flow field
#     - sym_quad    : flag for using the symmetry in the direction 2 of the field for the quadrant selection
#     - filvol      : volume for filtering the structures+
#     - shap_folder : folder of the shap values
#     - shap_folder : file of the shap values
#     - uv_folder   : folder of the uv structures
#     - uv_file     : file of the uv structures
#     - padding     : padding of the field
#     - data_type   : type of data used by the model
# -----------------------------------------------------------------------------------------------------------------------
field_ini         = sh_data.field_ini
field_fin         = sh_data.field_fin
field_delta       = sh_data.field_delta
Hperc             = 1.41
uvw_folder        = folders.uvw_folder
uvw_file          = folders.uvw_file
umean_file        = folders.umean_file
data_folder       = folders.data_folder
dx                = chd.dx
dy                = chd.dy
dz                = chd.dz
L_x               = chd.L_x
L_y               = chd.L_y
L_z               = chd.L_z
urms_file         = folders.urms_file
rey               = chd.rey
utau              = chd.utau
padding           = chd.padding
sym_quad          = True
filvol            = chd.filvol
shap_folder       = folders.shap_folder
shap_file         = folders.shap_file
uv_folder         = folders.uv_folder
uv_file           = folders.uv_file
streak_folder     = folders.streak_folder
streak_file       = folders.streak_file
chong_folder      = folders.chong_folder
chong_file        = folders.chong_file
padding           = chd.padding
data_type         = tr_data.data_type
plot_folder       = folders.plot_folder
shapfolder_uv     = folders.shapseg_uv_folder
shapfile_uv       = folders.shapseg_uv_file
shapfolder_streak = folders.shapseg_streak_folder
shapfile_streak   = folders.shapseg_streak_file
shapfolder_chong  = folders.shapseg_chong_folder
shapfile_chong    = folders.shapseg_chong_file
model_folder      = folders.model_folder
model_read        = folders.model_read
nfil              = tr_data.nfil
stride            = tr_data.stride
activation        = tr_data.activation
kernel            = tr_data.kernel
pooling           = tr_data.pooling


# -----------------------------------------------------------------------------------------------------------------------
# Read SHAP vs ens
# -----------------------------------------------------------------------------------------------------------------------
ff_freq = h5py.File(data_folder+"/"+file_frequency,"r")
grid_shapQ2              = np.array(ff_freq['grid_shapQ2'])*1e1
grid_ensQ2               = np.array(ff_freq['grid_ensQ2'])*1e3
grid_shapens_totQ2       = np.array(ff_freq['grid_shapens_totQ2'])
grid_shapQ4              = np.array(ff_freq['grid_shapQ4'])*1e1
grid_ensQ4               = np.array(ff_freq['grid_ensQ4'])*1e3
grid_shapens_totQ4       = np.array(ff_freq['grid_shapens_totQ4'])
grid_shap_streakh        = np.array(ff_freq['grid_shap_streakh'])*1e1
grid_ens_streakh         = np.array(ff_freq['grid_ens_streakh'])*1e3
grid_shapens_tot_streakh = np.array(ff_freq['grid_shapens_tot_streakh'])
grid_shap_streakl        = np.array(ff_freq['grid_shap_streakl'])*1e1
grid_ens_streakl         = np.array(ff_freq['grid_ens_streakl'])*1e3
grid_shapens_tot_streakl = np.array(ff_freq['grid_shapens_tot_streakl'])
grid_shap_chong          = np.array(ff_freq['grid_shap_chong'])*1e1
grid_ens_chong           = np.array(ff_freq['grid_ens_chong'])*1e3
grid_shapens_tot_chong   = np.array(ff_freq['grid_shapens_tot_chong'])
last_field               = np.array(ff_freq['last_field'])
n_fields                 = np.array(ff_freq['n_fields'])
ens_min                  = np.array(ff_freq['ens_min'])*1e3
ens_max                  = np.array(ff_freq['ens_max'])*1e3
shap_min                 = np.array(ff_freq['shap_min'])*1e1
shap_max                 = np.array(ff_freq['shap_max'])*1e1
ff_freq.close()


grid_shapens_totQ2                                   /= np.max(grid_shapens_totQ2)
grid_shapens_totQ2[grid_shapens_totQ2==0]             = 1e-20
grid_shapens_totQ4                                   /= np.max(grid_shapens_totQ4)
grid_shapens_totQ4[grid_shapens_totQ4==0]             = 1e-20
grid_shapens_tot_streakh                             /= np.max(grid_shapens_tot_streakh)
grid_shapens_tot_streakh[grid_shapens_tot_streakh==0] = 1e-20
grid_shapens_tot_streakl                             /= np.max(grid_shapens_tot_streakl)
grid_shapens_tot_streakl[grid_shapens_tot_streakl==0] = 1e-20
grid_shapens_tot_chong                               /= np.max(grid_shapens_tot_chong)
grid_shapens_tot_chong[grid_shapens_tot_chong==0]     = 1e-20
        

shap_min = 0
shap_max = 0.175*1e1
ens_min  = 0
ens_max  = 6.5e-3*1e3
# -----------------------------------------------------------------------------------------------------------------------
# Plot the SHAP uv
# -----------------------------------------------------------------------------------------------------------------------

plot_hist_SHAP_vol_type(data_in={"grid_shapvol":[grid_shapens_totQ2,grid_shapens_totQ4,grid_shapens_tot_streakl,
                                                 grid_shapens_tot_streakh,grid_shapens_tot_chong],
                                 "grid_vol":[grid_ensQ2,grid_ensQ4,grid_ens_streakl,grid_ens_streakh,grid_ens_chong],
                                 "grid_shap":[grid_shapQ2,grid_shapQ4,grid_shap_streakl,
                                              grid_shap_streakl,grid_shap_chong],
                                 "plot_folder":plot_folder,"xlabel":xlabel,"ylabel":ylabel,
                                 "zlabel":zlabel,"fontsize":fontsize,"figsize_x":figsize_x,"figsize_y":figsize_y,
                                 "colormap":colormap,"colornum":colornum,"fig_name":fig_name,"dpi":dpi,"cmap_flag":True,
                                 "shap_max":shap_max,"shap_min":shap_min,"vol_max":ens_max,"vol_min":ens_min,
                                 "padtext":[padtext_x,padtext_y,padtext_z],
                                 "lev_min":lev_min,"nlev":nlev,"linewidth":linewidth,"labels_pdf":labels_pdf,
                                 "colors_pdf":colors_pdf})


# -----------------------------------------------------------------------------------------------------------------------
# Read SHAP vol vs ens vol
# -----------------------------------------------------------------------------------------------------------------------
ff_freq = h5py.File(data_folder+"/"+file_frequencyvol,"r")
grid_shapvolQ2                  = np.array(ff_freq['grid_shapQ2'])*1e7
grid_ensvolQ2                   = np.array(ff_freq['grid_ensQ2'])*1e8
grid_shapvol_ensvol_totQ2       = np.array(ff_freq['grid_shapens_totQ2'])
grid_shapvolQ4                  = np.array(ff_freq['grid_shapQ4'])*1e7
grid_ensvolQ4                   = np.array(ff_freq['grid_ensQ4'])*1e8
grid_shapvol_ensvol_totQ4       = np.array(ff_freq['grid_shapens_totQ4'])
grid_shapvol_streakh            = np.array(ff_freq['grid_shap_streakh'])*1e7
grid_ensvol_streakh             = np.array(ff_freq['grid_ens_streakh'])*1e8
grid_shapvol_ensvol_tot_streakh = np.array(ff_freq['grid_shapens_tot_streakh'])
grid_shapvol_streakl            = np.array(ff_freq['grid_shap_streakl'])*1e7
grid_ensvol_streakl             = np.array(ff_freq['grid_ens_streakl'])*1e8
grid_shapvol_ensvol_tot_streakl = np.array(ff_freq['grid_shapens_tot_streakl'])
grid_shapvol_chong              = np.array(ff_freq['grid_shap_chong'])*1e7
grid_ensvol_chong               = np.array(ff_freq['grid_ens_chong'])*1e8
grid_shapvol_ensvol_tot_chong   = np.array(ff_freq['grid_shapens_tot_chong'])
last_field                      = np.array(ff_freq['last_field'])
n_fields                        = np.array(ff_freq['n_fields'])
ensvol_min                      = np.array(ff_freq['ens_min'])*1e8
ensvol_max                      = np.array(ff_freq['ens_max'])*1e8
shapvol_min                     = np.array(ff_freq['shap_min'])*1e7
shapvol_max                     = np.array(ff_freq['shap_max'])*1e7
ff_freq.close()

grid_shapvol_ensvol_totQ2                                          /= np.max(grid_shapvol_ensvol_totQ2)
grid_shapvol_ensvol_totQ2[grid_shapvol_ensvol_totQ2==0]             = 1e-20
grid_shapvol_ensvol_totQ4                                          /= np.max(grid_shapvol_ensvol_totQ4)
grid_shapvol_ensvol_totQ4[grid_shapvol_ensvol_totQ4==0]             = 1e-20
grid_shapvol_ensvol_tot_streakh                                    /= np.max(grid_shapvol_ensvol_tot_streakh)
grid_shapvol_ensvol_tot_streakh[grid_shapvol_ensvol_tot_streakh==0] = 1e-20
grid_shapvol_ensvol_tot_streakl                                    /= np.max(grid_shapvol_ensvol_tot_streakl)
grid_shapvol_ensvol_tot_streakl[grid_shapvol_ensvol_tot_streakl==0] = 1e-20
grid_shapvol_ensvol_tot_chong                                      /= np.max(grid_shapvol_ensvol_tot_chong)
grid_shapvol_ensvol_tot_chong[grid_shapvol_ensvol_tot_chong==0]     = 1e-20

shapvol_min = 0
shapvol_max = 6e-7*1e7
ensvol_min  = 0
ensvol_max  = 2.2e-8*1e8

plot_hist_SHAP_vol_type(data_in={"grid_shapvol":[grid_shapvol_ensvol_totQ2,grid_shapvol_ensvol_totQ4,
                                                  grid_shapvol_ensvol_tot_streakl,grid_shapvol_ensvol_tot_streakh,
                                                  grid_shapvol_ensvol_tot_chong],
                                  "grid_vol":[grid_ensvolQ2,grid_ensvolQ4,grid_ensvol_streakl,
                                              grid_ensvol_streakh,grid_ensvol_chong],
                                  "grid_shap":[grid_shapvolQ2,grid_shapvolQ4,grid_shapvol_streakl,
                                              grid_shapvol_streakl,grid_shapvol_chong],
                                  "plot_folder":plot_folder,"xlabel":xlabelvol,"ylabel":ylabelvol,
                                  "zlabel":zlabel,"fontsize":fontsize,"figsize_x":figsize_x,"figsize_y":figsize_y,
                                  "colormap":colormap,"colornum":colornum,"fig_name":fig_namevol,"dpi":dpi,"cmap_flag":True,
                                  "shap_max":shapvol_max,"shap_min":shapvol_min,"vol_max":ensvol_max,"vol_min":ensvol_min,
                                  "padtext":[padtext_x,padtext_y,padtext_z],
                                  "lev_min":lev_min,"nlev":nlev,"linewidth":linewidth,"labels_pdf":labels_pdf,
                                  "colors_pdf":colors_pdf})