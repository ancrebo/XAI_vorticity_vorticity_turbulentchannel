# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------------------------------------------------------------------
plot_uv_hist_shapvol.py
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
xlabel      = "$V^+$"
ylabel      = "$\phi$"
zlabel      = "$ $"
fontsize    = 18
figsize_x   = 8
figsize_y   = 7
colormap    = "viridis"
colornum    = 4
fig_name    = "uvstreak_hist_shapvol"
dpi         = 200
padtext_x   = 50
padtext_y   = 10
padtext_z   = 7
linewidth   = 2
lev_min     = 1e-3
bins        = 100
nlev        = 1
labels_pdf  = ["Ejections","Sweeps","Low-velocity\nstreaks","High-velocity\nstreaks"]
colors_pdf  = ['#440154','#0072B2','#F0E442','#E0115F'] #,'#009E73'
# -----------------------------------------------------------------------------------------------------------------------
# Define the names of the files containing the definitios of the parameters
# - folder_def : folder containing the files with the definitions required in the problem
# - chd_str    : file containing the data of the channel
# - folders    : file containing the folder and file structures
# - st_data    : file containing the data of the statistics
# -----------------------------------------------------------------------------------------------------------------------
folder_def  = "P125_83pi_240603_v0_definitions"
chd_str     = "channel_data"
folders_str = "folders_msi"
st_data_str = "stats_data"
sh_data_str = "shap_data"
tr_data_str = "training_data"

# -----------------------------------------------------------------------------------------------------------------------
# Import packages
# -----------------------------------------------------------------------------------------------------------------------
from py_bin.py_class.uv_structure import uv_structure
from py_bin.py_class.streak_structure import streak_structure
import os
from py_bin.py_plots.plot_hist_SHAP_magnitude import plot_hist_SHAP_vol_type
import numpy as np
from py_bin.py_class.shap_config import shap_config

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
field_ini         = 20001
field_fin         = 20002
field_delta       = 1
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
padding           = chd.padding
data_type         = tr_data.data_type
plot_folder       = folders.plot_folder
shapfolder_uv     = folders.shapseg_uv_folder
shapfile_uv       = folders.shapseg_uv_file
shapfolder_streak = folders.shapseg_streak_folder
shapfile_streak   = folders.shapseg_streak_file
model_folder      = folders.model_folder
model_read        = folders.model_read
nfil              = tr_data.nfil
stride            = tr_data.stride
activation        = tr_data.activation
kernel            = tr_data.kernel
pooling           = tr_data.pooling

# -----------------------------------------------------------------------------------------------------------------------
# Create the data of the structure
# -----------------------------------------------------------------------------------------------------------------------
for index in range(field_ini,field_fin,field_delta):

    # -------------------------------------------------------------------------------------------------------------------
    # Read the SHAP for the uv
    # -------------------------------------------------------------------------------------------------------------------
    print("Q events field "+str(index),flush=True)
    uv_struc = uv_structure(data_in={"uvw_folder":uvw_folder,"uvw_file":uvw_file,"Hperc":Hperc,"index":index,"dx":dx,
                                      "dy":dy,"dz":dz,"L_x":L_x,"L_y":L_y,"L_z":L_z,"rey":rey,"utau":utau,
                                      "padding":padding,"data_folder":data_folder,"umean_file":umean_file,
                                      "urms_file":urms_file,"sym_quad":True,"filvol":filvol,"shap_folder":shap_folder,
                                      "shap_file":shap_file,"folder":uv_folder,"file":uv_file,"padding":padding,
                                      "data_type":data_type})
    uv_struc.read_struc()
    uv_struc.detect_quadrant()
    data_shap_uv = {"shap_folder":shapfolder_uv,"shap_file":shapfile_uv,"uvw_folder":uvw_folder,"uvw_file":uvw_file,
                    "padding":padding,"dx":dx,"dy":dy,"dz":dz,"data_folder":data_folder,"umean_file":umean_file,
                    "unorm_file":"-","L_x":L_x,"L_z":L_z,"L_y":L_y,"rey":rey,"utau":utau,"ngpu":None,
                    "field_ini":0,"field_fin":0,"field_delta":1,"model_folder":model_folder,"model_read":model_read,
                    "nfil":nfil,"stride":stride,"activation":activation,"kernel":kernel,"pooling":pooling,"delta_pred":1,
                    "nsamples":0,"nsamples_max":0,"data_type":"float32",
                    "error_file":"-","umax_file":"-","urmspred_file":"-","mean_norm":False,"tfrecord_folder":'-',
                    "nrep_field":None,"shap_batch":0,"repeat_exist":False,"flag_model":True,"read_model":True}
    shconf_uv    = shap_config(data_in = data_shap_uv)
    stru_shap_uv = shconf_uv.read_shap_kernel(data_in = {"index":index})
    shap_uv      = abs(stru_shap_uv["SHAP"])
    shap_ind_uv  = np.array(stru_shap_uv["index_filtered"],dtype="int")
    
    # -------------------------------------------------------------------------------------------------------------------
    # Read the SHAP for the streak
    # -------------------------------------------------------------------------------------------------------------------
    print("Streaks field "+str(index),flush=True)
    streak_struc = streak_structure(data_in={"uvw_folder":uvw_folder,"uvw_file":uvw_file,"Hperc":Hperc,
                                              "index":index,"dx":dx,"dy":dy,"dz":dz,"L_x":L_x,"L_y":L_y,"L_z":L_z,
                                              "rey":rey,"utau":utau,"padding":padding,"data_folder":data_folder,
                                              "umean_file":umean_file,"urms_file":urms_file,"sym_quad":True,
                                              "filvol":filvol,"shap_folder":shap_folder,"shap_file":shap_file,
                                              "folder":streak_folder,"file":streak_file,"padding":padding,
                                              "data_type":data_type})
    streak_struc.read_struc()
    streak_struc.divide_high_low()
    data_shap_streak = {"shap_folder":shapfolder_streak,"shap_file":shapfile_streak,"uvw_folder":uvw_folder,
                        "uvw_file":uvw_file,"padding":padding,"dx":dx,"dy":dy,"dz":dz,"data_folder":data_folder,
                        "umean_file":umean_file,"unorm_file":"-","L_x":L_x,"L_z":L_z,"L_y":L_y,"rey":rey,"utau":utau,
                        "ngpu":None,"field_ini":0,"field_fin":0,"field_delta":1,"model_folder":model_folder,
                        "model_read":model_read,"nfil":nfil,"stride":stride,"activation":activation,"kernel":kernel,
                        "pooling":pooling,"delta_pred":1,"nsamples":0,"nsamples_max":0,"data_type":"float32",
                        "error_file":"-","umax_file":"-","urmspred_file":"-","mean_norm":False,"tfrecord_folder":'-',
                        "nrep_field":None,"shap_batch":0,"repeat_exist":False,"flag_model":True,"read_model":True}
    shconf_streak    = shap_config(data_in = data_shap_streak)
    stru_shap_streak = shconf_streak.read_shap_kernel(data_in = {"index":index})
    shap_streak      = abs(stru_shap_streak["SHAP"])
    shap_ind_streak  = np.array(stru_shap_streak["index_filtered"],dtype="int")
    
        
        
    # -------------------------------------------------------------------------------------------------------------------
    # Create the mesh for the pdf
    # -------------------------------------------------------------------------------------------------------------------
    print("Histograms  field "+str(index),flush=True)
    vol_uv      = uv_struc.structures.vol[shap_ind_uv]
    event_uv    = uv_struc.event[shap_ind_uv]
    shapQ2      = shap_uv[event_uv==2]
    volQ2       = vol_uv[event_uv==2]
    shapQ4      = shap_uv[event_uv==4]
    volQ4       = vol_uv[event_uv==4]
    shap_max_uv = np.max(shap_uv)
    shap_min_uv = np.min(shap_uv)
    vol_max_uv  = np.max(vol_uv) 
    vol_min_uv  = np.min(vol_uv)
    
    
    vol_streak      = streak_struc.structures.vol[shap_ind_streak]
    type_streak     = streak_struc.type[shap_ind_streak]
    shap_streakh    = shap_streak[type_streak==1]
    vol_streakh     = vol_streak[type_streak==1]
    shap_streakl    = shap_streak[type_streak==-1]
    vol_streakl     = vol_streak[type_streak==-1]
    shap_max_streak = np.max(shap_streak)
    shap_min_streak = np.min(shap_streak)
    vol_max_streak  = np.max(vol_streak) 
    vol_min_streak  = np.min(vol_streak)
    
    shap_min = np.min([shap_min_uv,shap_min_streak])
    shap_max = np.max([shap_max_uv,shap_max_streak])
    vol_min  = np.min([vol_min_uv,vol_min_streak])
    vol_max  = np.max([vol_max_uv,vol_max_streak])
    
    
    bins_shap = np.linspace(shap_min,shap_max,bins)
    bins_vol  = np.linspace(vol_min,vol_max,bins)
    
    # ----------------------------------------------------------------------------------------------------------------
    # Calculate the pdf for shap and volume
    # ----------------------------------------------------------------------------------------------------------------
    hist_shapvolQ2,hist_volQ2,hist_shapQ2 = np.histogram2d(volQ2,shapQ2,bins=(bins_vol,bins_shap))
    hist_volQ2                            = hist_volQ2[:-1]+np.diff(hist_volQ2)/2
    hist_shapQ2                           = hist_shapQ2[:-1]+np.diff(hist_shapQ2)/2
    grid_volQ2,grid_shapQ2                = np.meshgrid(hist_volQ2,hist_shapQ2)
    grid_shapvolQ2                        = hist_shapvolQ2.T.copy()
    if index == field_ini:
        grid_shapvol_totQ2  = grid_shapvolQ2
    else:
        grid_shapvol_totQ2 += grid_shapvolQ2
    hist_shapvolQ4,hist_volQ4,hist_shapQ4 = np.histogram2d(volQ4,shapQ4,bins=(bins_vol,bins_shap))
    hist_volQ4                            = hist_volQ4[:-1]+np.diff(hist_volQ4)/2
    hist_shapQ4                           = hist_shapQ4[:-1]+np.diff(hist_shapQ4)/2
    grid_volQ4,grid_shapQ4                = np.meshgrid(hist_volQ4,hist_shapQ4)
    grid_shapvolQ4                        = hist_shapvolQ4.T.copy()
    if index == field_ini:
        grid_shapvol_totQ4  = grid_shapvolQ4
    else:
        grid_shapvol_totQ4 += grid_shapvolQ4
    hist_shapvol_streakh,hist_vol_streakh,hist_shap_streakh = np.histogram2d(vol_streakh,shap_streakh,
                                                                             bins=(bins_vol,bins_shap))
    hist_vol_streakh                                        = hist_vol_streakh[:-1]+np.diff(hist_vol_streakh)/2
    hist_shap_streakh                                       = hist_shap_streakh[:-1]+np.diff(hist_shap_streakh)/2
    grid_vol_streakh,grid_shap_streakh                      = np.meshgrid(hist_vol_streakh,hist_shap_streakh)
    grid_shapvol_streakh                                    = hist_shapvol_streakh.T.copy()
    if index == field_ini:
        grid_shapvol_tot_streakh  = grid_shapvol_streakh
    else:
        grid_shapvol_tot_streakh += grid_shapvol_streakh
    hist_shapvol_streakl,hist_vol_streakl,hist_shap_streakl = np.histogram2d(vol_streakl,shap_streakl,
                                                                             bins=(bins_vol,bins_shap))
    hist_vol_streakl                                        = hist_vol_streakl[:-1]+np.diff(hist_vol_streakl)/2
    hist_shap_streakl                                       = hist_shap_streakl[:-1]+np.diff(hist_shap_streakl)/2
    grid_vol_streakl,grid_shap_streakl                      = np.meshgrid(hist_vol_streakl,hist_shap_streakl)
    grid_shapvol_streakl                                    = hist_shapvol_streakl.T.copy()
    if index == field_ini:
        grid_shapvol_tot_streakl  = grid_shapvol_streakl
    else:
        grid_shapvol_tot_streakl += grid_shapvol_streakl

grid_shapvol_totQ2                                   /= np.max(grid_shapvol_totQ2)
shapcontentQ2                                         = grid_shapvolQ2[np.where(grid_shapvol_totQ2>=lev_min)]
grid_shapvol_totQ2[grid_shapvol_totQ2==0]             = 1e-20
grid_shapvol_totQ4                                   /= np.max(grid_shapvol_totQ4)
shapcontentQ4                                         = grid_shapvolQ4[np.where(grid_shapvol_totQ4>=lev_min)]
grid_shapvol_totQ4[grid_shapvol_totQ4==0]             = 1e-20
grid_shapvol_tot_streakh                             /= np.max(grid_shapvol_tot_streakh)
shapcontent_streakh                                   = grid_shapvol_streakh[np.where(grid_shapvol_tot_streakh>=lev_min)]
grid_shapvol_tot_streakh[grid_shapvol_tot_streakh==0] = 1e-20
grid_shapvol_tot_streakl                             /= np.max(grid_shapvol_tot_streakl)
shapcontent_streakl                                   = grid_shapvol_streakl[np.where(grid_shapvol_tot_streakl>=lev_min)]
grid_shapvol_tot_streakl[grid_shapvol_tot_streakl==0] = 1e-20
        

# -----------------------------------------------------------------------------------------------------------------------
# Plot the SHAP uv
# -----------------------------------------------------------------------------------------------------------------------

plot_hist_SHAP_vol_type(data_in={"grid_shapvol":[grid_shapvol_totQ2,grid_shapvol_totQ4,grid_shapvol_tot_streakl,
                                                 grid_shapvol_tot_streakh],
                                 "grid_vol":[grid_volQ2,grid_volQ4,grid_vol_streakl,grid_vol_streakh],
                                 "grid_shap":[grid_shapQ2,grid_shapQ4,grid_shap_streakl,grid_shap_streakl],
                                 "plot_folder":plot_folder,"xlabel":xlabel,"ylabel":ylabel,
                                 "zlabel":zlabel,"fontsize":fontsize,"figsize_x":figsize_x,"figsize_y":figsize_y,
                                 "colormap":colormap,"colornum":colornum,"fig_name":fig_name,"dpi":dpi,"cmap_flag":True,
                                 "shap_max":shap_max,"shap_min":shap_min,"vol_max":vol_max,"vol_min":vol_min,
                                 "padtext":[padtext_x,padtext_y,padtext_z],
                                 "lev_min":lev_min,"nlev":nlev,"linewidth":linewidth,"labels_pdf":labels_pdf,
                                 "colors_pdf":colors_pdf})