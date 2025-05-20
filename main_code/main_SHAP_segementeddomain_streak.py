# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------------------------------------------------------------------
main_SHAP.py
-------------------------------------------------------------------------------------------------------------------------
Created on Mon Mar 18 10:12:32 2024

@author: Andres Cremades Botella

File to calculate the SHAP values of each grid-point. The SHAP values are calculated with the Expected Gradients
methdology. The link to the SHAP package and the reference paper are provided below:
    - GradientExplainer: https://shap.readthedocs.io/en/latest/generated/shap.GradientExplainer.html
    - Article: Erion, G., Janizek, J. D., Sturmfels, P., Lundberg, S. M., & Lee, S. I. (2021).
               Improving performance of deep learning models with axiomatic attribution priors and expected gradients.
               Nature machine intelligence, 3(7), 620-631.
In order to launch the code the following variables need to be fulfilled:
    - folder_def  : (str) name of the folder containing the files for configuring the case of analysis.
    - chd_str     : (str) name of the file containing the data of the channel.
    - folders_str : (str) name of the file containing the folders and files used in the problem.
    - tr_data_str : (str) name of the file containing the information required for the training.
    - sh_data_str : (str) name of the file containing the information required for calculating the SHAP values
"""
# ----------------------------------------------------------------------------------------------------------------------
# Define the names of the files containing the definitios of the parameters
#     - folder_def  : folder containing the files with the definitions required in the problem
#     - chd_str     : file containing the data of the channel
#     - folders_str : file containing the folder and file structures
#     - tr_data_str : file containing the data of the training
#     - sh_data_str : file containing the data for the SHAP values 
#     - struct_type : selection of the structure type (Qevent|streak|chong)
# ----------------------------------------------------------------------------------------------------------------------
folder_def  = "P125_83pi_240603_v0_definitions" #"d20240703_definitions"#
chd_str     = "channel_data"
folders_str = "folders" #"folders_local" #
tr_data_str = "training_data"
sh_data_str = "shap_data"
struct_type = "streak" #"Qevent"

# ----------------------------------------------------------------------------------------------------------------------
# Load the packages
# ----------------------------------------------------------------------------------------------------------------------
import py_bin.py_class.shap_config as sc
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
exec("from "+folder_def+" import "+tr_data_str+" as tr_data")
exec("from "+folder_def+" import "+sh_data_str+" as sh_data")

# ----------------------------------------------------------------------------------------------------------------------
# Load the channel data to import the information regarding the channel size and the friction Reynolds number 
# and velocity
#     - L_x     : Channel size in the streamwise direction
#     - L_z     : Channel size in the spanwise direction
#     - L_y     : Half of the channel width
#     - rey     : Friction Reynolds number
#     - utau    : Friction velocity
#     - dx      : Downsampling in x
#     - dy      : Downsampling in y
#     - dz      : Downsampling in z
#     - padding : Number of nodes of the padding
# ----------------------------------------------------------------------------------------------------------------------
L_x     = chd.L_x
L_z     = chd.L_z
L_y     = chd.L_y
rey     = chd.rey
utau    = chd.utau
dx      = chd.dx
dy      = chd.dy
dz      = chd.dz
padding = chd.padding

# ----------------------------------------------------------------------------------------------------------------------
# Define the data of the model definition: data, padding, downsampling...
#     - uvw_folder      : Folder of the velocity data
#     - uvw_file        : This file does not contain the file index
#     - data_folder     : Folder for storing the data of the model
#     - umean_file      : File for the mean velocity
#     - unorm_file      : File for the normalization of the velocity
#     - uvw_folder_tf   : Folder of the velocity data with tensorflow format
#     - uvw_folderii_tf : File of the velocity data with tensorflow format
# ----------------------------------------------------------------------------------------------------------------------
uvw_folder      = folders.uvw_folder
uvw_file        = folders.uvw_file
data_folder     = folders.data_folder
umean_file      = folders.umean_file
unorm_file      = folders.unorm_file
uvw_folder_tf   = folders.uvw_folder_tf
uvw_folderii_tf = folders.uvw_folderii_tf

# ----------------------------------------------------------------------------------------------------------------------
# Define the data for the training.
#     - ngpu            : Number of gpus
#     - field_ini       : Initial field of the shaps
#     - field_fin       : Final field of the shaps
#     - field_delta     : Separation between the files
#     - read_model      : Flag to define or read the model (False=define, True=read)
#     - model_folder    : Folder of the trained model files
#     - model_name      : Name of the trained model file
#     - nfil            : Number of filters of the first layer of the Unet
#     - stride          : Stride of the Unet
#     - activation      : Activation function
#     - kernel          : Kernel size of the unet
#     - pooling         : Size of the poolings of the unet
#     - delta_pred      : Number of fields to advance the prediction
#     - data_type       : Format of the data of the training
#     - error_file      : file to store the error
#     - umax_file       : file to store the maximum and minimum velocity
#     - urmspred_file   : file to store the rms predicted by the model 
#     - tfrecord_folder : folder of the tfrecods
#     - nrep_field      : number of repetitions of each field for calculating the SHAP values
#     - shap_batch      : batch size used for the shap
#     - repeat exist    : flag for repeating an existing file (True: recalculate, False: skip)
# ----------------------------------------------------------------------------------------------------------------------
ngpu            = tr_data.ngpu
field_ini       = sh_data.field_ini
field_fin       = sh_data.field_fin
field_delta     = sh_data.field_delta
read_model      = tr_data.read_model
model_folder    = folders.model_folder
model_read      = folders.model_read
nfil            = tr_data.nfil
stride          = tr_data.stride
activation      = tr_data.activation
kernel          = tr_data.kernel
pooling         = tr_data.pooling
delta_pred      = tr_data.delta_pred
nsamples        = sh_data.nsamples
nsamples_max    = sh_data.nsamples_max
data_type       = tr_data.data_type
error_file      = folders.error_file
umax_file       = folders.umax_file
urmspred_file   = folders.urmspred_file
mean_norm       = tr_data.mean_norm
tfrecord_folder = folders.tfrecord_folder
nrep_field      = sh_data.nrep_field
shap_batch      = sh_data.shap_batch
repeat_exist    = sh_data.repeat_exist
filvol          = chd.filvol
structurefolder = {"Qevent_folder":folders.uv_folder,"Qevent_file":folders.uv_file,"Qevent_fun":"uv_structure",
                   "Qevent_saveshap_folder":folders.shapseg_uv_folder,"Qevent_saveshap_file":folders.shapseg_uv_file,
                   "streak_folder":folders.streak_folder,"streak_file":folders.streak_file,
                   "streak_fun":"streak_structure","streak_saveshap_folder":folders.shapseg_streak_folder,
                   "streak_saveshap_file":folders.shapseg_streak_file,"chong_folder":folders.chong_folder,
                   "chong_file":folders.chong_file,"chong_fun":"chong_structure",
                   "chong_saveshap_folder":folders.shapseg_chong_folder,
                   "chong_saveshap_file":folders.shapseg_chong_file}

# ----------------------------------------------------------------------------------------------------------------------
# Define the data of the model definition: data, padding, downsampling...
#     - shap_folder     : Folder of the SHAP values
#     - shap_file       : File of the SHAP values
# ----------------------------------------------------------------------------------------------------------------------
shap_folder     = structurefolder[struct_type+"_saveshap_folder"]
shap_file       = structurefolder[struct_type+"_saveshap_file"]
# ----------------------------------------------------------------------------------------------------------------------
# Define dict containing the information needed for the shap model
# ----------------------------------------------------------------------------------------------------------------------
data_shap = {"shap_folder":shap_folder,"shap_file":shap_file,"uvw_folder":uvw_folder,"uvw_file":uvw_file,
             "padding":padding,"dx":dx,"dy":dy,"dz":dz,"data_folder":data_folder,"umean_file":umean_file,
             "unorm_file":unorm_file,"L_x":L_x,"L_z":L_z,"L_y":L_y,"rey":rey,"utau":utau,"ngpu":ngpu,
             "field_ini":field_ini,"field_fin":field_fin,"field_delta":field_delta,"model_folder":model_folder,
             "model_read":model_read,"nfil":nfil,"stride":stride,"activation":activation,"kernel":kernel,
             "pooling":pooling,"delta_pred":delta_pred,"nsamples":nsamples,"nsamples_max":nsamples_max,
             "data_type":data_type,"error_file":error_file,"umax_file":umax_file,"urmspred_file":urmspred_file,
             "mean_norm":mean_norm,"tfrecord_folder":tfrecord_folder,"nrep_field":nrep_field,"shap_batch":shap_batch,
             "repeat_exist":repeat_exist,"flag_model":True}
shap_model = sc.shap_config(data_in=data_shap)
shap_model.calc_kernelSHAP(data_in={"structurefolder":structurefolder,"struct_type":struct_type,"filvol":filvol})

