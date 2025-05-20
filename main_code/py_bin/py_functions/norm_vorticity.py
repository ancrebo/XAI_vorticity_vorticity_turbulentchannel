# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------------------------------------------------------------------
norm_velocity.py
-------------------------------------------------------------------------------------------------------------------------
Created on Fri Mar 22 11:59:50 2024

@author:  Andres Cremades Botella

File to normalize the velocity fields. The normalization generates values between 0 and 1 using the minimum and the 
maximum of the velocity values. The file contains the following functions:
    Functions:
        - norm_velocity : function for normalize the velocity
        - dim_velocity  : function for dimensionalize the velocity
"""

# -----------------------------------------------------------------------------------------------------------------------
# Import packages for all the functions
# -----------------------------------------------------------------------------------------------------------------------
import sys
import numpy as np

# -----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------
# Define the functions
# -----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------

def norm_vorticity(data_in={"vor_x":[],"vor_y":[],"vor_z":[],"folder_data":"Data","vornorm_file":"norm.txt",
                            "data_type":"float32","mean_norm":False}):
    """
    .....................................................................................................................
    # norm_vorticity: function for normalize the vorticity
    .....................................................................................................................
    Parameters
    ----------
    data_in : dict, optional
        Data for the normalization of the velocity.
        The default is {"vor_x":[],"vor_y":[],"vor_z":[],"folder_data":"Data","vornorm_file":"norm.txt",
                        "data_type":"float32","mean_norm":False}.
        Data:
            - vor_x        : streamwise vorticity
            - vor_y        : wall-normal vorticity
            - vor_z        : spanwise vorticity
            - folder_data  : path to the folder containing the normalization
            - vornorm_file : file with the normalization values
            - data_type    : type of the data (float32,float16...)
            - mean_norm    : choose normalizing with the mean

    Returns
    -------
    data_out : dict
        Normalized velocity. The velocity is configured in float16 to save memory during the training
        Data:
            - vor_x_norm : streamwise normalized vorticity
            - vor_y_norm : wall-normal normalized vorticity
            - vor_z_norm : spanwise normalized vorticity

    """    
    # -------------------------------------------------------------------------------------------------------------------
    # The normalization of the velocity fluctuations is generated in float format
    # -------------------------------------------------------------------------------------------------------------------
    vor_x     = np.array(data_in["vor_x"],dtype="float")               # velocity fluctuation in the streamwise direction.
    vor_y     = np.array(data_in["vor_y"],dtype="float")               # velocity fluctuation in the wall-normal direction.
    vor_z     = np.array(data_in["vor_z"],dtype="float")               # velocity fluctuation in the spanwise direction.
    mean_norm = bool(data_in["mean_norm"])
    
    # -------------------------------------------------------------------------------------------------------------------
    # import packages
    # -------------------------------------------------------------------------------------------------------------------
    if mean_norm:
        from py_bin.py_functions.normalization_normaldist import read_norm
    else:
        from py_bin.py_functions.normalization_vor import read_norm
        
    # -------------------------------------------------------------------------------------------------------------------
    # Check datatype
    # -------------------------------------------------------------------------------------------------------------------
    if "data_type" in data_in.keys():
        data_type = str(data_in["data_type"])                       # definition of the data type.
        if not (data_type=="float32" or data_type=="float16"):
            data_type = "float32"
    else:
        print("[trainvali_data.py:data_traintest_tf] Data type needs to be selected.")
        sys.exit()
    
    # -------------------------------------------------------------------------------------------------------------------
    # Read the normalization parameters
    # -------------------------------------------------------------------------------------------------------------------
    folder_data    = str(data_in["folder_data"])               # folder of the generated data
    vornorm_file   = str(data_in["vornorm_file"])                # file of the normalization data
    vornorm_data   = {"folder":folder_data,"file":vornorm_file}
    try:
        norm_param = read_norm(vornorm_data)
    except:
        print('Normalization file could not be located. Calculation is stopped...',flush=True)
        sys.exit()
    
    if mean_norm:
        print('Mean and std not implemented',flush=True)
        sys.exit()
    else:
        vor_x_max = float(norm_param["vor_x_max"])
        vor_y_max = float(norm_param["vor_y_max"])
        vor_z_max = float(norm_param["vor_z_max"])
        vor_x_min = float(norm_param["vor_x_min"])
        vor_y_min = float(norm_param["vor_y_min"])
        vor_z_min = float(norm_param["vor_z_min"])
        
        # ---------------------------------------------------------------------------------------------------------------
        # Define the normalized fields using float 16 format
        # ---------------------------------------------------------------------------------------------------------------
        vor_x_norm = np.array((vor_x-vor_x_min)/(vor_x_max-vor_x_min),dtype=data_type)
        vor_y_norm = np.array((vor_y-vor_y_min)/(vor_y_max-vor_y_min),dtype=data_type)
        vor_z_norm = np.array((vor_z-vor_z_min)/(vor_z_max-vor_z_min),dtype=data_type)
    data_out = {"vor_x_norm":vor_x_norm,"vor_y_norm":vor_y_norm,"vor_z_norm":vor_z_norm}
    return data_out

def dim_velocity(data_in={"vor_x_norm":[],"vor_y_norm":[],"vor_z_norm":[],"folder_data":"Data","vornorm_file":"norm.txt",
                          "data_type":"float32","mean_norm":False}):
    """
    .....................................................................................................................
    # dim_velocity: function for dimensionalize the velocity
    .....................................................................................................................
    Parameters
    ----------
    data_in : dict, optional
        Data for the non-normalization of the velocity.
        The default is {"vor_x_norm":[],"vor_y_norm":[],"vor_z_norm":[],"folder_data":"Data","vornorm_file":"norm.txt",
                        "data_type":"float32","mean_norm":False}.
        Data:
            - vor_x_norm   : streamwise vorticity
            - vor_x_norm   : wall-normal vorticity
            - vor_x_norm   : spanwise vorticity
            - folder_data  : path to the folder containing the normalization
            - vornorm_file : file with the normalization values
            - data_type    : type of the data (float32,float16...)
            - mean_norm    : choose normalizing with the mean

    Returns
    -------
    data_out : dict
        Normalized velocity. The velocity is configured in float16 to save memory during the training
        Data:
            vor_x : streamwise normalized vorticity
            vor_y : wall-normal normalized vorticity
            vor_z : spanwise normalized vorticity

    """
    # -------------------------------------------------------------------------------------------------------------------
    # Calculate the dimensional velocity
    # -------------------------------------------------------------------------------------------------------------------
    if "data_type" in data_in.keys():
        data_type = str(data_in["data_type"])                   # definition of the data type.
        if not (data_type=="float32" or data_type=="float16"):
            data_type = "float32"
    else:
        print("[trainvali_data.py:data_traintest_tf] Data type needs to be selected.")
        sys.exit()
    vor_x_norm = np.array(data_in["vor_x_norm"],dtype=data_type)      # velocity fluctuation in the streamwise direction.
    vor_y_norm = np.array(data_in["vor_y_norm"],dtype=data_type)      # velocity fluctuation in the wall-normal direction.
    vor_z_norm = np.array(data_in["vor_z_norm"],dtype=data_type)      # velocity fluctuation in the spanwise direction.
    mean_norm  = bool(data_in["mean_norm"])
    
    # -------------------------------------------------------------------------------------------------------------------
    # import packages
    # -------------------------------------------------------------------------------------------------------------------
    if mean_norm:
        from py_bin.py_functions.normalization_normaldist import read_norm
    else:
        from py_bin.py_functions.normalization import read_norm
    
    # -------------------------------------------------------------------------------------------------------------------
    # Read the normalization parameters, use the float 16 format
    # -------------------------------------------------------------------------------------------------------------------
    folder_data    = str(data_in["folder_data"])                # folder of the generated data
    vornorm_file   = str(data_in["vornorm_file"])                 # file of the normalization data
    vornorm_data   = {"folder":folder_data,"file":vornorm_file}
    try:
        norm_param = read_norm(vornorm_data)
    except:
        print('Normalization file could not be located. Calculation is stopped...',flush=True)
        sys.exit()
    
    if mean_norm:
        print('Mean and std not implemented',flush=True)
        sys.exit()
    else:
        vor_x_max = float(norm_param["vor_x_max"])
        vor_y_max = float(norm_param["vor_y_max"])
        vor_z_max = float(norm_param["vor_z_max"])
        vor_x_min = float(norm_param["vor_x_min"])
        vor_y_min = float(norm_param["vor_y_min"])
        vor_z_min = float(norm_param["vor_z_min"])
        
        # ---------------------------------------------------------------------------------------------------------------
        # Define the normalized fields using float 16 format
        # ---------------------------------------------------------------------------------------------------------------
        vor_x = np.array(vor_x_norm*(vor_x_max-vor_x_min)+vor_x_min,dtype=data_type)
        vor_y = np.array(vor_y_norm*(vor_y_max-vor_y_min)+vor_y_min,dtype=data_type)
        vor_z = np.array(vor_z_norm*(vor_z_max-vor_z_min)+vor_z_min,dtype=data_type)
    data_out = {"vor_x":vor_x,"vor_y":vor_y,"vor_z":vor_z}
    return data_out
