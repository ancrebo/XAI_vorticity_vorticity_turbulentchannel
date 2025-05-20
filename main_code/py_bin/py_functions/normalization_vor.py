# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------------------------------------------------------------------
normalization.py
-------------------------------------------------------------------------------------------------------------------------
Created on Fri Mar 22 12:33:24 2024

@author: Andres Cremades Botella

File to create the normalization values for the vorticity fields. The normalization generates 
values between 0 and 1 using the minimum and the maximum of the vorticity values. The file contains
the following functions:
    Functions:
        - save_norm : function for saving the normalization to a file
        - read_norm : function for reading the normalization file
        - calc_norm : function for calculating the normalization
"""

# ---------------------------------------------------------------------------------------------------------------------
# Import packages for all the functions
# ---------------------------------------------------------------------------------------------------------------------
import numpy as np

# ---------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------
# Define the functions
# ---------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------
    
def save_norm(data_in={"folder":"Data","file":"norm.txt","vor_x_max":0,"vor_y_max":0,"vor_z_max":0,"vor_x_min":0,
                       "vor_y_min":0,"vor_z_min":0}):
    """
    .....................................................................................................................
    # save_norm: function for saving the normalization to a file. The function saves the maximum and minimum values
                 of the vorticity components and stress components
    .....................................................................................................................
    Parameters
    ----------
    data_in : dict, optional
        Data for saving the normalization values. 
        The default is {"folder":data,"file":"norm.txt","vor_x_max":0,"vor_y_max":0,"vor_z_max":0,"vor_x_min":0,
                        "vor_y_min":0,"vor_z_min":0}.
        Data:
            - folder    : folder of the generated data
            - file      : file of the normalization data
            - vor_x_max : maximum streamwise vorticity
            - vor_y_max : maximum wall-normal vorticity
            - vor_z_max : maximum spanwise vorticity
            - vor_x_min : minimum streamwise vorticity
            - vor_y_min : minimum wall-normal vorticity
            - vor_z_min : minimum spanwise vorticity
    Returns
    -------
    None.
    """
    
    # -----------------------------------------------------------------------------------------------------------------
    # Read the data
    # -----------------------------------------------------------------------------------------------------------------
    folder    = str(data_in["folder"])    # folder of the normalization data
    file      = str(data_in["file"])      # file of the normalization data
    vor_x_max = float(data_in["vor_x_max"])   # maximum streamwise vorticity
    vor_y_max = float(data_in["vor_y_max"])   # maximum of the wall-normal vorticity 
    vor_z_max = float(data_in["vor_z_max"])   # maximum of the spanwise vorticity
    vor_x_min = float(data_in["vor_x_min"])   # minimum of the streamwise veloctity
    vor_y_min = float(data_in["vor_y_min"])   # minimum of the wall-normal veloctiy
    vor_z_min = float(data_in["vor_z_min"])   # minimum of the spanwise vorticity
    
    # -----------------------------------------------------------------------------------------------------------------
    # Save the data to a file
    # -----------------------------------------------------------------------------------------------------------------
    file_norm = folder+'/'+file
    file_save = open(file_norm, "w+")           
    content = str(vor_x_max)+'\n'
    file_save.write(content)    
    content = str(vor_y_max)+'\n'
    file_save.write(content)    
    content = str(vor_z_max)+'\n'
    file_save.write(content)          
    content = str(vor_x_min)+'\n'
    file_save.write(content)    
    content = str(vor_y_min)+'\n'
    file_save.write(content)    
    content = str(vor_z_min)+'\n'
    file_save.write(content) 
    

def read_norm(data_in={"folder":"Data","file":"norm.txt"}):
    """
    .....................................................................................................................
    # read_norm: function for reading the normalization file. The function reads the maximum and minimum values
                 of the vorticity components and stress components
    .....................................................................................................................
    Parameters
    ----------
    data_in : dict, optional
        Data for the normalization of the vorticity data. 
        The default is {folder:"Data",file:"norm.txt"}.
        Data:
            - folder : folder to read the data
            - file   : file to read the data

    Returns
    -------
    dict
        Data of the maximum and minimum values for the normalization.
        Data:
            - vor_x_max : maximum streamwise vorticity
            - vor_y_max : maximum wall-normal vorticity
            - vor_z_max : maximum spanwise vorticity
            - vor_x_min : minimum streamwise vorticity
            - vor_y_min : minimum wall-normal vorticity
            - vor_z_min : minimum spanwise vorticity

    """
    # -----------------------------------------------------------------------------------------------------------------
    # Read the data
    # -----------------------------------------------------------------------------------------------------------------
    folder = str(data_in["folder"]) # folder to read the normalization data
    file   = str(data_in["file"])   # file to read the normalization data
    
    # -----------------------------------------------------------------------------------------------------------------
    # Read the normalization file
    # -----------------------------------------------------------------------------------------------------------------
    file_norm = folder+'/'+file
    file_read = open(file_norm,"r")
    vor_x_max = np.array(file_read.readline().replace('[','').replace(']','').split(','),dtype='float')
    vor_y_max = np.array(file_read.readline().replace('[','').replace(']','').split(','),dtype='float')
    vor_z_max = np.array(file_read.readline().replace('[','').replace(']','').split(','),dtype='float')
    vor_x_min = np.array(file_read.readline().replace('[','').replace(']','').split(','),dtype='float')
    vor_y_min = np.array(file_read.readline().replace('[','').replace(']','').split(','),dtype='float')
    vor_z_min = np.array(file_read.readline().replace('[','').replace(']','').split(','),dtype='float')
    data_out = {"vor_x_max":vor_x_max,"vor_y_max":vor_y_max,"vor_z_max":vor_z_max,"vor_x_min":vor_x_min,
                "vor_y_min":vor_y_min,"vor_z_min":vor_z_min}
    return data_out

               
def calc_norm(data_in={"field_ini":1000,"field_fin":9999,"data_folder":"Data",
                       "dx":1,"dy":1,"dz":1,"folder":"../../P125_21pi_vu","file":"P125_21pi_vu.$INDEX$.h5.uvw",
                       "shpx":192,"shpy":201,"shpz":96,"save_file":True,"vornorm_file":"norm.txt",
                       "vormean_file":"vormean.txt"}):
    """
    .....................................................................................................................
    # calc_norm: function to calculate the normalization of the vorticity. The function calculates the maximum and 
                 minimum values of the vorticity components and stress components
    .....................................................................................................................
    Parameters
    ----------
    data_in : dict, optional
        Data for normalizing the vorticity.
        The default is {"field_ini":1000,"field_fin":9999,"data_folder":"Data","umean_file"="Umean.txt",
                        "dx":1,"dy":1,"dz":1,"folder":"../P125_21pi_vu","file":"P125_21pi_vu.$INDEX$.h5.uvw",
                        "shpx":192,"shpy":201,"shpz":96,"padding":15,"save_file":True,"vornorm_file":"norm.txt",
                        "vormean_file":"vormean.txt"}.
        Data:
            - field_ini    : initial field of the data to calculate the normalization
            - field_fin    : final field of the data to calculate the normalization
            - data_folder  : folder to store the data calculated by the code
            - dx           : downsampling of x direction
            - dy           : downsampling of y direction
            - dz           : downsampling of z direction
            - folder       : folder of the vorticity data
            - file         : file of the vorticity data without index
            - shpx         : shape of the tensor in x
            - shpy         : shape of the tensor in y
            - shpz         : shape of the tensor in z
            - save_file    : flag to save the normalization in a file
            - vornorm_file : file of the normalization data
            - vormean_file : file of the mean data
            

    Returns
    -------
    dict
        Data for the normalization. Only returns it in case of not saving a file
        Data:
            - vor_x_max : maximum streamwise vorticity
            - vor_y_max : maximum wall-normal vorticity
            - vor_z_max : maximum spanwise vorticity
            - vor_x_min : minimum streamwise vorticity
            - vor_y_min : minimum wall-normal vorticity
            - vor_z_min : minimum spanwise vorticity
    """
    # -----------------------------------------------------------------------------------------------------------------
    # Import packages
    # -----------------------------------------------------------------------------------------------------------------
    from py_bin.py_functions.read_vorticity import read_vorticity
    
    # -----------------------------------------------------------------------------------------------------------------
    # Read the data
    # -----------------------------------------------------------------------------------------------------------------
    field_ini    = int(data_in["field_ini"])   # initial field for calculating the normalization
    field_fin    = int(data_in["field_fin"])   # final field for calculating the normalization
    data_folder  = str(data_in["data_folder"]) # folder for reading the data
    folder       = str(data_in["folder"])      # folder to read the vorticity fields
    file         = str(data_in["file"])        # file to read the veloctity fields
    dx           = int(data_in["dx"])          # downsampling in x
    dy           = int(data_in["dy"])          # downsampling in y
    dz           = int(data_in["dz"])          # downsampling in z
    shpx         = int(data_in["shpx"])        # shape in the x direction
    shpy         = int(data_in["shpy"])        # shape in the y direction
    shpz         = int(data_in["shpz"])        # shape in the z direction
    save_file    = bool(data_in["save_file"])  # flag to decide if the normalization must be save in a file
    vornorm_file = str(data_in["vornorm_file"])  # file to save the normalization
    vormean_file = str(data_in["vormean_file"])
    
    # -----------------------------------------------------------------------------------------------------------------
    # In the loop
    #   - ii : index of the file that we are reading
    # -----------------------------------------------------------------------------------------------------------------
    for ii in range(field_ini,field_fin):
        print("Reading Field "+str(ii),flush=True)
        data_vorticity      = {"folder":folder,"file":file,"index":ii,"dx":dx,"dy":dy,"dz":dz,"shpx":shpx,
                               "shpy":shpy,"shpz":shpz,"padding":0,"data_folder":data_folder,
                               "vormean_file":vormean_file}            
        data_read_vorticity = read_vorticity(data_vorticity)
        vor_x_i0 = np.array(data_read_vorticity['vor_x'],dtype='float')
        vor_y_i0 = np.array(data_read_vorticity['vor_y'],dtype='float')
        vor_z_i0 = np.array(data_read_vorticity['vor_z'],dtype='float')
        if ii == field_ini:
            vor_x_max = np.max(vor_x_i0)
            vor_y_max = np.max(vor_y_i0)
            vor_z_max = np.max(vor_z_i0)
            vor_x_min = np.min(vor_x_i0)
            vor_y_min = np.min(vor_y_i0)
            vor_z_min = np.min(vor_z_i0)
        else:
            vor_x_max = np.max([vor_x_max,np.max(vor_x_i0)])
            vor_y_max = np.max([vor_y_max,np.max(vor_y_i0)])
            vor_z_max = np.max([vor_z_max,np.max(vor_z_i0)])
            vor_x_min = np.min([vor_x_min,np.min(vor_x_i0)])
            vor_y_min = np.min([vor_y_min,np.min(vor_y_i0)])
            vor_z_min = np.min([vor_z_min,np.min(vor_z_i0)])
            
    # -----------------------------------------------------------------------------------------------------------------
    # Save the normalization in a file or return the values of the normalization
    # -----------------------------------------------------------------------------------------------------------------
    if save_file:
        data_norm_save = {"folder":data_folder,"file":vornorm_file,"vor_x_max":vor_x_max,"vor_y_max":vor_y_max,
                          "vor_z_max":vor_z_max,"vor_x_min":vor_x_min,"vor_y_min":vor_y_min,"vor_z_min":vor_z_min}
        save_norm(data_in=data_norm_save)
    else:
        data_out = {"vor_x_max":vor_x_max,"vor_y_max":vor_y_max,"vor_z_max":vor_z_max,"vor_x_min":vor_x_min,
                    "vor_y_min":vor_y_min,"vor_z_min":vor_z_min}
        return data_out
    
        