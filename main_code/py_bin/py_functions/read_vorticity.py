# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------------------------------------------------------------------
read_velocity.py
-------------------------------------------------------------------------------------------------------------------------
Created on Thu Mar 21 15:18:38 2024

@author: Andres Cremades Botella

File to read the data of the velocity fields. The file contains the following functions:
    Functions:
        - read_velocity : file to read the velocity
"""
# -----------------------------------------------------------------------------------------------------------------------
# Read the packages for all the functions
# -----------------------------------------------------------------------------------------------------------------------
import sys
import numpy as np

# -----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------
# Define the functions
# -----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------

def read_vorticity(data_in={"folder":"../../P125_21pi_vu","file":"P125_21pi_vu.$INDEX$.h5.uvw","index":1000,
                            "dx":1,"dy":1,"dz":1,"shpx":196,"shpy":201,"shpz":96,"padding":15,
                            "data_folder":"Data","vormean_file":"vormean.txt"}):
    """
    .....................................................................................................................
    # read_velocity: Function to read the velocity, the function generates three arrays containing the velocities in
                     the streamwise (uu), wall-normal (vv) and spanwise (ww) directions
    .....................................................................................................................
    Parameters
    ----------
    data_in : TYPE, optional
        DESCRIPTION. The default is {"folder":"../P125_21pi_vu","file":"P125_21pi_vu.$INDEX$.h5.uvw","index":1000,
                                     "dx":1,"dy":1,"dz":1,"shpx":196,"shpy":201,"shpz":96,"padding":15
                                     "data_folder":"Data","vormean_file":"vormean.txt"}.
        Data:
            - folder       : folder of the velocity data
            - file         : file of the velocity data without the index
            - index        : index of the velocity data file
            - dx           : downsampling in x
            - dy           : downsampling in y
            - dz           : downsampling in z
            - shpx         : shape in x of the tensors
            - shpy         : shape in y of the tensors
            - shpz         : shape in z of the tensors
            - padding      : padding of the fields
            - data_folder  : folder to store generated data
            - vormean_file : file containing the mean vorticity
    Returns
    -------
    dict
        Velocity fluctuation in the streamwise, wall-normal and spanwise directions.
        Data:
            - vor_x : vorticity fluctuation in the streamwise direction
            - vor_y : vorticity fluctuation in the wall-normal direction
            - vor_z : vorticity fluctuation in the spanwise direction

    """
    # -------------------------------------------------------------------------------------------------------------------
    # Import packages
    # -------------------------------------------------------------------------------------------------------------------
    import h5py
    from py_bin.py_functions.vormean import read_vormean
    from py_bin.py_functions.padding_field import padding_field
    
    # -------------------------------------------------------------------------------------------------------------------
    # Read the data
    # -------------------------------------------------------------------------------------------------------------------
    folder        = str(data_in["folder"])      # folder of the data
    file          = str(data_in["file"])        # file of the data
    index         = str(int(data_in["index"]))  # index of the file
    dx            = int(data_in["dx"])          # downsampling in x
    dy            = int(data_in["dy"])          # downsampling in y
    dz            = int(data_in["dz"])          # downsampling in z
    shpx          = int(data_in["shpx"])        # shape of the tensors in the x direction
    shpy          = int(data_in["shpy"])        # shape of the tensors in the y direction
    shpz          = int(data_in["shpz"])        # shape of the tensors in the z direction
    padding       = int(data_in["padding"])     # padding of the fields
    data_folder   = str(data_in["data_folder"]) # folder to store generated data
    vormean_file  = str(data_in["vormean_file"]) # file of the mean vorticity
    file_complete = folder+'/'+file
    file_ii       = file_complete.replace("$INDEX$",index) 
    try:
        datavormean  = {"folder":data_folder,"file":vormean_file,"dy":dy}
        meanvor_data = read_vormean(datavormean)
    except:
        print("Mean velocity file needs to be provided. Breaking calculation...",flush=True)
        sys.exit()
        
    # -------------------------------------------------------------------------------------------------------------------
    # Read the mean velocity in the streamwise direction, other directions have null mean velocity
    # -------------------------------------------------------------------------------------------------------------------
    vor_x_mean  = meanvor_data["vor_x_mean"]
    vor_y_mean  = meanvor_data["vor_y_mean"]
    vor_z_mean  = meanvor_data["vor_z_mean"]
    print('Reading field:' + str(file_ii),flush=True)
        
    
    # -------------------------------------------------------------------------------------------------------------------
    # Read the information from the files. The information requires the mean values in the wall-normal directions
    # This values should be stored in the data folder, in the case of missing the file, the software will calculate it.
    # -------------------------------------------------------------------------------------------------------------------
    file = h5py.File(file_ii,'r')    
    VOR_x = np.array(file['vorx'])[::dy,::dz,::dx]
    vor_x = VOR_x-vor_x_mean.reshape(-1,1,1)
    VOR_y = np.array(file['vory'])[::dy,::dz,::dx]
    vor_y = VOR_y-vor_y_mean.reshape(-1,1,1)
    VOR_z = np.array(file['vorz'])[::dy,::dz,::dx]
    vor_z = VOR_z-vor_z_mean.reshape(-1,1,1)
    
    # -------------------------------------------------------------------------------------------------------------------
    # Apply the padding if necessary. The padding takes the variable padding to add that number of nodes in both sizes
    # of the channel in the streamwise and the spanwise directions. The idea is to preserve the periodicity of the 
    # channel
    # -------------------------------------------------------------------------------------------------------------------
    if padding > 0:
        vor_x = padding_field(data_in={"field":vor_x,"shpx":shpx,"shpy":shpy,"shpz":shpz,"padding":padding})["field"]
        vor_y = padding_field(data_in={"field":vor_y,"shpx":shpx,"shpy":shpy,"shpz":shpz,"padding":padding})["field"]
        vor_z = padding_field(data_in={"field":vor_z,"shpx":shpx,"shpy":shpy,"shpz":shpz,"padding":padding})["field"]
    data_output = {"vor_x":vor_x,"vor_y":vor_y,"vor_z":vor_z}
    return data_output