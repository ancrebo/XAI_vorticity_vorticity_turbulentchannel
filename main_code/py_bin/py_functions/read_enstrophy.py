# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------------------------------------------------------------------
read_enstrophy.py
-------------------------------------------------------------------------------------------------------------------------
Created on Thu Mar 21 15:18:38 2024

@author: Andres Cremades Botella

File to read the data of the velocity fields. The file contains the following functions:
    Functions:
        - read_enstrophy : file to read the enstrophy
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

def read_enstrophy(data_in={"folder":"../../P125_21pi_vu","file":"P125_21pi_vu.$INDEX$.h5.uvw","index":1000,
                            "dx":1,"dy":1,"dz":1,"shpx":196,"shpy":201,"shpz":96,"padding":15}):
    """
    .....................................................................................................................
    # read_enstrophy: Function to read the enstrophy, the function generates one array containing the enstrophy
    .....................................................................................................................
    Parameters
    ----------
    data_in : TYPE, optional
        DESCRIPTION. The default is {"folder":"../../P125_21pi_vu","file":"P125_21pi_vu.$INDEX$.h5.uvw","index":1000,
                                     "dx":1,"dy":1,"dz":1,"shpx":196,"shpy":201,"shpz":96,"padding":15}.
        Data:
            - folder      : folder of the velocity data
            - file        : file of the velocity data without the index
            - index       : index of the velocity data file
            - dx          : downsampling in x
            - dy          : downsampling in y
            - dz          : downsampling in z
            - shpx        : shape in x of the tensors
            - shpy        : shape in y of the tensors
            - shpz        : shape in z of the tensors
            - padding     : padding of the fields
    Returns
    -------
    dict
        Velocity fluctuation in the streamwise, wall-normal and spanwise directions.
        Data:
            - ens : enstrophy data

    """
    # -------------------------------------------------------------------------------------------------------------------
    # Import packages
    # -------------------------------------------------------------------------------------------------------------------
    import h5py
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
    file_complete = folder+'/'+file
    file_ii        = file_complete.replace("$INDEX$",index) 
        
    
    # -------------------------------------------------------------------------------------------------------------------
    # Read the information from the files. 
    # -------------------------------------------------------------------------------------------------------------------
    file = h5py.File(file_ii,'r')    
    ens  = np.array(file['ens'])[::dy,::dz,::dx]
    
    # -------------------------------------------------------------------------------------------------------------------
    # Apply the padding if necessary. The padding takes the variable padding to add that number of nodes in both sizes
    # of the channel in the streamwise and the spanwise directions. The idea is to preserve the periodicity of the 
    # channel
    # -------------------------------------------------------------------------------------------------------------------
    if padding > 0:
        ens = padding_field(data_in={"field":ens,"shpx":shpx,"shpy":shpy,"shpz":shpz,"padding":padding})["field"]
    data_output = {"ens":ens}
    return data_output