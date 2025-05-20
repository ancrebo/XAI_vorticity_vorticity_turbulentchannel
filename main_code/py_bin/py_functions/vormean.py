# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------------------------------------------------------------------
umean.py
-------------------------------------------------------------------------------------------------------------------------
Created on Fri Mar 22 09:23:59 2024

@author: Andres Cremades Botella

File for calculating and reading the mean velocity. The file contains the following functions:
    Functions:
        - save_Umean : function for saving the mean velocity
        - read_Umean : function for reading the mean velocity
        - calc_Umean : function for calculating the mean velocity
"""

# -----------------------------------------------------------------------------------------------------------------------
# Import packages for all the functions
# ----------------------------------------------------------------------------------------------------------------------- 
import numpy as np

# -----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------
# Define the functions
# -----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------

def save_vormean(data_in={"folder":"Data","file":"Umean.txt","vor_x_mean":[],"vor_y_mean":[],"vor_z_mean":[]}):
    """
    .....................................................................................................................
    # save_vormean: Function for saving the mean vorticity
    .....................................................................................................................
    Parameters
    ----------
    data_in : dict, optional
        Information for saving the mean velocity.
        The default is {"folder":"Data","file":"Umean.txt","vor_x_mean":[],"vor_y_mean":[],"vor_z_mean":[]}.
        Data:
            - folder     : folder to save the information
            - file       : file to save the information
            - vor_x_mean : mean streamwise velocity
            - vor_y_mean : mean wall-normal velocity
            - vor_z_mean : mean spanwise velocity

    Returns
    -------
    None.

    """
    # -------------------------------------------------------------------------------------------------------------------
    # Read data
    # -------------------------------------------------------------------------------------------------------------------
    folder     = str(data_in["folder"])                    # Folder to save teh mean velocity
    file       = str(data_in["file"])                      # File to save the mean velocity
    vor_x_mean = np.array(data_in["vor_x_mean"],dtype='float') # Mean velocity in the streamwise direction
    vor_y_mean = np.array(data_in["vor_y_mean"],dtype='float') # Mean velocity in the wall-normal direction
    vor_z_mean = np.array(data_in["vor_z_mean"],dtype='float') # Mean velocity in the spanwise direction 

    # -------------------------------------------------------------------------------------------------------------------
    # Save the information in a file
    # -------------------------------------------------------------------------------------------------------------------    
    file_vormean = folder+'/'+file                     
    file_save    = open(file_vormean, "w+")           
    content      = str(vor_x_mean.tolist())+'\n'
    file_save.write(content)          
    content      = str(vor_y_mean.tolist())+'\n'
    file_save.write(content)          
    content      = str(vor_z_mean.tolist())+'\n'
    file_save.write(content)
    file_save.close()
    
    
def read_vormean(data_in={"folder":"Data","file":"vormean.txt","dy":1}):
    """ 
    .....................................................................................................................   
    # read_vormean: Function for read the mean vorticity
    .....................................................................................................................
    Parameters
    ----------
    data_in : dict, optional
        Data to read the mean velocity.
        The default is {"folder":"Data","file":"vormean.txt","dy":1}.
        Data:
            - folder : folder to read the calculated data
            - file   : file of the mean velocity
            - dy     : downsampling in the y direction

    Returns
    -------
    dict
        Mean velocity in the streamwise, wall-normal and spanwise directions is returned.
        Data:
            - vor_x_mean : streamwise mean vorticity as a function of the wall-normal distance
            - vor_y_mean : wall-normal mean vorticity as a function of the wall-normal distance
            - vor_z_mean : spanwise mean vorticity as a function of the wall-normal distance

    """
    # -------------------------------------------------------------------------------------------------------------------
    # Read data
    # -------------------------------------------------------------------------------------------------------------------
    folder    = str(data_in["folder"]) # folder for reading the mean velocity data
    file      = str(data_in["file"])   # file for reading the mean velocity data
    dy        = int(data_in["dy"])     # downsampling in the wall-normal direction 
    
    # -------------------------------------------------------------------------------------------------------------------
    # Read the velocity from file
    # -------------------------------------------------------------------------------------------------------------------
    filevormean = folder+'/'+file
    file_read   = open(filevormean,"r")
    vor_x_mean  = np.array(file_read.readline().replace('[','').replace(']','').split(','),dtype='float')[::dy]
    vor_y_mean  = np.array(file_read.readline().replace('[','').replace(']','').split(','),dtype='float')[::dy]
    vor_z_mean  = np.array(file_read.readline().replace('[','').replace(']','').split(','),dtype='float')[::dy]
    data_out  = {"vor_x_mean":vor_x_mean,"vor_y_mean":vor_y_mean,"vor_z_mean":vor_z_mean}
    return data_out
    
    
def calc_vormean(data_in={"field_ini":1000,"field_fin":9999,"folder":"../../P125_21pi_vu",\
                          "file":"P125_21pi_vu.$INDEX$.h5.uvw","save_file":True,"vormean_file":"vormean.txt",
                          "data_folder":"Data","shpx":192,"shpy":201,"shpz":96}):
    """
    .....................................................................................................................
    # calc_vormean: Function for calculating the mean vorticity
    .....................................................................................................................
    Parameters
    ----------
    data_in : dict, optional
        Data for the calculation of the mean velocity. 
        The default is {"field_ini":1000,"field_fin":9999,"folder":"../P125_21pi_vu",
                        "file":"P125_21pi_vu.$INDEX$.h5.uvw","save_file":True,"vormean_file":"vormean.txt",
                        "data_folder":"Data"}.
        Data:
            - field_ini    : index of the initial field used for calculating the mean velocity
            - field_fin    : index of the final field used for calculating the mean velocity
            - folder       : path of the folder to read the velocity data base
            - file         : name of the file to read the velocity
            - save_file    : flag to save the mean velocity in a file
            - vormean_file : file for saving the mean velocity
            - data_folder  : path of the folder of the data calculated by the code
            - shpx         : shape of the tensors in the streamwise direction
            - shpy         : shape of the tensors in the wall-normal direction
            - shpz         : shape of the tensors in the spanwise direction

    Returns
    -------
    dict
        Mean vorticity in the streamwise, wall-normal and spanwise directions is returned. Only used when the saving
        option is not active.
        Data:
            - vor_x_mean : streamwise mean vorticity as a function of the wall-normal distance
            - vor_y_mean : wall-normal mean vorticity as a function of the wall-normal distance
            - vor_z_mean : spanwise mean vorticity as a function of the wall-normal distance

    """
    # -------------------------------------------------------------------------------------------------------------------
    # Load packages
    # -------------------------------------------------------------------------------------------------------------------
    import h5py
    
    # -------------------------------------------------------------------------------------------------------------------
    # Read data
    # -------------------------------------------------------------------------------------------------------------------
    field_ini    = int(data_in["field_ini"])    # index of the initial field used for calculating the mean velocity
    field_fin    = int(data_in["field_fin"])    # index of the final field used for calculating the mean velocity
    folder       = str(data_in["folder"])       # path to the folder for reading the velocity data base
    file         = str(data_in["file"])         # name of the file containing the velocity data
    file_comp    = folder+'/'+file
    save_file    = bool(data_in["save_file"])   # flag for saving the file
    vormean_file = str(data_in["vormean_file"]) # file for the mean velocity
    data_folder  = str(data_in["data_folder"])  # folder of the data calculated by the code
    shpx         = int(data_in["shpx"])         # shape of the tensors in the streamwise direction
    shpy         = int(data_in["shpy"])         # shape of the tensors in the wall-normal direction
    shpz         = int(data_in["shpz"])         # shape of the tensors in the spanwise direction
     
    # -------------------------------------------------------------------------------------------------------------------
    # In the following lines:
    #     - ii : index related to the field
    #     - file_ii : file of the velocity field including the index
    #     - UU      : streamwise velocity
    #     - VV      : wall-normal velocity
    #     - WW      : spanwise velocity
    #     - UU_cum  : cumulative streamwise velocity
    #     - VV_cum  : cumulative wall-normal velocity
    #     - WW_cum  : cumulative spanwise velocity
    #     - nn_cum  : number of gridpoints used for calculate the mean velocity
    # -------------------------------------------------------------------------------------------------------------------
    for ii in range(field_ini,field_fin):            
        file_ii = file_comp.replace("$INDEX$",str(ii))
        print('Mean vorticity calculation:' + str(file_ii),flush=True)
        file_r = h5py.File(file_ii,'r')
        vor_x = np.array(file_r['vorx'])
        vor_y = np.array(file_r['vory'])
        vor_z = np.array(file_r['vorz'])
        if ii == field_ini:
            vor_x_cum = np.sum(vor_x,axis=(1,2))
            vor_y_cum = np.sum(vor_y,axis=(1,2))
            vor_z_cum = np.sum(vor_z,axis=(1,2))
            nn_cum    = np.ones((shpy,))*shpx*shpz
        else:
            vor_x_cum += np.sum(vor_x,axis=(1,2))
            vor_y_cum += np.sum(vor_y,axis=(1,2))
            vor_z_cum += np.sum(vor_z,axis=(1,2))
            nn_cum += np.ones((shpy,))*shpx*shpz
            
    # -------------------------------------------------------------------------------------------------------------------
    # Calculate the mean velocity
    #     - UUmean : mean velocity in the streamwise velocity
    #     - VVmean : mean velocity in the wall-normal velocity
    #     - WWmean : mean velocity in the spanwise veloctiy
    # -------------------------------------------------------------------------------------------------------------------
    vor_x_mean = np.divide(vor_x_cum,nn_cum)
    vor_y_mean = np.divide(vor_y_cum,nn_cum)
    vor_z_mean = np.divide(vor_z_cum,nn_cum)
    
    # -------------------------------------------------------------------------------------------------------------------
    # Save the mean in a file or return the mean velocity
    # -------------------------------------------------------------------------------------------------------------------
    if save_file:
        data_vormean = {"folder":data_folder,"file":vormean_file,"vor_x_mean":vor_x_mean,"vor_y_mean":vor_y_mean,
                        "vor_z_mean":vor_z_mean}
        save_vormean(data_in=data_vormean)
    else:
        data_out  = {"vor_x_mean":vor_x_mean,"vor_y_mean":vor_y_mean,"vor_z_mean":vor_z_mean}
        return data_out
            
    