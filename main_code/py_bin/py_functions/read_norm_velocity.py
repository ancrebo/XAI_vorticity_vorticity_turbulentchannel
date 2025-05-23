# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------------------------------------------------------------------
read_norm_velocity.py
-------------------------------------------------------------------------------------------------------------------------
Created on Fri Apr  5 10:39:47 2024

@author: Andres Cremades Botella

File to read and normalize the flow field. 
    Functions:
        - read_norm_velocity : reads the flow and normalize the values
"""
import time

def read_norm_velocity(data_in={"folder":"../../P125_21pi_vu","file":"P125_21pi_vu.$INDEX$.h5.uvw",
                                "padding":15,"shpx":1,"shpy":1,"shpz":1,"dx":1,"dy":1,"dz":1,
                                "data_folder":"Data","umean_file":"Umean.txt","unorm_file":"Unorm.txt",
                                "index":7000,"data_type":"float32","mean_norm":False}):
    """
    .....................................................................................................................
    # read_norm_velocity: Function to read and normalize the velocity
    .....................................................................................................................
    Parameters
    ----------
    data_in : dict, optional
        Data required for reading and normalizing the velocity.
        The default is {"folder":"../../P125_21pi_vu","file":"P125_21pi_vu.$INDEX$.h5.uvw",
                                "padding":15,"shpx":1,"shpy":1,"shpz":1,"dx":1,"dy":1,"dz":1,
                                "data_folder":"Data","umean_file":"Umean.txt","unorm_file":"Unorm.txt",
                                "index":7000,"data_type":"float32","mean_norm":False}.
        Data:
            - folder      : folder to read the data of the velocity fields
            - file        : file to read the data of the velocity fields
            - padding     : padding of the fields
            - shpx        : shape of the fields in x
            - shpy        : shape of the fields in y
            - shpz        : shape of the fields in z
            - dx          : downsampling in x
            - dy          : downsampling in y
            - dz          : downsampling in z
            - data_folder : folder to store the data generated by the code
            - umean_file  : mean velocity file
            - unorm_file  : file for the normalization of the velocity
            - index       : index of the velocity field to read
            - data_type   : definition of the type of data (float32, float16)
            - mean_norm   : flag for normalizing with the mean (True: use mean and std
                                                                False: use min and max)

    Returns
    -------
    dict
        Data containing the normalized velocity and the time required for reading the field and calculating
        the normalization.
        Data:
            - norm_velocity : normalized velocity
            - time_read     : time for reading the file
            - time_norm     : time for normalizing the field

    """
    # -------------------------------------------------------------------------------------------------------------------
    # Load packages
    # -------------------------------------------------------------------------------------------------------------------
    from py_bin.py_functions.read_velocity import read_velocity
    from py_bin.py_functions.norm_velocity import norm_velocity
    
    # -------------------------------------------------------------------------------------------------------------------
    # Read the data
    # -------------------------------------------------------------------------------------------------------------------
    folder      = str(data_in["folder"])                       # folder to read the data of the velocity fields
    file        = str(data_in["file"])                         # file to read the data of the velocity fields
    padding     = int(data_in["padding"])                      # padding of the fields
    shpx        = int(data_in["shpx"])                         # shape of the fields in x
    shpy        = int(data_in["shpy"])                         # shape of the fields in y
    shpz        = int(data_in["shpz"])                         # shape of the fields in z
    dx          = int(data_in["dx"])                           # downsamplin in x
    dy          = int(data_in["dy"])                           # downsamplin in y
    dz          = int(data_in["dz"])                           # downsamplin in z
    data_folder = str(data_in["data_folder"])                  # folder to store generated data
    umean_file  = str(data_in["umean_file"])                   # file to read the mean velocity
    unorm_file  = str(data_in["unorm_file"])                   # file to read the normalization values
    index       = int(data_in["index"])                        # index of the field to read
    mean_norm   = bool(data_in["mean_norm"])                   # flag for using mean and std for normalizing data
    if "data_type" in data_in.keys():
        data_type = str(data_in["data_type"])                  # definition of the data type.
        if not (data_type=="float32" or data_type=="float16"):
            data_type = "float32"
    else:
        print("[trainvali_data.py:data_traintest_tf] Data type needs to be selected.")
        sys.exit()
    
    
    # -------------------------------------------------------------------------------------------------------------------
    # Read and normalize the velocity
    # -------------------------------------------------------------------------------------------------------------------
    tstart             = time.time()
    data_velocity      = {"folder":folder,"file":file,"index":index,"dx":dx,"dy":dy,"dz":dz,
                          "shpx":shpx,"shpy":shpy,"shpz":shpz,"padding":padding,"data_folder":data_folder,
                          "umean_file":umean_file}
    data_read_velocity = read_velocity(data_velocity)
    tread              = time.time()
    data_norm          = {"uu":data_read_velocity["uu"],"vv":data_read_velocity["vv"],
                          "ww":data_read_velocity["ww"],"folder_data":data_folder,"unorm_file":unorm_file,
                          "data_type":data_type,"mean_norm":mean_norm}
    norm_velocity      = norm_velocity(data_norm)
    tnorm              = time.time()
    time_read          = tread-tstart
    time_norm          = tnorm-tread
    
    # -------------------------------------------------------------------------------------------------------------------
    # Return the output
    # -------------------------------------------------------------------------------------------------------------------
    data_out = {"norm_velocity":norm_velocity,"time_read":time_read,"time_norm":time_norm}
    return data_out