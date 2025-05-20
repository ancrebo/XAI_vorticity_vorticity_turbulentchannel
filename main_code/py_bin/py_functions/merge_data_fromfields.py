"""
-------------------------------------------------------------------------------------------------------------------------
merge_data.py
-------------------------------------------------------------------------------------------------------------------------
Created on Wed May 15 11:58:23 2024

@author: Andres Cremades Botella - Inherited from Xuan Gu (Berzelius support)

File containing the function to improve the reading of the data in Berzelius. The file contains the following functions:
    Functions:
        - load_datasets : function for loading the data
        - read_data_tf           : function to read the data with the tensorflow format
        - read_inout_notprepared : function to read the data for the training and the test directly from the flow
                                   files
        - data_traintest_tf      : function to separe training and test data from the flow files and convert them 
                                   into the tensorflow format
"""
# -----------------------------------------------------------------------------------------------------------------------
# Import packages
# -----------------------------------------------------------------------------------------------------------------------
import sys
import tensorflow as tf
import os
import numpy as np
from tqdm import tqdm

def create_datasets(data_in={"folder_uvw":"../../P125_21pi_vu/","file_uvw":'P125_21pi_vu.1000.h5.uvw',
                             "folder_vor":"../../P125_21pi_vu/","file_vor":'P125_21pi_vu.1000.h5.uvw',"elem_spec":[],
                             "padding":15,"shpx":1,"shpy":1,"shpz":1,"dx":1,"dy":1,"dz":1,"data_folder":"Data",
                             "umean_file":"umean.txt","unorm_file":"norm.txt","vornorm_file":"norm.txt",
                             "data_type":"float32","mean_norm":False,
                             "index":1000,"delta_pred":1,"vormean_file":"vormean.txt"}):
    """
    .....................................................................................................................
    # create_datasets: Function for reading the data
    .....................................................................................................................
    Parameters
    ----------
    data_in : dict, optional
        data required for selecting the fields for the training process.
        The default is {"folder_uvw":"../../P125_21pi_vu/","file_uvw":'P125_21pi_vu.1000.h5.uvw',
                        "folder_vor":"../../P125_21pi_vu/","file_vor":'P125_21pi_vu.1000.h5.uvw',"elem_spec":[],
                        "padding":15,"shpx":1,"shpy":1,"shpz":1,"dx":1,"dy":1,"dz":1,"data_folder":"Data",
                        "umean_file":"umean.txt","unorm_file":"norm.txt","vornorm_file":"norm.txt",
                        "data_type":"float32","mean_norm":False,
                        "index":1000,"delta_pred":1,"vormean_file":"vormean.txt"}.
        Data:
            - folder_uvw   : path to the folder to read the fields
            - file_uvw     : path to the file to read the fields
            - folder_vor   : path to the folder to read the vorticity fields
            - file_vor     : path to the file to read the vorticity fields
            - elem_spec    : specification of the format of the elements
            - version_tf   : version of the tensorflow
            - padding      : padding applied to the fields
            - shpx         : shape in the streamwise direction
            - shpy         : shape in the wall-normal direction
            - shpz         : shape in the spanwise direction
            - dx           : downsampling in the streamwise direction
            - dy           : downsampling in the wall-normal direction
            - dz           : downsampling in the spanwise direction
            - data_folder  : folder to store the data of the model
            - umean_file   : file storing the mean velocity
            - unorm_file   : file storing the normalization
            - vornorm_file : file storing the normalization of the vorticity
            - data_type    : type of data used (float32,float16)
            - mean_norm    : flag for choosing between standarization and normalization (True: standarization with
                                                                                         mean and standard deviation,
                                                                                         False: normalization 
                                                                                         between minimum and maximum)
            - delta_pred   : number of fields to advance in the predictions
            - vormean_file : file of the mean vorticity
    
    Returns
    -------
    dict
        Structure containing the data for the training and the data for the test
        Data:
            - data_XY : input and output data
    """
    # -------------------------------------------------------------------------------------------------------------------
    # Import modules
    # -------------------------------------------------------------------------------------------------------------------
    from py_bin.py_functions.read_norm_velocity import read_norm_velocity
    from py_bin.py_functions.read_norm_vorticity import read_norm_vorticity
    
    # -------------------------------------------------------------------------------------------------------------------
    # Read data
    # -------------------------------------------------------------------------------------------------------------------
    file_uvw     = str(data_in["file_uvw"])
    folder_uvw   = str(data_in["folder_uvw"])
    file_vor     = str(data_in["file_vor"])
    folder_vor   = str(data_in["folder_vor"])
    elem_spec    = data_in["elem_spec"]
    shpx         = int(data_in["shpx"])
    shpy         = int(data_in["shpy"])
    shpz         = int(data_in["shpz"])
    dx           = int(data_in["dx"])
    dy           = int(data_in["dy"])
    dz           = int(data_in["dz"])
    padding      = int(data_in["padding"])
    data_folder  = str(data_in["data_folder"])
    umean_file   = str(data_in["umean_file"])
    unorm_file   = str(data_in["unorm_file"])
    vornorm_file = str(data_in["vornorm_file"])
    mean_norm    = bool(data_in["mean_norm"])
    index        = int(data_in["index"])
    delta_pred   = int(data_in["delta_pred"])
    data_type    = str(data_in["data_type"])
    vormean_file = str(data_in["vormean_file"])
    version_tf   = np.array(tf.__version__.split('.'),dtype="int")
    data_X       = np.zeros((1,shpy,shpz+2*padding,shpx+2*padding,3),dtype=data_type)
    data_Y       = np.zeros((1,shpy,shpz,shpx,3),dtype=data_type)
    
    # ---------------------------------------------------------------------------------------------------------------
    # Read the input field
    # ---------------------------------------------------------------------------------------------------------------
    data_norm_X      = {"folder":folder_vor,"file":file_vor,"padding":padding,"shpx":shpx,
                        "shpy":shpy,"shpz":shpz,"dx":dx,"dy":dy,"dz":dz,"data_folder":data_folder,
                        "vornorm_file":vornorm_file,"index":index,"data_type":data_type,
                        "mean_norm":mean_norm,"vormean_file":vormean_file}
    data_vorti_norm  = read_norm_vorticity(data_in=data_norm_X)
    norm_vorticity_X = data_vorti_norm["norm_vorticity"]
    print("Time for reading the field: "+str(data_vorti_norm["time_read"]),flush=True)
    print("Time for normalizing the field: "+str(data_vorti_norm["time_norm"]),flush=True)
    del data_norm_X,data_vorti_norm
    
    # ---------------------------------------------------------------------------------------------------------------
    # Read the output field
    # ---------------------------------------------------------------------------------------------------------------
    data_norm_Y      = {"folder":folder_vor,"file":file_vor,"padding":0,"shpx":shpx,
                        "shpy":shpy,"shpz":shpz,"dx":dx,"dy":dy,"dz":dz,"data_folder":data_folder,
                        "vornorm_file":vornorm_file,"index":index+delta_pred,"data_type":data_type,
                        "mean_norm":mean_norm,"vormean_file":vormean_file}
    data_vorti_norm  = read_norm_vorticity(data_in=data_norm_Y)
    norm_vorticity_Y = data_vorti_norm["norm_vorticity"]
    print("Time for reading the field: "+str(data_vorti_norm["time_read"]),flush=True)
    print("Time for normalizing the field: "+str(data_vorti_norm["time_norm"]),flush=True)
    del data_norm_Y,data_vorti_norm
    print("-"*100,flush=True)
    
    # ---------------------------------------------------------------------------------------------------------------
    # Store the input and output data
    # ---------------------------------------------------------------------------------------------------------------
    data_X[0,:,:,:,0] = norm_vorticity_X['vor_x_norm']
    data_X[0,:,:,:,1] = norm_vorticity_X['vor_y_norm']
    data_X[0,:,:,:,2] = norm_vorticity_X['vor_z_norm']
    del norm_vorticity_X
    data_Y[0,:,:,:,0] = norm_vorticity_Y['vor_x_norm']
    data_Y[0,:,:,:,1] = norm_vorticity_Y['vor_y_norm']
    data_Y[0,:,:,:,2] = norm_vorticity_Y['vor_z_norm']
    del norm_velocity_Y
    print("-"*100,flush=True)
    
    data_XY = tf.data.Dataset.from_tensor_slices((data_X,data_Y))
    data_out = {"data_XY":data_XY}
    return data_out

def serialize_example(data_in={"feature":[],"label":[]}):
    """
    .....................................................................................................................
    # _serialize_example: Create a tf.train.Example message ready to be written to a file.
    .....................................................................................................................
    Parameters
    ----------
    data_in : dict, optional
        data containing the input information.
        The default is {"feature":[],"label":[]}.
        Data:
            - feature : features of the tfrecord
            - label   : name of the tfrecord
    
    Returns
    -------
    dict
        Structure containing mapped example
        Data:
            - example : Structure containing mapped example
    """
    # -------------------------------------------------------------------------------------------------------------------
    # Read the data
    # -------------------------------------------------------------------------------------------------------------------
    feature = data_in["feature"]
    label   = data_in["label"]
    
    # -------------------------------------------------------------------------------------------------------------------
    # Create a dictionary mapping the feature name to the tf.train.Example-compatible data type.
    # -------------------------------------------------------------------------------------------------------------------
    feature = {'feature': tf.train.Feature(float_list=tf.train.FloatList(value=feature.flatten())),
               'label': tf.train.Feature(float_list=tf.train.FloatList(value=label.flatten()))}

    # -------------------------------------------------------------------------------------------------------------------
    # Create a Features message using tf.train.Example.
    # -------------------------------------------------------------------------------------------------------------------
    example_proto = tf.train.Example(features=tf.train.Features(feature=feature))
    data_out      = {"example":example_proto.SerializeToString()}
    return data_out

def write_tfrecord(data_in={"output_path":'/tfrecord/merged_vali_data/dataset_0000.tfrecord',"dataset":[]}):
    """
    ---------------------------------------------------------------------------------------------------------------------
    # write_tfrecord : Writes dataset to a TFRecord file.
    ---------------------------------------------------------------------------------------------------------------------
    Parameters
    ----------
    data_in : TYPE, optional
        DESCRIPTION. The default is {"file_path":'/tfrecord/merged_vali_data/dataset_0000.tfrecord',
                                     "dataset":[]}.
        Data:
            - file      : path to the file to write
            - dataset   : Data to save

    Returns
    -------
    None.

    """
    output_path = str(data_in["output_path"])
    dataset     = data_in["dataset"]
    with tf.io.TFRecordWriter(output_path) as writer:
        for features, labels in tqdm(dataset, desc=f"Writing {output_path}"):
            data_serialize = {"feature":features.numpy(),"label":labels.numpy()}
            example        = serialize_example(data_in=data_serialize)["example"]
            writer.write(example)
    

def create_tfrecords(data_in={"output_directory":'/tfrecord/',"base_directory_uvw":"../../P125_21pi_vu_tf_float32/",
                              "base_file_uvw":"P125_21pi_vu.$INDEX$.h5.uvw",
                              "base_directory_vor":"../../P125_21pi_vu_tf_float32/",
                              "base_file_vor":"P125_21pi_vu.$INDEX$.h5.uvw","datasets":[],"elem_spec":[],"shpx":1,
                              "shpy":1,"shpz":1,"dx":1,"dy":1,"dz":1,"padding":15,"data_type":"float32",
                              "data_folder":"Data","umean_file":"umean.txt","unorm_file":"norm.txt",
                              "vornorm_file":"norm.txt","mean_norm":False,
                              "delta_pred":1,"vormean_file":"vormean.txt"}):
    """
    ---------------------------------------------------------------------------------------------------------------------
    # create_tfrecords : Function to create the tfrecords.
    ---------------------------------------------------------------------------------------------------------------------
    Parameters
    ----------
    data_in : TYPE, optional
        DESCRIPTION. The default is {"output_directory":'/tfrecord/',"base_directory_uvw":"../../P125_21pi_vu_tf_float32/",
                                      "base_file_uvw":"P125_21pi_vu.$INDEX$.h5.uvw",
                                      "base_directory_vor":"../../P125_21pi_vu_tf_float32/",
                                      "base_file_vor":"P125_21pi_vu.$INDEX$.h5.uvw","datasets":[],"elem_spec":[],"shpx":1,
                                      "shpy":1,"shpz":1,"dx":1,"dy":1,"dz":1,"padding":15,"data_type":"float32",
                                      "data_folder":"Data","umean_file":"umean.txt","unorm_file":"norm.txt",
                                      "vornorm_file":"norm.txt","mean_norm":False,
                                      "delta_pred":1,"vormean_file":"vormean.txt"}.
        Data:
            - output_directory   : directory to store the tfrecord files
            - base_directory_uvw : directory of the files to store in the output directory
            - base_file_uvw      : name of the files stored in the base_directory_uvw
            - base_directory_vor : directory of the files to store in the output directory
            - base_file_vor      : name of the files stored in the base_directory_vor
            - datasets           : indices of the fields to read
            - elem_spec          : specification of the elements of the dataset
            - shpx               : shape in the streamwise direction
            - shpy               : shape in the wall-normal direction
            - shpz               : shape in the spanwise direction
            - dx                 : downsampling in the streamwise direction
            - dy                 : downsampling in the wall-normal direction
            - dz                 : downsampling in the spanwise direction
            - padding            : padding applied to the fields
            - data_type          : type of data used (float32,float16)
            - data_folder        : folder to store the data of the model
            - umean_file         : file storing the mean velocity
            - unorm_file         : file storing the normalization
            - vornorm_file       : file storing the normalization of the vorticity
            - mean_norm          : flag for choosing between standarization and normalization (True: standarization with
                                                                                               mean and standard deviation,
                                                                                               False: normalization 
                                                                                               between minimum and maximum)
            - delta_pred         : number of fields to advance in the predictions
            - vormean_file       : file of the mean vorticity

    Returns
    -------
    None.

    """
    # -------------------------------------------------------------------------------------------------------------------
    # Read data
    # -------------------------------------------------------------------------------------------------------------------
    output_directory   = str(data_in["output_directory"])
    base_directory_uvw = str(data_in["base_directory_uvw"])
    base_file_uvw      = str(data_in["base_file_uvw"])
    base_directory_vor = str(data_in["base_directory_vor"])
    base_file_vor      = str(data_in["base_file_vor"])
    datasets           = np.array(data_in["datasets"],dtype="int")
    elem_spec          = data_in["elem_spec"]
    shpx               = int(data_in["shpx"])
    shpy               = int(data_in["shpy"])
    shpz               = int(data_in["shpz"])
    dx                 = int(data_in["dx"])
    dy                 = int(data_in["dy"])
    dz                 = int(data_in["dz"])
    padding            = int(data_in["padding"])
    data_folder        = str(data_in["data_folder"])
    umean_file         = str(data_in["umean_file"])
    unorm_file         = str(data_in["unorm_file"])
    vornorm_file       = str(data_in["vornorm_file"])
    mean_norm          = bool(data_in["mean_norm"])
    delta_pred         = int(data_in["delta_pred"])
    data_type          = str(data_in["data_type"])
    vormean_file       = str(data_in["vormean_file"])
    num_datasets       = len(datasets)
    
    # -------------------------------------------------------------------------------------------------------------------
    # Create output directory if it does not exist
    # -------------------------------------------------------------------------------------------------------------------
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    
    # -------------------------------------------------------------------------------------------------------------------
    # Process each dataset file one by one and write the tfrecords
    # -------------------------------------------------------------------------------------------------------------------
    for idx in np.arange(len(datasets)):
        # ---------------------------------------------------------------------------------------------------------------
        # Create the data to load to the tfrecords
        # ---------------------------------------------------------------------------------------------------------------
        index       = int(datasets[idx])
        dataset_in  = {"folder_uvw":base_directory_uvw,"file_uvw":base_file_uvw,
                       "folder_vor":base_directory_vor,"file_vor":base_file_vor,
                       "elem_spec":elem_spec,"shpx":shpx,"shpy":shpy,
                       "shpz":shpz,"dx":dx,"dy":dy,"dz":dz,"padding":padding,"data_type":data_type,
                       "data_folder":data_folder,"umean_file":umean_file,"unorm_file":unorm_file,
                       "vornorm_file":vornorm_file,"mean_norm":mean_norm,
                       "index":index,"delta_pred":delta_pred,"vormean_file":vormean_file}
        dataset     = create_datasets(dataset_in)["data_XY"]
        # ---------------------------------------------------------------------------------------------------------------
        # Create the tfrecords
        # ---------------------------------------------------------------------------------------------------------------
        output_path = os.path.join(output_directory, f'dataset_{index}.tfrecord')
        data_write  = {"output_path":output_path,"dataset":dataset}
        write_tfrecord(data_in=data_write)
        print(f"Dataset {idx} written to {output_path}",flush=True)
    

def merge_data(data_in={"base_directory_uvw":"../../P125_21pi_vu_tf_float32/",
                        "base_file_uvw":"P125_21pi_vu.$INDEX$.h5.uvw",
                        "base_directory_vor":"../../P125_21pi_vu_tf_float32/",
                        "base_file_vor":"P125_21pi_vu.$INDEX$.h5.uvw",
                        "padding":15,"shpx":1,"shpy":1,"shpz":1,"dx":1,"dy":1,"dz":1,"data_type":"float32",
                        "datasets":[],"output_path":'/tfrecord/merged_vali_data/',"data_folder":"Data",
                        "umean_file":"umean.txt","unorm_file":"norm.txt","vornorm_file":"norm.txt",
                        "mean_norm":False,"delta_pred":1,"vormean_file":"vormean.txt"}):
    """
    .....................................................................................................................
    # merge_data: Function for reading the training and validation data in the tensorflow format for proper usage of 
                  Berzelius.
    .....................................................................................................................
    Parameters
    ----------
    data_in : dict, optional
        data required for selecting the fields for the training process.
        The default is {"base_directory_uvw":"../../P125_21pi_vu_tf_float32/",
                        "base_file_uvw":"P125_21pi_vu.$INDEX$.h5.uvw",
                        "base_directory_vor":"../../P125_21pi_vu_tf_float32/",
                        "base_file_vor":"P125_21pi_vu.$INDEX$.h5.uvw",
                        "padding":15,"shpx":1,"shpy":1,"shpz":1,"dx":1,"dy":1,"dz":1,"data_type":"float32",
                        "datasets":[],"output_path":'/tfrecord/merged_vali_data/',"data_folder":"Data",
                        "umean_file":"umean.txt","unorm_file":"norm.txt","vornorm_file":"norm.txt",
                        "mean_norm":False,"delta_pred":1,"vormean_file":"vormean.txt"}.
        Data:
            - base_directory_uvw : folder in which the data is stored
            - base_file_uvw      : file in which the data is stored
            - base_directory_vor : folder in which the data is stored
            - base_file_vor      : file in which the data is stored
            - padding            : padding of the fields
            - shpx               : shape in the streamwise direction
            - shpy               : shape in the wall-normal direction
            - shpz               : shape in the spanwise direction
            - dx                 : downsampling in the streamwise direction
            - dy                 : downsampling in the wall-normal direction
            - dz                 : downsampling in the spanwise direction
            - data_type          : type of data of the folder
            - datasets           : number of fields to read
            - output_path        : path to store the file
            - data_folder        : folder to store the data generated in the code
            - umean_file         : file containing the mean velocity
            - unorm_file         : file containing the normalization
            - mean_norm          : flag for defining the type of normalization (True: standarize with mean and standard
                                                                                deviation, False: normalize between 
                                                                                minimum and maximum)
            - delta_pred         : number of files to advance the predictions
            - vormean_file       : file of the mean vorticity
  
    Returns
    -------
    None
    """
    # -------------------------------------------------------------------------------------------------------------------
    # Read the data
    # -------------------------------------------------------------------------------------------------------------------
    base_directory_uvw = str(data_in["base_directory_uvw"])
    base_file_uvw      = str(data_in["base_file_uvw"])
    base_directory_vor = str(data_in["base_directory_vor"])
    base_file_vor      = str(data_in["base_file_vor"])
    padding            = int(data_in["padding"])
    shpx               = int(data_in["shpx"])
    shpy               = int(data_in["shpy"])
    shpz               = int(data_in["shpz"])
    dx                 = int(data_in["dx"])
    dy                 = int(data_in["dy"])
    dz                 = int(data_in["dz"])
    data_type          = str(data_in["data_type"])
    datasets           = np.array(data_in["datasets"],dtype='int')
    output_path        = str(data_in["output_path"])
    data_folder        = str(data_in["data_folder"])
    umean_file         = str(data_in["umean_file"])
    unorm_file         = str(data_in["unorm_file"])
    vornorm_file       = str(data_in["vornorm_file"])
    mean_norm          = bool(data_in["mean_norm"])
    delta_pred         = int(data_in["delta_pred"])
    vormean_file       = str(data_in["vormean_file"])
    if data_type == "float32":
        dtype = tf.float32
    elif data_type == "float16":
        dtype = tf.float16
    else:
        print("Exit the calculation due to invalid data type.",flush=True)
        sys.exit()
    
    # -------------------------------------------------------------------------------------------------------------------
    # Define your element spec here based on the data structure
    # -------------------------------------------------------------------------------------------------------------------
    elem_spec = (tf.TensorSpec(shape=(shpy,shpz+2*padding,shpx+2*padding,3),dtype=dtype),
                 tf.TensorSpec(shape=(shpy,shpz,shpx,3),dtype=dtype))
    
    # -------------------------------------------------------------------------------------------------------------------
    # Create the tfrecords
    # -------------------------------------------------------------------------------------------------------------------
    data_tfrecords = {"output_directory":output_path,"base_directory_uvw":base_directory_uvw,
                      "base_file_uvw":base_file_uvw,"base_directory_vor":base_directory_vor,
                      "base_file_vor":base_file_vor,
                      "datasets":datasets,"elem_spec":elem_spec,"shpx":shpx,"shpy":shpy,"shpz":shpz,"dx":dx,"dy":dy,
                      "dz":dz,"padding":padding,"data_type":data_type,"data_folder":data_folder,"umean_file":umean_file,
                      "unorm_file":unorm_file,"vornorm_file":vornorm_file,"mean_norm":mean_norm,"delta_pred":delta_pred,
                      "vormean_file":vormean_file}
    create_tfrecords(data_in=data_tfrecords)

