# -*- coding: utf-8 -*-
"""
Created on Wed Mar 19 07:49:39 2025

@author: andre
"""

from py_bin.py_functions.read_tfrecord import read_tfrecord_test
from py_bin.py_functions.read_norm_velocity import read_norm_velocity
import numpy as np

index = 10000
pad   = 15
shpx  = 384
shpy  = 201
shpz  = 288
deltapred = 25

data = read_tfrecord_test(data_in={"tfrecord_folder":'../../tfrecord/',"interval":[index],
                                   "test_size":0.2,"padding":pad,"shpx":shpx,
                                   "shpy":shpy,"shpz":shpz,"data_type":"float32"})["data"]
numpy_data = [(x, y) for x, y in data.as_numpy_iterator()]

tfr_in = numpy_data[0][0]
tfr_out = numpy_data[0][1]


veli = read_norm_velocity(data_in={"folder":"../../P125_83pi","file":"P125_83pi.$INDEX$.h5.uvw",
                                   "padding":pad,"shpx":shpx,"shpy":shpy,"shpz":shpz,"dx":1,"dy":1,"dz":1,
                                   "data_folder":"../../P125_83pi_data/d20240703_Data/","umean_file":"Umean.txt",
                                   "unorm_file":"norm.txt","index":index,"data_type":"float32",
                                   "mean_norm":False})["norm_velocity"]

velo_ini = np.zeros((201,288+30,384+30,3))
velo_ini[:,:,:,0] = veli["unorm"]
velo_ini[:,:,:,1] = veli["vnorm"]
velo_ini[:,:,:,2] = veli["wnorm"]


velo = read_norm_velocity(data_in={"folder":"../../P125_83pi","file":"P125_83pi.$INDEX$.h5.uvw",
                                       "padding":0,"shpx":shpx,"shpy":shpy,"shpz":shpz,"dx":1,"dy":1,"dz":1,
                                       "data_folder":"../../P125_83pi_data/d20240703_Data/","umean_file":"Umean.txt",
                                       "unorm_file":"norm.txt","index":index+deltapred,"data_type":"float32",
                                       "mean_norm":False})["norm_velocity"]

velo_out = np.zeros((201,288,384,3))
velo_out[:,:,:,0] = velo["unorm"]
velo_out[:,:,:,1] = velo["vnorm"]
velo_out[:,:,:,2] = velo["wnorm"]


print(f"error in {np.mean(abs(velo_ini-tfr_in))}")
print(f"error out {np.mean(abs(velo_out-tfr_out))}")