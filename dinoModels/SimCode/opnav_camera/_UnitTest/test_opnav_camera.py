#!/usr/bin/python

''' '''
'''
 ISC License

 Copyright (c) 2016-2017, Autonomous Vehicle Systems Lab, University of Colorado at Boulder

 Permission to use, copy, modify, and/or distribute this software for any
 purpose with or without fee is hereby granted, provided that the above
 copyright notice and this permission notice appear in all copies.

 THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

'''
#
#   Integrated Unit Test Script
#   Purpose:  Run a test of the star tracker module
#   Author:  John Alcorn
#   Creation Date:  October 12, 2016
#

import pytest
import sys, os, inspect
import numpy as np
import ctypes
import math
import csv
import logging

filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))
splitPath = path.split('Basilisk')
sys.path.append(splitPath[0]+'/Basilisk/modules')
sys.path.append(splitPath[0]+'/Basilisk/PythonModules')
sys.path.append(splitPath[0]+'/Basilisk/ADCSAlgorithms/sensorInterfaces/STSensorData')

# import MessagingAccess
# import SimulationBaseClass
# import unitTestSupport  # general support file with common unit test functions
# import matplotlib.pyplot as plt
# import macros
# import star_tracker
# import sim_model
# import RigidBodyKinematics as rbk
# import spice_interface
import opnav_camera
from datetime import datetime
import pdb


def camera_test(camera,name,datafile):
	start_time = datetime.now()
	deltas = []
	print(name + " Star Field Initialize Camera: " + str(datetime.now()-start_time))
	start_time = datetime.now()

	camera.updateState()
	print(name + " Star Field First Update: " + str(datetime.now()-start_time))
	start_time = datetime.now()

	for i in range(0,10):
	    camera.updateState()
	    delta = datetime.now()-start_time
	    print("Update " + str(i) + ": " + str(delta))
	    start_time = datetime.now()
	    deltas.append((delta).microseconds)


	print(name + " Star Field Min: " + str(min(deltas)/1000000.) + "s")
	print(name + " Star Field Max: " + str(max(deltas)/1000000.) + "s")

	print(name + " Star Field Mean: " + str(np.mean(deltas)/1000000.) + "s")
	print(name + " Star Field STD: " + str(np.std(deltas)/1000000.) + "s")
	Verified_Data = np.load(datafile)

	print(max(abs(Verified_Data['pixel'][0:2000]-camera.OpnavMessage.pixel_out)))
	print(max(abs(Verified_Data['line'][0:2000]-camera.OpnavMessage.line_out)))
	print(max(camera.OpnavMessage.mag_out))
camera = opnav_camera.OpnavCamera()
camera.alpha = np.deg2rad(85);
camera.beta = np.deg2rad(0);
camera.gamma = np.deg2rad(72);
camera.f = 3;
camera.a = 1.5;
camera.b = 2.5;
camera.alpha_resolution = 512;
camera.beta_resolution = 1024;
camera_test(camera,"Orion","Orion_Data.npz")

camera = opnav_camera.OpnavCamera()
camera.alpha = np.deg2rad(232.5);
camera.beta = np.deg2rad(-80);
camera.gamma = np.deg2rad(30);
camera.f = 5;
camera.a = 2.5;
camera.b = 2.5;
camera.alpha_resolution = 1024;
camera.beta_resolution = 1024;
camera_test(camera,'UMi',"UMi_Data.npz")

camera = opnav_camera.OpnavCamera()
camera.alpha = np.deg2rad(187);
camera.beta = np.deg2rad(59);
camera.gamma = np.deg2rad(0);
camera.f = 4;
camera.a = 2.3;
camera.b = 2.3;
camera.alpha_resolution = 512;
camera.beta_resolution = 1024;
camera_test(camera,'Crux',"Crux_Data.npz")


#Manual Verification
Crux_Data = np.load("Crux_Data.npz")
ind = np.extract(Crux_Data['VTmag'], Crux_Data['VTmag'] < 4)
Crux_RA_list = Crux_Data['RA'][ind][[18,21,17,20,8,11]]
Crux_DE_list = Crux_Data['DE'][ind][[18,21,17,20,8,11]]
Crux_VTmag_list = Crux_Data['VTmag'][ind][[18,21,17,20,8,11]]

Orion_Data = np.load("Orion_Data.npz")
ind = np.extract(Orion_Data['VTmag'], Orion_Data['VTmag'] < 3.57)
Ori_RA_list = 
Ori_DE_list = DE[[7,6,9,10,12,8,4,1,0,3]]
Ori_VTmag_list = 
!$
UMi_Data = np.load("UMi_Data.npz")
ind = np.extract(UMi_Data['VTmag'], UMi_Data['VTmag'] < 5)
UMi_RA_list = 
UMi_DE_list = DE[[14,2,1,12,8,10,17,0,9,11]]
UMi_VTmag_list = 
pdb.set_trace()








