import math
import matplotlib.pyplot as plt
import numpy as np
import sys, os, inspect

sys.path.append('../dependencies/')

import imageProcessingFunctions as imfunc
import searchLocationFunctions as locfunc
import objectIDFunctions as idfunc
import dynamicFunctions as dyn
import imageProcessingExecutive as ip


##################################################
##################################################
# Inputs to Image Processing Module
#
# cam_res = (512, 512)
# cam_pixel_size = (39E-6, 39E-6)  # horizontal, vertical [m]
# cam_focal_length = .05  # [m]
# cam_sensor_size = (cam_res[0] * cam_pixel_size[0], cam_res[1] * cam_pixel_size[1])  # [m]
# cam_fov = (2 * math.degrees(math.atan2(cam_sensor_size[0] / 2., cam_focal_length)),
#            2 * math.degrees(math.atan2(cam_sensor_size[1] / 2., cam_focal_length)))

file_in = np.load('test_cases/testcase0_2.npz')

ex_image = file_in['imageMap']

print file_in['imageMap']

ex_image = ex_image.reshape(512, 512)

BN_dcm_cam = file_in['scDCM']

r_N_cam = file_in['scPos'][0:3]

# beacon list IDs and radius (defined at module initialization)
beaconIDs = ('Earth', 'Moon')
beaconRadius = (6378.1, 1737.)
r_N_beacons = (file_in['earthPos'][0:3],
                file_in['moonPos'][0:3])


cameraParam = {}
cameraParam['resolution'] = file_in['resolution']
cameraParam['focal length'] = file_in['focalLength']
cameraParam['sensor size'] = file_in['sensorSize']
cameraParam['field of view'] = file_in['FoV']
cameraParam['pixel size'] = file_in['pixelSize']


##################################################
##################################################
# Run image processing executive function

# only include detected beacons and PL coordinates
idsOutput, pixelLineOutput, sigma_BN_output = ip.imageProcessing(
    ex_image, cameraParam, r_N_cam, BN_dcm_cam,
    r_N_beacons, beaconIDs, beaconRadius, makePlots=True, debugMode=True)


















