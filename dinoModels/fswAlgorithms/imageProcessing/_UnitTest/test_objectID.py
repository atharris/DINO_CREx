# Title   : find_centroid_functions.py
# Author  :

import sys, os, inspect
sys.path.append('../dependencies/')

import math         # common math functions
import numpy as np  # matrix algebra
import searchLocationFunctions as locfunc
import objectIDFunctions as idfunc
import dynamicFunctions as dyn

#############################################################
#############################################################

ROI_parameters = {}
ROI_parameters['signal_threshold'] = 1.5
ROI_parameters['noise_threshold'] = 1E-6
ROI_parameters['ROI_size'] = 100
ROI_parameters['ROI_border_width'] = 2
ROI_parameters['max_search_dist'] = 50

imageProcessingParam ={}
imageProcessingParam['voteCountMinRatio'] = .5
imageProcessingParam['dthetaMax'] = 15
imageProcessingParam['filenameSearchCatalog'] = '../../../../external/tycho_mag_cutoff.db'
imageProcessingParam['filenameObjectIDCatalog'] = '../../../../external/objectID_catalog.db'
imageProcessingParam['dthetaError'] = 5E-6

fname_catalog = '../../../../external/tycho_mag_cutoff.db'

#############################################################
#############################################################

do_CDR_stars = True
if do_CDR_stars:
    file_in = np.load('test_cases/stars_only_cdr.npz')

    pos_sc = file_in['sc_pos']
    attde_sc = file_in['sc_dcm']
    # attde_sc = np.matmul(attde_sc, dyn.eulerDCM_313(math.radians(1), math.radians(1), 0))

    ex_image = file_in['detector_array']
    ex_image = (ex_image/np.amax(ex_image)) * 255
    ex_image = ex_image.reshape(512, 512)

    cam_res = (512, 512)
    cam_pixel_size = (39E-6, 39E-6)         # horizontal, vertical [m]
    cam_focal_length = .05                  # [m]
    cam_sensor_size = (cam_res[0] * cam_pixel_size[0], cam_res[1] * cam_pixel_size[1])  # [m]
    cam_fov = ( 2 * math.degrees(math.atan2(cam_sensor_size[0]/2., cam_focal_length)),
                2 * math.degrees(math.atan2(cam_sensor_size[1]/2., cam_focal_length)))

    camParam = {}
    camParam['focal_length'] = cam_focal_length
    camParam['resolution'] = cam_res
    camParam['pixel_size'] = (39E-6, 39E-6)

    # fig1 = plt.figure(1)
    # plt.imshow(ex_image)
    # fig1.suptitle('Original Image', fontsize=12, fontweight='bold')
    # plt.show()

    cameraParam = {}
    cameraParam['resolution'] = cam_res
    cameraParam['focal length'] = cam_focal_length
    cameraParam['sensor size'] = cam_sensor_size
    cameraParam['field of view'] = cam_fov
    cameraParam['pixel size'] = cam_pixel_size

# generate pixel line estimates for stars in camera field of view
pixelLineEstimate, catalogID = locfunc.initial_stars_estimate(
    attde_sc, cameraParam, fname_catalog)

# num_obj = 5
num_obj = len(pixelLineEstimate)

# add random noise to pixel and line coordinates
doNoise = False
if doNoise:
    sd = 1E-5
    noise = np.random.normal(0, sd, num_obj)
    tempPixel = []
    tempLine = []
    for ind in range(num_obj):
        tempPixel.append(pixelLineEstimate[ind]+noise[ind])
        tempLine.append(pixelLineEstimate[ind]+noise[ind])
    pixel_star = tempPixel
    line_star = tempLine
    print '\nNoise added (Gaussian with simga=', sd, ')'
    print noise

print 'Number of Initial Estimates: ', num_obj

dthetaMax = 15.
voteMin = round(num_obj/2.)

# Object ID
objectID = idfunc.objectIDStars(pixelLineEstimate[0:num_obj],
                                 imageProcessingParam,
                                 cameraParam)

print '\nObject ID Results: '

print '\nStar #', '\t\t', \
    "%6s" % 'Pixel ', '\t\t', \
    "%6s" % 'Line  ', '\t\t', \
    "%6s" % 'Ref ID', '\t\t', \
    "%6s" % 'Obj ID', '\t\t'

for ind in range(num_obj):
    print 'Star #', ind, '\t:', \
        "%6s" % round(pixelLineEstimate[ind][0], 2), '\t', \
        "%6s" % round(pixelLineEstimate[ind][1], 2), '\t\t', \
        "%6s" % catalogID[ind], '\t\t', \
        "%6s" % objectID[ind]

N_radec = idfunc.catalogIDsToRaDec(objectID, imageProcessingParam['filenameSearchCatalog'])
#
print '\nInertial right ascension and declination for identified objects (reference catalog lookup):'
for indStar in range(num_obj): print objectID[indStar], '\t', N_radec[indStar]














