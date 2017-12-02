# Title   : find_centroid_functions.py
# Author  :

import os.path
import math         # common math functions
import numpy as np  # matrix algebra
import search_location_functions as locfunc
import object_id_functions as idfunc
import dynamics as dyn

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
imageProcessingParam['filenameSearchCatalog'] = 'star_catalog/tycho_BTmag_cutoff.db'
imageProcessingParam['filenameObjectIDCatalog'] = 'star_catalog/objectID_catalog.db'
imageProcessingParam['dthetaError'] = 5E-6


fname_catalog = 'star_catalog/tycho_BTmag_cutoff.db'

#############################################################
#############################################################

do_CDR_stars = True
if do_CDR_stars:
    file_in = np.load('CDR_save_files/stars_only_cdr.npz')

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
pixel_star, line_star, star_catalog = locfunc.initial_stars_estimate(
    attde_sc, cameraParam, fname_catalog)

# num_obj = 5
num_obj = len(pixel_star)

# add random noise to pixel and line coordinates
doNoise = True
if doNoise:
    sd = 1E-5
    noise = np.random.normal(0, sd, num_obj)
    tempPixel = []
    tempLine = []
    for ind in range(num_obj):
        tempPixel.append(pixel_star[ind]+noise[ind])
        tempLine.append(line_star[ind]+noise[ind])
    pixel_star = tempPixel
    line_star = tempLine
    print '\nNoise added (Gaussian with simga=', sd, ')'
    print noise

print 'Number of Initial Estimates: ', len(pixel_star)
print 'Pixel, Line, and Catalog ID'
for ind in range(num_obj):
    print 'Star #', ind, '\t:', \
        "%6s" % round(pixel_star[ind], 2), '\t', \
        "%6s" % round(line_star[ind], 2), '\t\t', \
        "%6s" % star_catalog[ind]

dthetaMax = 15.
voteMin = round(num_obj/2)

# Object ID
objectID = idfunc.objectID_stars((pixel_star[0:num_obj], line_star[0:num_obj]),
                                                    attde_sc, imageProcessingParam,
                                                    cameraParam)

print '\nObject ID Results:'
for indStar in range(num_obj):
    print objectID[indStar]

N_radec = idfunc.catalogIDsToRaDec(objectID, imageProcessingParam['filenameSearchCatalog'])
#
print '\nInertial right ascension and declination for identified objects (reference catalog lookup):'
for indStar in range(num_obj): print objectID[indStar], '\t', N_radec[indStar]

# List containers for unit vectors (1x3) numpy arrays
B_ehat = []
N_ehat = []

print '\nBody and Inertial Unit Vectors'
for indStar in range(num_obj):

    B_ehatTemp =idfunc.pixelline_to_ehat((pixel_star[indStar], line_star[indStar]), camParam)
    B_ehat.append(B_ehatTemp)

    N_ehatTemp = dyn.radec_to_unitvector(N_radec[indStar])
    N_ehat.append(N_ehatTemp)

    # verify results with actual s/c attitude
    B_ehatCheck = np.matmul(attde_sc, N_ehatTemp)

    print B_ehatTemp, '\t', N_ehatTemp, '\t', B_ehatCheck - B_ehatTemp


# Compute s/c attitude from unit vectors using QUEST
# Test case using solution from ASEN5010 HW3
# v1_b = np.array([.8273, .5541, -.0920])
# v2_b = np.array([-.8285, .5522, -.0955])
# v3_b = np.array([.2155, .5522, .8022])
# v4_b = np.array([.5570, -.7442, -.2884])
# v1_n = np.array([-.1517, -.9669, .2050])
# v2_n = np.array([-.8393, .4494, -.3044])
# v3_n = np.array([-.0886, -.5856, -.80])
# v4_n = np.array([.8814, -.0303, .5202])
# B_ehat = (v1_b, v2_b, v3_b, v4_b)
# N_ehat = (v1_n, v2_n, v3_n, v4_n)


# DCM = dyn.quest(B_ehat, N_ehat, (2,1,1,1))
dcmQUEST = dyn.quest(B_ehat, N_ehat)
mrpQUEST = dyn.dcm2mrp(dcmQUEST)

print '\nQUEST Attitude Solution: '
print dcmQUEST

print '\nS/C Attitude Truth: '
print attde_sc

print '\nAttitude Error:'
print dcmQUEST-attde_sc

print '\nMRP Solution: ', mrpQUEST















