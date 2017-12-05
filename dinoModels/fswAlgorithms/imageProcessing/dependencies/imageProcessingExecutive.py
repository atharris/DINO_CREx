import math
import matplotlib.pyplot as plt
import numpy as np
# from scipy import misc

import image_processing_functions as imfunc
import search_location_functions as locfunc
import object_id_functions as idfunc
import dynamics as dyn


##################################################
##################################################
# Inputs to Image Processing Module

# beacon list IDs and radius (defined at module initialization)
beaconIDs = ('Earth', 'Moon')
radius_beacons = (6378.1, 1737.)
numBeacons = len(beaconIDs)

# read in inputs for image processing module
# file_in = np.load('CDR_save_files/90_deg_orig.npz')
file_in = np.load('CDR_save_files/stars_only_cdr.npz')

# image map input from image generation module
ex_image = file_in['detector_array']
ex_image = (ex_image / max(ex_image)) * 2**12
ex_image = ex_image.reshape(512, 512)

timestampPic = 1.E9

# camera parameters from DRM definition
cam_res = (512, 512)
cam_pixel_size = (39E-6, 39E-6)  # horizontal, vertical [m]
cam_focal_length = .05  # [m]
cam_sensor_size = (cam_res[0] * cam_pixel_size[0], cam_res[1] * cam_pixel_size[1])  # [m]
cam_fov = (2 * math.degrees(math.atan2(cam_sensor_size[0] / 2., cam_focal_length)),
           2 * math.degrees(math.atan2(cam_sensor_size[1] / 2., cam_focal_length)))

# most recent beacon position from nav module
# N_r_beacons = (file_in['earth_pos'], file_in['moon_pos'])
# use below to put beacons out of field of view (if only interested in image processing stars)
N_r_beacons = (np.array([-1E3, -1E3, 0]), np.array([-1E3, -1E3, -1E3]))

# spacecraft state from nav module
N_r_cam = file_in['sc_pos']
BN_dcm_cam = file_in['sc_dcm']

print '\nImage Processing Module Inputs:'
print 'Spacecraft Position: ', N_r_cam
print 'Spacecraft Attitude:'
print BN_dcm_cam
print 'Beacon IDs, Radius, and Positions: '
for ind in range(numBeacons):
    print beaconIDs[ind], radius_beacons[ind], N_r_beacons[ind]

print '\nCamera Parameters'
print 'Camera Resolution: ', cam_res
print 'Camera Pixel Size [m]: ', cam_pixel_size
print 'Camera Focal Length [m]: ', cam_focal_length
print 'Camera Sensor Size [m]: ', cam_sensor_size
print 'Camera Field of View [deg]: ', round(cam_fov[0],4), round(cam_fov[1],4)


##################################################
##################################################
# Parameter definition for image processing module

cameraParam = {}
cameraParam['resolution'] = cam_res
cameraParam['focal length'] = cam_focal_length
cameraParam['sensor size'] = cam_sensor_size
cameraParam['field of view'] = cam_fov
cameraParam['pixel size'] = cam_pixel_size

ROI_parameters = {}
ROI_parameters['signal_threshold'] = 1.5
ROI_parameters['noise_threshold'] = 1E-6
ROI_parameters['ROI_size'] = 100
ROI_parameters['ROI_border_width'] = 1
ROI_parameters['max_search_dist'] = 50

imageProcessingParam ={}
imageProcessingParam['voteCountMinRatio'] = .5      # minimum ratio of positive matches out of possible matches
imageProcessingParam['dthetaMax'] = 15.     #[deg] dependent on object ID reference catalog
imageProcessingParam['filenameSearchCatalog'] = 'star_catalog/tycho_BTmag_cutoff.db'
imageProcessingParam['filenameObjectIDCatalog'] = 'star_catalog/objectID_catalog.db'
imageProcessingParam['dthetaError'] = 5E-6

pixelSizeMin = 15.      # number of pixels required for a beacon to be considered a resolved body
pixelDimenMin = min(cam_pixel_size)
singlePixelAngSize = math.degrees(math.atan(pixelDimenMin/cam_focal_length))
angularSizeMin = singlePixelAngSize * pixelSizeMin

print 'Beacon Angular Size Requirements: ', angularSizeMin


##################################################
##################################################
# Check Current status of target beacons
# 0 for out of field of view, 1 for point source, 2 for resolved body

beaconStatus = imfunc.checkBeaconStatus(N_r_beacons, N_r_cam, BN_dcm_cam, cameraParam['field of view'],
                                        radius_beacons, angularSizeMin)

print '\nBeacon Status: ', beaconStatus

# Separate visible beacons from non-visible
N_r_beaconsVisible = []
N_r_beaconsPtSource = []
N_r_beaconsResolved = []
beaconIDsVisible = []
beaconIDsPtSource = []
beaconIDsResolved = []

for ind in range(numBeacons):

    if beaconStatus[ind] != 0:
        N_r_beaconsVisible.append(N_r_beacons[ind])
        beaconIDsVisible.append(beaconIDs[ind])

    if beaconStatus[ind] == 1:
        N_r_beaconsPtSource.append(N_r_beacons[ind])
        beaconIDsPtSource.append(beaconIDs[ind])

    if beaconStatus[ind] == 2:
        N_r_beaconsResolved.append(N_r_beacons[ind])
        beaconIDsResolved.append(beaconIDs[ind])

numBeaconsVisible = len(N_r_beaconsVisible)
numBeaconsPtSource = len(N_r_beaconsPtSource)
numBeaconsResolved = len(N_r_beaconsResolved)

print '\nVisible Beacons: ', beaconIDsVisible
print 'Pt Source Beacons: ', beaconIDsPtSource
print 'Resolved Beacons: ', beaconIDsResolved


##################################################
##################################################
# Generate Pixel Line Estimates for Stars and Beacons

pixel_line_ptSource_i, tempID = locfunc.initial_stars_estimate(
    BN_dcm_cam, cameraParam, imageProcessingParam['filenameSearchCatalog'])

if numBeaconsResolved > 0:
    pixel_line_resolved_i = locfunc.initial_beacons_estimate(
        N_r_beaconsResolved, N_r_cam, BN_dcm_cam, cameraParam)

if numBeaconsPtSource > 0:
    pixel_line_beacon_ptSource_i = locfunc.initial_beacons_estimate(
        N_r_beaconsPtSource, N_r_cam, BN_dcm_cam, cameraParam)
    pixel_line_ptSource_i = np.vstack((pixel_line_ptSource_i, pixel_line_beacon_ptSource_i))


print '\n# of Initial Estimates of Pt Sources: ', len(pixel_line_ptSource_i)
# for currentPL in pixel_line_ptSource_i:
#     print currentPL

if numBeaconsResolved > 0:
    print '\nInitial Estimate of Resolved Bodies'

    for currentPL in pixel_line_resolved_i:
        print currentPL

##################################################
##################################################
# Determine Center and Centroid Locations

if numBeaconsResolved > 0:
    pixel_line_center, DN = imfunc.find_center_resolved_body(ex_image, pixel_line_resolved_i, ROI_parameters)
else:
    pixel_line_center, DN = imfunc.find_centroid_point_source(ex_image, pixel_line_ptSource_i, ROI_parameters)

numObjects = len(pixel_line_center)

print '\nIdentified Pixel and Line Center Locations: '
for ind in range(numObjects):
    print pixel_line_center[ind]


##################################################
##################################################
# Run Object ID Algorithm
#
# if numBeaconsResolved == 0:
#     objectID = idfunc.objectID_stars(pixel_line_center, BN_dcm_cam, imageProcessingParam, cameraParam)
#
#
#     print '\nObject ID Results:'
#     for indStar in range(numBeacons):
#         print objectID[indStar]
#
# N_radec = idfunc.catalogIDsToRaDec(objectID, imageProcessingParam['filenameSearchCatalog'])
# #
# print '\nInertial right ascension and declination for identified objects (reference catalog lookup):'
# for indStar in range(numBeacons): print objectID[indStar], '\t', N_radec[indStar]
#
# # List containers for unit vectors (1x3) numpy arrays
# B_ehat = []
# N_ehat = []
#
# print '\nBody and Inertial Unit Vectors'
# for indStar in range(numBeacons):
#
#     B_ehatTemp =idfunc.pixelline_to_ehat((pixel_star[indStar], line_star[indStar]), camParam)
#     B_ehat.append(B_ehatTemp)
#
#     N_ehatTemp = dyn.radec_to_unitvector(N_radec[indStar])
#     N_ehat.append(N_ehatTemp)
#
#     # verify results with actual s/c attitude
#     B_ehatCheck = np.matmul(BN_dcm_cam, N_ehatTemp)
#
#     print B_ehatTemp, '\t', N_ehatTemp, '\t', B_ehatCheck - B_ehatTemp
#
#
# # DCM = dyn.quest(B_ehat, N_ehat, (2,1,1,1))
# dcmQUEST = dyn.quest(B_ehat, N_ehat)
# mrpQUEST = dyn.dcm2mrp(dcmQUEST)
#
# print '\nQUEST Attitude Solution: '
# print dcmQUEST
#
# print '\nS/C Attitude Truth: '
# print BN_dcm_cam
#
# print '\nAttitude Error:'
# print dcmQUEST-BN_dcm_cam
#
# print '\nMRP Solution: ', mrpQUEST















