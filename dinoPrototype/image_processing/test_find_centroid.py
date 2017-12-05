#	Title   : test_find_centroid.py
#	Author  : Joe Park
#	Date    : 08/22/2017
#	Synopsis: Tests centroid finding functions as used in DINO C-REx.

import numpy as np
from scipy import misc
import math
import matplotlib.pyplot as plt
import image_processing_functions as imfunc
import search_location_functions as locfunc


##################################################
##################################################

# signal_threshold, noise_threshold, ROI_size (n x n pixel border), single side ROI_border_width
ROI_parameters = {}
ROI_parameters['signal_threshold'] = 1.5
ROI_parameters['noise_threshold'] = 1E-6
ROI_parameters['ROI_size'] = 100
ROI_parameters['ROI_border_width'] = 1
ROI_parameters['max_search_dist'] = 50

cam_res = (512, 512)
cam_pixel_size = (39E-6, 39E-6)  # horizontal, vertical [m]
cam_focal_length = .05  # [m]
cam_sensor_size = (cam_res[0] * cam_pixel_size[0], cam_res[1] * cam_pixel_size[1])  # [m]
cam_fov = (2 * math.degrees(math.atan2(cam_sensor_size[0] / 2., cam_focal_length)),
           2 * math.degrees(math.atan2(cam_sensor_size[1] / 2., cam_focal_length)))

cameraParam = {}
cameraParam['resolution'] = cam_res
cameraParam['focal length'] = cam_focal_length
cameraParam['sensor size'] = cam_sensor_size
cameraParam['field of view'] = cam_fov
cameraParam['pixel size'] = cam_pixel_size

# reference star catalog filename
fname_catalog = 'star_catalog/tycho_BTmag_cutoff.db'

# select which scenario to run
do_CDR_stars = False
do_CDR_beacon = True

doplot_centroid = True
doplot_isearch = True

##################################################
##################################################

# note: image map is in the (line, pixel) format ... same as (row, column) .... (y, x)

if do_CDR_beacon:
    file_in = np.load('CDR_save_files/90_deg.npz')
    pos_beacon = np.vstack((file_in['earth_pos'], file_in['moon_pos']))
    pos_sc = file_in['sc_pos']
    attde_sc = file_in['sc_dcm']
    # attde_sc = np.matmul(attde_sc, dyn.eulerDCM_313(math.radians(1), math.radians(1), 0))
    print attde_sc
    ex_image = file_in['detector_array']
    ex_image = (ex_image / max(ex_image)) * 512
    ex_image = ex_image.reshape(512, 512)
    plt.imshow(ex_image)
    plt.show()

    # generate pixel line estimates for beacons in camera field of view
    pixel_line_beacon_i = locfunc.initial_beacons_estimate(
        pos_beacon, pos_sc, attde_sc, cameraParam)

    # pass in three initial estimates to test ROI generation logic
    pixel_truth = np.array([256, 394.99206801, 240])
    line_truth = np.array([256, 256, 245])
    ROI_estimates = []
    for ind in range(len(pixel_truth)):
        ROI_estimates.append((pixel_truth[ind], line_truth[ind]))

    pixel_line_center, DN = imfunc.find_center_resolved_body(ex_image, ROI_estimates, ROI_parameters)


if do_CDR_stars:
    file_in = np.load('CDR_save_files/stars_only_cdr.npz')

    pos_sc = file_in['sc_pos']
    attde_sc = file_in['sc_dcm']
    # attde_sc = np.matmul(attde_sc, dyn.eulerDCM_313(math.radians(1), math.radians(1), 0))

    ex_image = file_in['detector_array']
    ex_image = (ex_image/np.amax(ex_image)) * 255
    ex_image = ex_image.reshape(512, 512)

    # generate pixel line estimates for stars in camera field of view
    # pixel_truth, line_star, star_catalog = locfunc.initial_stars_estimate(
    #     attde_sc, cameraParam, fname_catalog)

    # find centroid
    pixel_truth = np.array([443.71657484, 493.96318093, 95.3500384])
    line_truth = np.array([103.55964507, 105.8875924, 152.30559197])
    ROI_estimates = []
    for ind in range(len(pixel_truth)):
        ROI_estimates.append((pixel_truth[ind], line_truth[ind]))
    pixel_line_center, DN = imfunc.find_centroid_point_source(ex_image, ROI_estimates, ROI_parameters)

    print pixel_line_center

##################################################
##################################################

# Print all parameters

print '\nCamera Parameters:'
print 'resolution: ', cam_res
print 'focal length [mm]: ', cam_focal_length
print 'pixel size [mm]: ', cam_pixel_size
print 'sensor size [mm]: ', cam_sensor_size
print 'field of view [deg]: ', cam_fov
print '\nInitial Image Shape'
print ex_image.shape
print 'S/C Attitude (BN DCM): \n', attde_sc

if do_CDR_beacon:

    print '\nS/C and Beacon Parameters'
    print 'S/C Position (HCI [km]): ', pos_sc
    print 'Beacon Position (HCI [km]): '

    for current_beacon in pos_beacon:
        print current_beacon

    print '\n\nInitial Beacon Location Estimate'
    print 'Beacon Location (pixel/line coord.): ', pixel_line_beacon_i

if do_CDR_stars:

    print '\n\nInitial Star Location Estimate:'
    print pixel_truth
    n_star = len(pixel_truth)
    print 'stars in field of view: ', n_star
    print 'min/max pixel coord.: ', min(pixel_truth), max(pixel_truth)
    print 'min/max line coord.: ', min(line_truth), max(line_truth)
    print '\nPixel Line Centroid Locations:'
    print pixel_line_center
    print pixel_line_center.shape


##################################################
##################################################

if doplot_isearch == True:
    fig1 = plt.figure(1)
    # plt.imshow(ex_image, cmap='gray')
    # plt.scatter(loc_centroid[0], loc_centroid[1], color='r', marker='x', s=50)
    plt.imshow(ex_image)
    # plt.colorbar()

    if do_CDR_beacon:

        plt.scatter(pixel_line_beacon_i[0][0], pixel_line_beacon_i[0][1], color='r', marker='+', s=30)
        plt.scatter(pixel_line_beacon_i[1][0], pixel_line_beacon_i[1][1], color='r', marker='+', s=30)

        #plt.savefig('CDR_save_files/90_deg_orig_initial_estimate.png')

    if do_CDR_stars:

        # for ind in range(pixel_line_truth_i[0].shape[1]):
        for ind in range(3):
            plt.scatter(pixel_truth[ind], line_truth[ind], color='r', marker='+', s=30)
            fig1.suptitle('Original Image with Initial Beacon Location Estimate',
                          fontsize=12, fontweight='bold')
            #plt.savefig('CDR_save_files/stars_only_initial_estimate.png')

    fig1.suptitle('Original Image with Initial Beacon Location Estimate', fontsize=12, fontweight='bold')
    plt.show()


if doplot_centroid == True and do_CDR_stars == True:

    for ind_cent in range(n_star):

        if pixel_line_center is not None:

            crop_center = (int(round(pixel_line_center[ind_cent][0]/1, 1)),
                               int(round(pixel_line_center[ind_cent][1]/1, 1)))
            print '\nSubplot Center Coordinate'
            print crop_center
            crop_box = 5
            plt.imshow(ex_image[crop_center[1]-crop_box:crop_center[1]+crop_box,
                      crop_center[0]-crop_box:crop_center[0]+crop_box],
                      interpolation='none', cmap='viridis')

            # plot measured centroid location
            plt.scatter(pixel_line_center[ind_cent][0]-crop_center[0]+crop_box,
                       pixel_line_center[ind_cent][1]-crop_center[1]+crop_box,
                            color='r', marker='x', s=75)

            # plot truth value
            plt.scatter(pixel_truth[ind_cent]-crop_center[0]+crop_box-.5,
                            line_truth[ind_cent]-crop_center[1]+crop_box-.5,
                            color='b', marker='+', s=75)
            # plt.scatter(0, 0, color='w', marker='+', s=75)
            plt.show()


if doplot_centroid == True and do_CDR_beacon == True:

    for ind_cent in range(len(pixel_line_center)):

        if pixel_line_center[ind_cent] is not None:

            crop_center = (int(round(pixel_line_center[ind_cent][0], 1)), \
                           int(round(pixel_line_center[ind_cent][1], 1)))

            print '\nSubplot Center Coordinate'
            print crop_center
            crop_box = 30

            plt.imshow(ex_image[crop_center[1]-crop_box:crop_center[1]+crop_box,
                       crop_center[0]-crop_box:crop_center[0]+crop_box],
                       interpolation='none', cmap='viridis')

            # plot measured center
            plt.scatter(pixel_line_center[ind_cent][0]-crop_center[0]+crop_box,
                        pixel_line_center[ind_cent][1]-crop_center[1]+crop_box,
                        color='r', marker='x', s=75)

            # plot truth value
            plt.scatter(pixel_truth[ind_cent]-crop_center[0]+crop_box-.5,
                        line_truth[ind_cent]-crop_center[1]+crop_box-.5,
                        color='b', marker='o', s=75)
            # plt.scatter(1, 1, color='w', marker='o', s=75)

            plt.show()
