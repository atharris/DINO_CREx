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
# import object_id_functions as idfunc
import dynamics as dyn


##################################################
##################################################

# inputs to measurement generation module

# camera parameters
cam_res = (2000, 1000)            # horizontal, vertical [pixels]
cam_focal_length = 100.         # [mm]
cam_pixel_size = (.2, .2)       # horizontal, vertical [mm]

cam_sensor_size = (cam_res[0]*cam_pixel_size[0], cam_res[1]*cam_pixel_size[1])  # [mm]
cam_fov = (2 * math.degrees(math.atan2(cam_sensor_size[0]/2., cam_focal_length)),
           2 * math.degrees(math.atan2(cam_sensor_size[1]/2., cam_focal_length)))

# Nav module inputs
pos_beacon = np.array([[10100, 10000, 0]])
pos_sc = np.array([10000, 0, 0])

# Euler 321 DCM (90 deg, 0, 22.5)
attde_sc = dyn.eulerDCM_321(math.pi/2, 0, math.pi/8)

# Euler 313 DCM (90 deg, 90 deg, 45 deg) (midway polar example)
#attde_sc = dyn.eulerDCM_313(math.pi/2, math.pi/2, math.pi/8)

# 90 Deg Rotation about 3rd axis
#attde_sc = dyn.eulerDCM_313(math.pi/2, 0, 0)

# BCT NSC1 example image (may be -65.2 for final rotation)
attde_sc = dyn.eulerDCM_321(math.radians(251.115), math.radians(-32.842), math.radians(65.2))


##################################################
##################################################

# load example image to use
# note: image map is in the (line, pixel) format ... same as (row, column) .... (y, x)

# load example image from image generation module
do_CDR_beacon = True
if do_CDR_beacon:
    ROI_parameters = {}
    ROI_parameters['signal_threshold'] = .01
    ROI_parameters['noise_threshold'] = .001
    ROI_parameters['ROI_size'] = 100
    ROI_parameters['ROI_border_width'] = 2
    ROI_parameters['max_search_dist'] = 50

    ROI_estimates = ((256, 256),)# (270, 238))# (390, 256))#, (240, 245))

    file_in = np.load('CDR_save_files/0_deg_orig.npz')
    pos_beacon = np.vstack((file_in['earth_pos'], file_in['moon_pos']))
    pos_sc = file_in['sc_pos']
    attde_sc = file_in['sc_dcm']
    # attde_sc = np.matmul(attde_sc, dyn.eulerDCM_313(math.radians(1), math.radians(1), 0))
    print attde_sc
    ex_image = file_in['detector_array']
    ex_image = (ex_image / max(ex_image)) * 256
    ex_image = ex_image.reshape(512, 512)

    for i in range (379, 399):
        for j in range (245, 267):
            if (ex_image[j][i] >= .01):
                ex_image[j-15][i-116] = ex_image[j][i] + ex_image[j-15][i-116]

    # plt.imshow(ex_image)
    # plt.show()
    cam_res = (512, 512)
    cam_pixel_size = (39E-6, 39E-6)         # horizontal, vertical [m]
    cam_focal_length = .05                  # [m]

    cam_sensor_size = (cam_res[0] * cam_pixel_size[0], cam_res[1] * cam_pixel_size[1])  # [m]
    cam_fov = ( 2 * math.degrees(math.atan2(cam_sensor_size[0]/2., cam_focal_length)),
                2 * math.degrees(math.atan2(cam_sensor_size[1]/2., cam_focal_length)))

do_CDR_stars = False
if do_CDR_stars:
    ROI_parameters = {}
    ROI_parameters['signal_threshold'] = 1.5
    ROI_parameters['noise_threshold'] = 1E-6
    ROI_parameters['ROI_size'] = 100
    ROI_parameters['ROI_border_width'] = 1
    ROI_parameters['max_search_dist'] = 50


    #ROI_estimates = [(282, 122), (330, 112), (23, 29)]
    ROI_estimates = [(1, 1), (510, 509)]

    file_in = np.load('CDR_save_files/stars_only.npz')

    pos_sc = file_in['sc_pos']
    attde_sc = file_in['sc_dcm']
    # attde_sc = np.matmul(attde_sc, dyn.eulerDCM_313(math.radians(1), math.radians(1), 0))

    ex_image = file_in['detector_array']
    ex_image = (ex_image/np.amax(ex_image)) * 255
    ex_image = ex_image.reshape(512,512)

    for i in range(280, 285):
        for j in range(121, 126):
            if (ex_image[j][i] > 5):
                ex_image[j - 121][i - 280] = ex_image[j][i]
            elif (ex_image[j][i] >= 1.5):
                ex_image[j - 121][i - 280] = 5

    for i in range(318, 322):
        for j in range(122, 126):
            ex_image[j + 386][i + 190] = ex_image[j][i] * 5

    cam_res = (512, 512)
    cam_pixel_size = (39E-6, 39E-6)         # horizontal, vertical [m]
    cam_focal_length = .05                  # [m]
    cam_sensor_size = (cam_res[0] * cam_pixel_size[0], cam_res[1] * cam_pixel_size[1])  # [m]
    cam_fov = ( 2 * math.degrees(math.atan2(cam_sensor_size[0]/2., cam_focal_length)),
                2 * math.degrees(math.atan2(cam_sensor_size[1]/2., cam_focal_length)))
    # fig1 = plt.figure(1)
    # plt.imshow(ex_image)
    # plt.colorbar()
    # fig1.suptitle('Original Image', fontsize=12, fontweight='bold')
    # plt.show()
    # plt.savefig('CDR_save_files/stars_only.png')


do_CDR_stars2 = False
if do_CDR_stars2:
    ROI_parameters = {}
    ROI_parameters['signal_threshold'] = 1
    ROI_parameters['noise_threshold'] = .001
    ROI_parameters['ROI_size'] = 100
    ROI_parameters['ROI_border_width'] = 2
    ROI_parameters['max_search_dist'] = 50

    ROI_estimates = [(282, 122), (330, 112), (23, 29)]

    ex_image = misc.imread('psf_examples/bct_nsc1.png')
    print ex_image.shape

    attde_sc = dyn.eulerDCM_321(math.radians(251.115), math.radians(-32.842), math.radians(-65.2))

    cam_res = (868, 725)
    cam_fov = (14.4, 12)
    cam_pixel_size = (1E-6, 1E-6)         # horizontal, vertical [m]
    cam_focal_length = (cam_res[0] * cam_pixel_size[0]) / math.tan(math.radians(cam_fov[0]))

    cam_sensor_size = (2 * cam_focal_length * math.tan(cam_fov[0]/2.),
                       2 * cam_focal_length * math.tan(cam_fov[1]/2.))
    cam_fov = ( 2 * math.degrees(math.atan2(cam_sensor_size[0]/2., cam_focal_length)),
                2 * math.degrees(math.atan2(cam_sensor_size[1]/2., cam_focal_length)))

    print cam_fov

# star catalog filename
fname_catalog = 'star_catalog/tycho_mag7cutoff.db'

# TODO calculate/select ROI parameters <-- unclear if handled by image gen already
# signal_threshold, noise_threshold, ROI_size (n x n pixel border), single side ROI_border_width
#ROI_parameters = {}
#ROI_parameters['signal_threshold'] = 1
#ROI_parameters['noise_threshold'] = .001
#ROI_parameters['ROI_size'] = 100
#ROI_parameters['ROI_border_width'] = 2
#ROI_parameters['max_search_dist'] = 50


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

print '\nS/C and Beacon Parameters'
print 'S/C Position (HCI [km]): ', pos_sc
print 'Beacon Position (HCI [km]): '
for current_beacon in pos_beacon:
    print current_beacon
print 'S/C Attitude (BN DCM): \n', attde_sc


##################################################
##################################################

# Run Image Processing Module

# generate pixel line estimates for beacons in camera field of view
pixel_line_beacon_i = locfunc.initial_beacons_estimate(
    pos_beacon, pos_sc, attde_sc, cam_res, cam_focal_length, cam_pixel_size)

print '\nInitial Beacon Location Estimate'
print 'Beacon Location (pixel/line coord.): ',pixel_line_beacon_i

# placeholder below places initial estimate in the center of the image
#ex_res = ex_image.shape
#pixel_line_beacon_i = (int(round(ex_res[0]*0.5)), int(round(ex_res[1]*.5)))
#pixel_line_beacon_i = (400, 800)

# geenrate pixel line estimates for stars in camera field of view
pixel_star, line_star, star_catalog = locfunc.initial_stars_estimate(
    attde_sc, cam_res, cam_focal_length, cam_pixel_size, fname_catalog)

print '\nInitial Star Location Estimate:'
n_star = pixel_star.shape[1]
print 'stars in field of view: ', n_star
print 'min/max pixel coord.: ', np.amin(pixel_star, 1), np.amax(pixel_star, 1)
print 'min/max line coord.: ', np.amin(line_star, 1), np.amax(line_star, 1)
print '3 Brightest Stars: '
for ind in range(4):
    print pixel_star[0,ind], line_star[0,ind], star_catalog[0][ind], star_catalog[1][ind], star_catalog[2][ind]

# pixel_line_star_i = (5,10)
# pixel_line_estimate_i = np.vstack((pixel_line_beacon_i, pixel_line_star_i))


# PLACEHOLDER ---- ONLY DO ONE BEACON AT A TIME
# pixel_line_beacon_i[0]

# find centroid
# pixel_line_centroid, DN = imfunc.find_centroid_point_source(ex_image, ((101, 453),(282, 122), (330, 112)), ROI_parameters, 3)
# pixel_line_centroid, DN = imfunc.find_centroid_point_source(ex_image, ( (pixel_star[0, 0], line_star[0, 0]),
#                                                                         (pixel_star[0, 1], line_star[0, 1]),
#                                                                         (pixel_star[0, 2], line_star[0, 2])),
#                                                                         ROI_parameters, 3)

# PLACEHOLDER ---- Included (240, 245) to check for finding same beacon twice
pixel_line_centroid, DN = imfunc.find_center_resolved_body(ex_image, ROI_estimates, ROI_parameters)

print '\nPixel Line Centroid Locations:'
# pixel_line_centroid = (0, 0)
# DN = 0
print pixel_line_centroid
print pixel_line_centroid.shape

##################################################
##################################################

# Display Results

#print '\nCentroid Location: ','%.2f' % pixel_line_centroid[0], '%.2f' % pixel_line_centroid[1]
#print 'DN Value: ', DN

doplot_isearch = True
if doplot_isearch == True:
    fig1 = plt.figure(1)
    # plt.imshow(ex_image, cmap='gray')
    # plt.scatter(loc_centroid[0], loc_centroid[1], color='r', marker='x', s=50)
    plt.imshow(ex_image)
    # plt.colorbar()

    if do_CDR_beacon:
        plt.scatter(pixel_line_beacon_i[0][0], pixel_line_beacon_i[1][0], color='r', marker='+', s=30)
        plt.scatter(390, 256, color='r', marker='+', s=30)

        #plt.savefig('SER_output/earth.png')

    if do_CDR_stars or do_CDR_stars2:

        # for ind in range(pixel_line_star_i[0].shape[1]):
        for ind in range(20):
            plt.scatter(pixel_star[0, ind], line_star[0, ind], color='r', marker='+', s=30)
            fig1.suptitle('Original Image with Initial Beacon Location Estimate', fontsize=12, fontweight='bold')
            #plt.savefig('CDR_save_files/stars_only_initial_estimate.png')

    fig1.suptitle('Original Image with Initial Beacon Location Estimate', fontsize=12, fontweight='bold')
    plt.show()

    # roi1 = ex_image[450:455, 98:105]
    # print '\nROI 1 Post-Processed Plot Shape:'
    # print roi1.shape
    # plt.imshow(roi1, interpolation='none')
    # # plt.scatter(pixel_star[0, 0], line_star[0, 1], color='r', marker='+', s=30)
    # plt.scatter(2.7, 1.8, color='r', marker='+', s=30)
    # plt.show()


doplot_centroid = False
if doplot_centroid == True:
    fig1 = plt.figure()
    # plt.imshow(ex_image, cmap='gray')
    # plt.scatter(loc_centroid[0], loc_centroid[1], color='r', marker='x', s=50)
    plt.imshow(ex_image)

    for ind_cent in range(len(pixel_line_centroid)):
        plt.scatter(pixel_line_centroid[ind_cent][0], pixel_line_centroid[ind_cent][1], color='r', marker='x', s=75)

    plt.show()
    fig1.suptitle('Original Image with Centroid Marked', fontsize=12, fontweight='bold')

#plt.figure(2)
#plt.imshow(ex_rast)
#plt.show()

# circ = imfunc.hough_circles(ex_image, center_dist=50, blur=1, accum=5, show_img=False)
# print '\nFull Image Hough Transform Results: ', circ
