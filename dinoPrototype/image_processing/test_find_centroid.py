#	Title   : test_find_centroid.py
#	Author  : Joe Park
#	Date    : 08/22/2017
#	Synopsis: Tests centroid finding functions as used in DINO C-REx.

import numpy as np
from scipy import misc
import matplotlib.pyplot as plt
import image_processing_functions as imfunc




##################################################
##################################################

# inputs to measurement generation module

# camera parameters
cam_res = (1024, 1024)      # horizontal, vertical [pixels]
cam_fov = (30, 30)          # horizontal, vertical [deg]

# Nav module inputs
pos_beacon = np.array([1E3, 1E3, 0])
pos_sc = np.array([1E3, 0, 0])
attde_sc = np.array([[0,1,0],
                     [-1,0,0],
                     [0,0,1]])


##################################################
##################################################

# load example image to use
# note: image map is in the (line, pixel) format ... same as (row, column) .... (y, x)

#ex_image = misc.imread('psf_examples/psf_lower_right2.png', mode='L')
ex_image = misc.imread('psf_examples/psf_example_landscape.png', mode='L')


# TODO calculate/select ROI parameters <-- unclear if handled by image gen already
# signal_threshold, noise_threshold, ROI_size (n x n pixel border), single side ROI_border_width
ROI_parameters = {}
ROI_parameters['signal_threshold'] = 250
ROI_parameters['noise_threshold'] = 10
ROI_parameters['ROI_size'] = 250
ROI_parameters['ROI_border_width'] = 2

# convert output from Nav module to estimate of target beacon pixel/line coordinates
# Note: pixel_line_beacon_i is in (line, pixel) format
pixel_line_beacon_i = imfunc.initial_beacon_estimate(pos_beacon, pos_sc, attde_sc, cam_res, cam_fov)
# placeholder below places initial estimate in the center of the image
ex_res = ex_image.shape
pixel_line_beacon_i = (int(round(ex_res[0]*0.5)), int(round(ex_res[1]*.5)))
pixel_line_beacon_i = (400, 800)

print '\nInitial Image Resolution and Initial Estimate Location'
print ex_res, pixel_line_beacon_i

# find centroid
pixel_line_centroid, DN = imfunc.find_centroid_point_source(ex_image, pixel_line_beacon_i, ROI_parameters)


#################################################
# display results

print 'Initial Image Shape'
print ex_image.shape
print '\nCentroid Location: ','%.2f' % pixel_line_centroid[0], '%.2f' % pixel_line_centroid[1]
print 'DN Value: ', DN

doplot = True
if doplot == True:
    fig1 = plt.figure(1)
    # plt.imshow(ex_image, cmap='gray')
    # plt.scatter(loc_centroid[0], loc_centroid[1], color='r', marker='x', s=50)
    plt.imshow(ex_image)
    plt.colorbar()
    plt.scatter(pixel_line_centroid[1], pixel_line_centroid[0], color='w', marker='x', s=50)
    plt.savefig('saved_output/psf_centroid_landscape.png')
    fig1.suptitle('Original Image with Centroid Marked', fontsize=12, fontweight='bold')
    plt.show()

#plt.figure(2)
#plt.imshow(ex_rast)
#plt.show()