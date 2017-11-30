#	Title   : find_centroid_functions.py
#	Author  : Joe Park
#	Date    : 08/22/2017
#	Synopsis: Finds centroid of a point source from a pixel map image. DINO C-REx module.

import math         #common math functions
import numpy as np  #matrix algebra
import matplotlib.pyplot as plt
import cv2




##################################################
##################################################

# function to convert s/c state and beacon position to initial estimate of beacon centroid pixel/line location
# may not be necessary if pixel line location can be directly provided by navigation module

def initial_beacon_estimate(pos_beacon, pos_sc, attde_sc, cam_res, cam_fov):

    # TODO update function to convert nav module output into pixel/line estimate
    pixel = 10
    line = 10

    pl_beacon = (pixel, line)
    return pl_beacon



##################################################
##################################################

# Function to generate a ROI window around closest grouping of above threshold pixels near initial beacon estimate
# input:    pixel_map               original image map, [n x m] numpy array
#           pixel_line_beacon_i     tuple of initial centroid estimate, (row, col)
#           ROI_parameters          {signal_threshold, noise_threshold, ROI_size, ROI_border_width}
# output:   corner_ROI              lower left corner location of ROI
#           image_ROI               image map of ROI only [p x p] numpy array

def generate_point_source_ROI(pixel_map, pixel_line_beacon_i, ROI_parameters):

    signal_threshold = ROI_parameters['signal_threshold']
    noise_threhsold = ROI_parameters['noise_threshold']
    ROI_size = ROI_parameters['ROI_size']
    ROI_border_width = ROI_parameters['ROI_border_width']
    pixel_map_height, pixel_map_width = pixel_map.shape

    # calculate starting pixel-line value to initiate ROI
    pixel_line_beacon_i = (int(round(pixel_line_beacon_i[0])), int(round(pixel_line_beacon_i[1])))

    # find initial pixel location above signal threshold
    center_ROI_array = np.empty((0,2), dtype=int)

    if pixel_map[pixel_line_beacon_i[0], pixel_line_beacon_i[1]] >= signal_threshold:
        center_ROI_array = np.vstack(center_ROI_array, pixel_line_beacon_i[0], pixel_line_beacon_i[1])

    # if initial pixel is not above threshold begin expanding search outward to find a pixel above signal threshold
    else:

        # set limits of search perimeter
        search_row_min = pixel_line_beacon_i[0] -1
        search_row_max = pixel_line_beacon_i[0] +1
        search_col_min = pixel_line_beacon_i[1] -1
        search_col_max = pixel_line_beacon_i[1] +1

        signal_found = False

        # continue searching for a pixel above threshold or until entire image is searched
        while signal_found == False:

            # cycle through each row of search region
            for current_row in range(search_row_min, search_row_max+1):

                # cycle through horizontal edges of search border
                if current_row == search_row_min or current_row == search_row_max:

                    for current_col in range(search_col_min, search_col_max+1):
                        current_value = pixel_map[current_row, current_col]

                        if current_value >= signal_threshold:
                            center_ROI_array = np.vstack((center_ROI_array, np.array([current_row, current_col])))

                # cycle through vertical edges of search border
                else:
                    current_value_left = pixel_map[current_row, search_col_min]
                    current_value_right = pixel_map[current_row, search_col_max]

                    if current_value_left >= signal_threshold:
                        center_ROI_array = np.vstack((center_ROI_array, np.array([current_row, search_col_min])))
                    if current_value_right >= signal_threshold:
                        center_ROI_array = np.vstack((center_ROI_array, np.array([current_row, search_col_max])))

            # if no pixels above threshold found expand ROI
            if center_ROI_array.size > 0:
                signal_found = True
            else:
                search_row_min = search_row_min -1
                search_row_max = search_row_max +1
                search_col_min = search_col_min -1
                search_col_max = search_col_max +1

            # exit if entire image is searched without any signals found
            if search_row_min < 0 and search_col_min < 0 \
                and search_row_max > pixel_map_height-1 and search_col_max > pixel_map_width-1:
                signal_found = True
                print "\nNo signal found in image"

            # make sure search border is within image resolution
            if search_row_min < 0:
                search_row_min = 0
            if search_col_min < 0:
                search_col_min = 0
            if search_row_max > pixel_map_height-1:
                search_row_max = pixel_map_height-1
            if search_col_max > pixel_map_width-1:
                search_col_max = pixel_map_width-1

    # create a pixel_map of the region of interest only
    if center_ROI_array.size > 0:

        # TODO - add method of choosing center of ROI if multiple pixels above threshold are found in border
        # placeholder currently selects first entry

        center_ROI = np.empty([2], dtype=int)
        center_ROI[0] = center_ROI_array[0,0]
        center_ROI[1] = center_ROI_array[0,1]

        # TODO - add method of dynamically creating ROI around identified grouping of pixels above threshold
        # placeholder creates preset sized ROI centered around identified pixel

        ROI_min_row = center_ROI[0] - int(math.floor(ROI_size/2)) -1
        ROI_min_col = center_ROI[1] - int(math.floor(ROI_size/2)) -1
        ROI_max_row = center_ROI[0] + int(math.floor(ROI_size/2)) -1
        ROI_max_col = center_ROI[1] + int(math.floor(ROI_size/2)) -1

        # make sure ROI corner does not exceed image dimensions
        if ROI_max_row > pixel_map_width -1 - ROI_border_width:
            ROI_max_row = pixel_map_width -1 - ROI_border_width
        #if ROI_max_col > pixel_map_height -1 - ROI_border_width:
        #    ROI_max_col = pixel_map_height -1 - ROI_border_width
        if ROI_min_row < ROI_border_width:
            ROI_min_row = ROI_border_width
        #if ROI_min_col < ROI_border_width:
        #    ROI_min_col = ROI_border_width

        image_ROI = pixel_map[int(ROI_min_row):int(ROI_max_row+1), int(ROI_min_col):int(ROI_max_col+1)]
        print '\nROI center location and value'
        print center_ROI, pixel_map[int(center_ROI[0]), int(center_ROI[1])]
        print 'ROI limits'
        print (ROI_min_row, ROI_max_row), (ROI_min_col, ROI_max_col)

        corner_ROI = (ROI_min_row, ROI_min_col)

    else:
        image_ROI = pixel_map
        corner_ROI = (0,0)

    do_plots = False
    if do_plots == True:
        plt.figure(299)
        plt.imshow(image_ROI, extent=[0, ROI_size, 0, ROI_size])
        plt.colorbar()

    return (corner_ROI, image_ROI)



##################################################
##################################################

# Function to remove ROI border background noise value from ROI
# input:    pixel_map               pixel map for centroid calculations, [n x m] numpy array
#                                   (limited to ROI with border adjustment already applied)
#           ROI_parameters          {signal_threshold, noise_threshold, ROI_size, ROI_border_width}
# output:   pixel_map_corrected     pixel map with average value of ROI border subtracted from all values [n x m] array

def apply_ROI_border(pixel_map, ROI_parameters):

    border_ROI = ROI_parameters['ROI_border_width']
    num_rows, num_cols = pixel_map.shape
    num_pixels_border = border_ROI*2*num_rows + border_ROI*2*num_cols - 4*border_ROI**2

    ##############################################
    # find total value of all pixels in ROI border
    sum_pixels_border = 0

    for ind_col in range(1, num_cols+1):
        if ind_col <= border_ROI or ind_col > num_cols-border_ROI:
            for ind_row in range(1, border_ROI*2+1):
                sum_pixels_border = sum_pixels_border + pixel_map[ind_row-1, ind_col-1]
        else:
            for ind_row in range(1, border_ROI+1):
                sum_pixels_border = sum_pixels_border + pixel_map[ind_row - 1, ind_col - 1]
            for ind_row in range(num_rows - border_ROI+1, num_rows):
                sum_pixels_border = sum_pixels_border + pixel_map[ind_row - 1, ind_col - 1]

    ave_border = sum_pixels_border/num_pixels_border

    pixel_map_corrected = np.empty([num_rows, num_cols])
    for ind_col in range(0, num_cols):
        for ind_row in range(0, num_rows):
            pixel_map_corrected[ind_row, ind_col] = pixel_map[ind_row, ind_col] - ave_border

    return pixel_map_corrected



##################################################
##################################################

# Function to find the centroid of a pixel map using expected value
# input:    pixel_map               pixel map for centroid calculations, [n x m] numpy array
#                                   (limited to ROI with border adjustment already applied)
#           ROI_parameters          {signal_threshold, noise_threshold, ROI_size, ROI_border_width}
# output:   loc_centroid            location of the centroid in original image map coordinates (row, col)
#           DN                      total signal in ROI

def find_centroid(pixel_map, corner_ROI, ROI_parameters):

    border_ROI = ROI_parameters['ROI_border_width']
    num_rows, num_col = pixel_map.shape

    pixel_map_ROI = pixel_map[border_ROI:num_rows-border_ROI , border_ROI:num_col-border_ROI]
    num_rows_ROI = num_rows - border_ROI*2
    num_col_ROI = num_col - border_ROI*2

    # calculate total pixel value and centroid location
    DN = 0
    DN_row = 0
    DN_col = 0

    for ind_row in range(num_rows_ROI):
        for ind_col in range(num_col_ROI):
            DN = DN + pixel_map_ROI[ind_row, ind_col]
            DN_row = DN_row + pixel_map_ROI[ind_row, ind_col]*(ind_row+1)
            DN_col = DN_col + pixel_map_ROI[ind_row, ind_col]*(ind_col+1)

    loc_centroid_row = DN_col/DN
    loc_centroid_col = DN_row/DN
    loc_centroid = (loc_centroid_col+corner_ROI[0], loc_centroid_row+corner_ROI[1 ])

    do_plots = True
    if do_plots == True:
        fig298 = plt.figure(298)
        plt.imshow(pixel_map_ROI)
        plt.colorbar()
        plt.scatter(loc_centroid_row, loc_centroid_col, marker='x', s=100, c='w')
        fig298.suptitle('Search Region with Centroid Marked', fontsize=12, fontweight='bold')

    return loc_centroid, DN



##################################################
##################################################

def find_centroid_point_source(pixel_map, pixel_line_beacon_i, ROI_parameters):

    # crop origiinal image to an ROI based on initial
    corner_ROI, image_ROI = generate_point_source_ROI(pixel_map, pixel_line_beacon_i, ROI_parameters)

    # determine average value of region of interest border, subtract from rest of pixel map
    image_ROI = apply_ROI_border(image_ROI, ROI_parameters)

    # calculate centroid position and ROI brigtness value <-- centroid location is pixel number not index value (starts with 1)
    loc_centroid, DN = find_centroid(image_ROI, corner_ROI, ROI_parameters)

    return loc_centroid, DN


def hough_circles(img, blur=5, canny_thresh=200, dp=1, center_dist=200, accum=18, min_rad=0, max_rad=0, show_img=False):
    """This function attempts to find the center of curves in a given image
    If the minimum and/or maximum circle radius is unknown leave as 0
    @param img The image is the only required argument
    @param blur The kernel size
    @param canny_thresh The threshold for the hysteresis procedure in the canny function
    @param dp Inverse ratio of the accumulator resolution to the image resolution, recommend leaving as 1
    @param center_dist The minimum distance between centers of detected circles
    @param accum The accumulator threshold for circle centers. Smaller value gives more error.
    @param min_rad The minimum radius for the circles
    @oaram max_rad The maximum radius for the circles
    @param show_img Set to show the images at each stage in the processing
    @return An array of circles that should have the best matches first if there are multiple
    """
    if show_img:
        cv2.imshow('before', img)
        cv2.waitKey(0)
        cv2.destroyAllWindows()

    orrig_img = cv2.cvtColor(img, cv2.COLOR_GRAY2BGR)
    img = cv2.medianBlur(img, blur)
    img = cv2.GaussianBlur(img, (blur, blur), 0)

    if show_img:
        cv2.imshow('blur', img)
        cv2.waitKey(0)
        cv2.destroyAllWindows()

    if show_img:
        canny_img = cv2.Canny(img, canny_thresh, canny_thresh / 15)

        cv2.imshow('canny', canny_img)
        cv2.waitKey(0)
        cv2.destroyAllWindows()

    circles = cv2.HoughCircles(img, cv2.HOUGH_GRADIENT, dp, center_dist, param1=canny_thresh, param2=accum, minRadius=min_rad, maxRadius=max_rad)

    if show_img:
        tmp = np.uint16(np.around(circles))

        average = np.uint16(np.mean(tmp[0, :], axis=0))

        for i in circles[0, :]:
            cv2.circle(orrig_img, (i[0], i[1]), i[2], (0, 255, 0), 2)
            cv2.circle(orrig_img, (i[0], i[1]), 2, (0, 0, 255), 3)

        cv2.imshow('Detected circles', orrig_img)
        cv2.waitKey(0)
        cv2.destroyAllWindows()

    return circles
