#	Title   : find_centroid_functions.py
#	Author  : Joe Park
#	Date    : 08/22/2017
#	Synopsis: Finds centroid of a point source from a pixel map image. DINO C-REx module.

import math         #common math functions
import numpy as np  #matrix algebra
import matplotlib.pyplot as plt

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
#           pixel_line_beacon_i     Array of tuples of initial centroid estimate, (row, col)
#           ROI_parameters          {signal_threshold, noise_threshold, ROI_size, ROI_border_width}
#           num_beacons             Number of beacons to find
# output:   corner_ROI              Array of lower left corner locations of ROI
#           image_ROI               Array of image maps of ROI only [p x p] numpy array

def generatePointSourceROI(pixel_map, pixel_line_beacon_i, ROI_parameters):

    # Pull out threshold parameters
    signal_threshold = ROI_parameters['signal_threshold']
    noise_threshold = ROI_parameters['noise_threshold']
    ROI_size = ROI_parameters['ROI_size']
    ROI_border_width = ROI_parameters['ROI_border_width']
    pixel_map_height, pixel_map_width = pixel_map.shape

    # Initialize return parameters
    corner_ROI = []
    image_ROI = []

    # Find the corner ROI and image ROI for each beacon
    for i in pixel_line_beacon_i:
        # calculate starting pixel-line value to initiate ROI
        pixel_line_beacon_i = (int(round(i[0])), int(round(i[1])))

        # find initial pixel location above signal threshold
        center_ROI_array = np.empty((0, 2), dtype=int)

        # If the initial pixel is above the threshold, use it.
        if pixel_map[pixel_line_beacon_i[0], pixel_line_beacon_i[1]] >= signal_threshold:
            center_ROI_array = np.vstack(center_ROI_array, pixel_line_beacon_i[0], pixel_line_beacon_i[1])

        # if initial pixel is not above threshold begin expanding search outward to find a pixel above signal threshold
        else:

            # set limits of search perimeter
            search_row_min = pixel_line_beacon_i[0] - 1
            search_row_max = pixel_line_beacon_i[0] + 1
            search_col_min = pixel_line_beacon_i[1] - 1
            search_col_max = pixel_line_beacon_i[1] + 1

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
                    search_row_min = search_row_min - 1
                    search_row_max = search_row_max + 1
                    search_col_min = search_col_min - 1
                    search_col_max = search_col_max + 1

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
            center_ROI[0] = center_ROI_array[0, 0]
            center_ROI[1] = center_ROI_array[0, 1]

            # TODO - add method of dynamically creating ROI around identified grouping of pixels above threshold
            # placeholder creates preset sized ROI centered around identified pixel

            ROI_min_row = center_ROI[0] - int(math.floor(ROI_size/2)) - 1
            ROI_min_col = center_ROI[1] - int(math.floor(ROI_size/2)) - 1
            ROI_max_row = center_ROI[0] + int(math.floor(ROI_size/2)) - 1
            ROI_max_col = center_ROI[1] + int(math.floor(ROI_size/2)) - 1

            # make sure ROI corner does not exceed image dimensions
            if ROI_max_row > pixel_map_width - 1 - ROI_border_width:
                ROI_max_row = pixel_map_width - 1 - ROI_border_width
            #if ROI_max_col > pixel_map_height - 1 - ROI_border_width:
            #    ROI_max_col = pixel_map_height - 1 - ROI_border_width
            if ROI_min_row < ROI_border_width:
                ROI_min_row = ROI_border_width
            #if ROI_min_col < ROI_border_width:
            #    ROI_min_col = ROI_border_width

            image_ROI.append(pixel_map[int(ROI_min_row):int(ROI_max_row+1), int(ROI_min_col):int(ROI_max_col+1)])
            print '\nROI center location and value'
            print center_ROI, pixel_map[int(center_ROI[0]), int(center_ROI[1])]
            print 'ROI limits'
            print (ROI_min_row, ROI_max_row), (ROI_min_col, ROI_max_col)

            corner_ROI.append([ROI_min_row, ROI_min_col])

        else:
            image_ROI.append(pixel_map)
            corner_ROI.append([0, 0])

        do_plots = False

        if do_plots == True:
            plt.figure(299)
            plt.imshow(image_ROI[i], extent=[0, ROI_size, 0, ROI_size])
            plt.colorbar()

    return (corner_ROI, image_ROI)

##################################################
##################################################

# Function to remove ROI border background noise value from ROI
# input:    pixel_map               pixel map for centroid calculations, [n x m] numpy array
#                                   (limited to ROI with border adjustment already applied)
#           ROI_parameters          {signal_threshold, noise_threshold, ROI_size, ROI_border_width}
# output:   pixel_map_corrected     pixel map with average value of ROI border subtracted from all values [n x m] array

def applyROIBborder(pixel_map, ROI_parameters):

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

def find_centroid_point_source(pixel_map, pixel_line_beacon_i, ROI_parameters, num_beacons):

    loc_centroid = np.empty([num_beacons], dtype=tuple)
    DN = np.empty([num_beacons], dtype=int)

    # crop original image to an ROI based on initial
    corner_ROI, image_ROI = generatePointSourceROI(pixel_map, pixel_line_beacon_i, ROI_parameters, num_beacons)

    for i in range(0, num_beacons):
        # determine average value of region of interest border, subtract from rest of pixel map
        image_ROI[i] = applyROIBborder(image_ROI[i], ROI_parameters)

        # calculate centroid position and ROI brigtness value <-- centroid location is pixel number not index value (starts with 1)
        loc_centroid[i], DN[i] = find_centroid(image_ROI[i], corner_ROI[i], ROI_parameters)

    return loc_centroid, DN

##################################################
##################################################

#def neighbors_are_on(pixel_loc):
    # NDY
    #return False

##################################################
##################################################

#def find_upper_edge(initial_pixel_loc):
    # NDY
    #return (0,0)

##################################################
##################################################

#def find_lower_edge(initial_pixel_loc):
    # NDY
    #return (0,0)

##################################################
##################################################

#def find_right_edge(initial_pixel_loc):
    # NDY
    #return (0,0)

##################################################
##################################################

#def find_left_edge(initial_pixel_loc):
    # NDY
    #return (0,0)