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

# Function to generate a ROI window around closest grouping of above threshold pixels near initial beacon estimate
# input:    pixel_map               original image map, [n x m] numpy array
#           pixel_line_beacon_i     Array of tuples of initial centroid estimate, (row, col)
#           ROI_parameters          {signal_threshold, noise_threshold, ROI_size, ROI_border_width}
#           num_beacons             Number of beacons to find
# output:   corner_ROI              Array of lower left corner locations of ROI
#           image_ROI               Array of image maps of ROI only [p x p] numpy array

def generate_point_source_ROI(pixel_map, pixel_line_beacon_i, ROI_parameters):

    # Pull out threshold parameters
    signal_threshold = ROI_parameters['signal_threshold']
    noise_threshold = ROI_parameters['noise_threshold']
    ROI_size = ROI_parameters['ROI_size']
    ROI_border_width = ROI_parameters['ROI_border_width']
    max_search_pixels = ROI_parameters['max_search_dist']
    pixel_map_height, pixel_map_width = pixel_map.shape

    # Initialize return parameters
    corner_ROI = []
    image_ROI = []

    # Find the corner ROI and image ROI for each beacon
    for i in range(0, len(pixel_line_beacon_i)):
        # calculate starting pixel-line value to initiate ROI
        current_beacon = (int(round(pixel_line_beacon_i[i][0])), int(pixel_line_beacon_i[i][1]))

        # find initial pixel location above signal threshold
        center_ROI_array = np.empty((0, 2), dtype=int)

        # If the initial pixel is above the threshold, use it.
        if pixel_map[current_beacon[1], current_beacon[0]] >= signal_threshold \
                and not already_found_ROI((current_beacon[1], current_beacon[0]), corner_ROI, image_ROI):
            center_ROI_array = np.vstack((center_ROI_array, current_beacon))
            signal_found = True

        # if initial pixel is not above threshold begin expanding search outward to find a pixel above signal threshold
        else:

            # set limits of search perimeter
            search_row_min = current_beacon[1] - 1
            search_row_max = current_beacon[1] + 1
            search_col_min = current_beacon[0] - 1
            search_col_max = current_beacon[0] + 1

            signal_found = False
            search_iterations = 0

            # continue searching for a pixel above threshold or until max number of searches reached
            while (not signal_found) and (search_iterations < max_search_pixels):

                # cycle through each row of search region
                for current_row in range(search_row_min, search_row_max+1):

                    # cycle through horizontal edges of search border
                    if current_row == search_row_min or current_row == search_row_max:

                        for current_col in range(search_col_min, search_col_max+1):
                            current_value = pixel_map[current_row, current_col]

                            if current_value >= signal_threshold \
                                    and not already_found_ROI((current_row, current_col), corner_ROI, image_ROI):
                                center_ROI_array = np.vstack((center_ROI_array, np.array([current_col, current_row])))
                                signal_found = True

                    # cycle through vertical edges of search border
                    else:
                        current_value_left = pixel_map[current_row, search_col_min]
                        current_value_right = pixel_map[current_row, search_col_max]

                        if current_value_left >= signal_threshold \
                                    and not already_found_ROI((current_row, search_col_min), corner_ROI, image_ROI):
                            center_ROI_array = np.vstack((center_ROI_array, np.array([search_col_min, current_row])))
                            signal_found = True

                        elif current_value_right >= signal_threshold \
                                    and not already_found_ROI((current_row, search_col_max), corner_ROI, image_ROI):
                            center_ROI_array = np.vstack((center_ROI_array, np.array([search_col_max, current_row])))
                            signal_found = True

                # if no pixels above threshold found expand ROI
                #if center_ROI_array.size > 0:
                #    signal_found = True
                #else:
                if not signal_found:
                    search_row_min = search_row_min - 1
                    search_row_max = search_row_max + 1
                    search_col_min = search_col_min - 1
                    search_col_max = search_col_max + 1

                # make sure search border is within image resolution
                if search_row_min < 0:
                    search_row_min = 0
                if search_col_min < 0:
                    search_col_min = 0
                if search_row_max > pixel_map_height-1:
                    search_row_max = pixel_map_height-1
                if search_col_max > pixel_map_width-1:
                    search_col_max = pixel_map_width-1

                # Made another search.  Increment search_iterations.
                search_iterations = search_iterations + 1

            # If signal was not found within allowed search iterations, print error message.
            if not signal_found:
                print "\nNo signal found in image"

        # create a pixel_map of the region of interest only if we found a signal
        if signal_found:

            center_ROI = np.empty([2], dtype=int)
            center_ROI[0] = center_ROI_array[0, 0]
            center_ROI[1] = center_ROI_array[0, 1]

            ROI_min_row = find_lowest_pixel(center_ROI_array[0], pixel_map, signal_threshold)
            ROI_min_col = find_left_pixel(center_ROI_array[0], pixel_map, signal_threshold)
            ROI_max_row = find_highest_pixel(center_ROI_array[0], pixel_map, signal_threshold)
            ROI_max_col = find_right_pixel(center_ROI_array[0], pixel_map, signal_threshold)

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

            # print '\nROI center location and value'
            # print center_ROI, pixel_map[int(center_ROI[1]), int(center_ROI[0])]
            # print 'ROI limits'
            # print (ROI_min_row, ROI_max_row), (ROI_min_col, ROI_max_col)

            corner_ROI.append([ROI_min_row, ROI_min_col])

        do_plots = False
        if do_plots == True:
            #plt.figure(299)
            plt.imshow(image_ROI[i], extent=[0, ROI_size, 0, ROI_size], interpolation='none', cmap='viridis')
            plt.show()

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

    # print '\nCurrent ROI Pixel Map'
    # print pixel_map

    ##############################################
    # find total value of all pixels in ROI border

    if pixel_map.shape[0] >= 2*border_ROI or pixel_map.shape[1] >= 2*border_ROI:

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
                newvalue = pixel_map[ind_row, ind_col] - ave_border
                if newvalue >= 0:
                    pixel_map_corrected[ind_row, ind_col] = newvalue
                else:
                    pixel_map_corrected[ind_row, ind_col] = 0

    else:
        pixel_map_corrected = pixel_map

    return pixel_map_corrected

##################################################
##################################################

# Function to find the centroid of a pixel map using expected value
# input:    pixel_map               pixel map for centroid calculations, [n x m] numpy array
#                                   (limited to ROI with border adjustment already applied)
#           ROI_parameters          {signal_threshold, noise_threshold, ROI_size, ROI_border_width}
# output:   loc_centroid            location of the centroid in original image map coordinates (x, y)
#           DN                      total signal in ROI

def find_centroid(pixel_map, corner_ROI, ROI_parameters):

    border_ROI = ROI_parameters['ROI_border_width']
    num_rows, num_col = pixel_map.shape

    #pixel_map_ROI = pixel_map[border_ROI:num_rows-border_ROI, border_ROI:num_col-border_ROI]
    pixel_map_ROI = pixel_map

    num_rows_ROI = num_rows
    num_col_ROI = num_col
    #num_rows_ROI = num_rows - border_ROI*2
    #num_col_ROI = num_col - border_ROI*2

    # calculate total pixel value and centroid location
    DN = 0
    DN_row = 0
    DN_col = 0

    for ind_row in range(num_rows_ROI):
        for ind_col in range(num_col_ROI):
            DN = DN + pixel_map_ROI[ind_row, ind_col]
            DN_row = DN_row + pixel_map_ROI[ind_row, ind_col] * (ind_row)
            DN_col = DN_col + pixel_map_ROI[ind_row, ind_col] * (ind_col)

    loc_centroid_row = DN_col/DN
    loc_centroid_col = DN_row/DN
    loc_centroid = (loc_centroid_row + corner_ROI[1], loc_centroid_col + corner_ROI[0])

    do_plots = False
    if do_plots == True:
        # print '\nCurrent pixel_map_ROI Size: ', pixel_map_ROI.shape
        # print pixel_map_ROI
        # print loc_centroid_row, loc_centroid_col
        plt.imshow(pixel_map_ROI, interpolation='none', cmap='viridis')
        # plt.scatter(loc_centroid_row, loc_centroid_col, marker='x', s=150, linewidth=2, c='r')
        plt.show()
        #plt.suptitle('Search Region with Centroid Marked', fontsize=12, fontweight='bold')

    return loc_centroid, DN


##################################################
##################################################

# Function to determine if enough neighbors are on for a pixel to be considered a hit
# input:    pixel_loc               The location of the pixel to check
#           pixel_map               Map of pixels to examine
#           threshold               Minimum value to be considered a hit
# output:   enough_neighbors        True if there are enough neighbors on for the pixel to be considered a hit and false
#                                   otherwise
def neighbors_are_on(pixel_loc, pixel_map, threshold):

    # Initialize neighbors_on to 0 and pull off the height and width of hte pixel map.
    # Also set the minimum number of neighbors needed to be on to consider it a hit.  Number includes pixel_loc.
    neighbors_on = 0
    max_height, max_width = pixel_map.shape
    MIN_NEIGHBORS_ON = 3

    # Iterate through left, same column, and right of pixel_loc in each row to check if that pixel is on.
    for i in range(-1, 1):

        # Only check if the row is between 0 and the maximum width.  Otherwise, it will be out of the range of the map.
        if pixel_loc[0] + i >= 0 & pixel_loc[0] + i < max_width:

            # Iterate through from below, same row, and above pixel_loc to check if each neighbor is on.
            for j in range (-1, 1):

                # Only check if the column is between 0 and the max height.  Otherwise, it will be out of the range of
                # the map.
                if pixel_loc[1] + j >= 0 & pixel_loc[1] + j < max_height:

                    # If the pixel at the current location to check is on, increment neighbors_on.
                    if pixel_map[pixel_loc[1] + j, pixel_loc[0] + i] > threshold:
                        neighbors_on = neighbors_on + 1

    # If there are at least the minimum required neighbors on, consider this a hit.  Otherwise, it is not a hit.
    # This count includes the pixel we are currently checking.
    return neighbors_on >= MIN_NEIGHBORS_ON

##################################################
##################################################

# Function to find the uppermost pixel of the beacon at the initial estimate
# input:    initial_pixel_loc       The initial location of the beacon
#           pixel_map               Map of pixels to examine
#           threshold               Minimum value to be considered a hit
# output:   current_y               The y coordinate of the highest pixel of the beacon
def find_highest_pixel(initial_pixel_loc, pixel_map, threshold):

    # Pull off starting x and y locations
    current_x = initial_pixel_loc[0]
    current_y = initial_pixel_loc[1]
    max_y, max_x = pixel_map.shape
    found_upper_pixel = False

    # Continue trying to find a pixel farther up until we have found the uppermost pixel and next_y is not at the top
    # of the map
    while (not found_upper_pixel) and current_y < max_y-1:

        # Save the next values because we now know they are on and farther up than the previous values.
        next_y = current_y + 1

        # First try the next y.  If that is greater than the threshold and there are enough neighbors on, set up for
        # the next pixel down.
        if (pixel_map[next_y, current_x] > threshold) \
                and neighbors_are_on((current_x, next_y), pixel_map, threshold):
            current_y = next_y

        # Otherwise, if there is a pixel up and to the left, it is greater than the threshold, and it has enough
        # neighbors on, set up for that to be the next pixel to check.  current_y + 1 must be in the map because loop
        # would have ended if it was at the top.
        elif (current_x - 1 >= 0) and (pixel_map[next_y, current_x - 1] > threshold) \
                and neighbors_are_on((current_x - 1, next_y), pixel_map, threshold):
            current_y = next_y
            current_x = current_x - 1

        # Otherwise, if there is a pixel up and to the right, it is greater than the threshold, and it has enough
        # neighbors on, set up for that to be the next pixel to check.  current_y + 1 must be in the map because loop
        # would have ended if it was at the top.
        elif (current_x + 1 < max_x) and (pixel_map[next_y, current_x + 1] > threshold) \
                and neighbors_are_on((current_x + 1, next_y), pixel_map, threshold):
            current_y = next_y
            current_x = current_x + 1

        # Otherwise, we don't have any other pixels that meet the threshold.  We have found the upper pixel.
        else:
            found_upper_pixel = True

    return current_y

##################################################
##################################################

# Function to find the lowest pixel of the beacon at the initial estimate
# input:    initial_pixel_loc       The initial location of the beacon
#           pixel_map               Map of pixels to examine
#           threshold               Minimum value to be considered a hit
# output:   current_y               The y coordinate of the lowest pixel of the beacon
def find_lowest_pixel(initial_pixel_loc, pixel_map, threshold):

    # Pull off starting x and y locations
    current_x = initial_pixel_loc[0]
    current_y = initial_pixel_loc[1]
    next_y = current_y
    max_y, max_x = pixel_map.shape
    found_lower_pixel = False

    # Continue trying to find a pixel farther down until we have found the lowermost pixel and next_y is not at the
    # bottom of the map
    while (not found_lower_pixel) and current_y > 0:
        # Save the next values because we now know they are on and farther up than the previous values.
        next_y = current_y - 1

        # First try the current y.  If that is greater than the threshold and there are enough neighbors on, set up for
        # the next pixel down.
        if (pixel_map[next_y, current_x] > threshold) \
                and neighbors_are_on((current_x, next_y), pixel_map, threshold):
            current_y = next_y

        # Otherwise, if there is a pixel down and to the left, it is greater than the threshold, and it has enough
        # neighbors on, set up for that to be the next pixel to check.  current_y - 1 must be in the map because loop
        # would have ended if it was at the bottom.
        elif (current_x - 1 >= 0) and (pixel_map[next_y, current_x - 1] > threshold) \
                and neighbors_are_on((current_x - 1, next_y), pixel_map, threshold):
            current_y = next_y
            current_x = current_x - 1

        # Otherwise, if there is a pixel down and to the right, it is greater than the threshold, and it has enough
        # neighbors on, set up for that to be the next pixel to check.  current_y - 1 must be in the map because loop
        # would have ended if it was at the bottom.
        elif (current_x + 1 < max_x) and (pixel_map[next_y, current_x + 1] > threshold) \
                and neighbors_are_on((current_x + 1, next_y), pixel_map, threshold):
            current_y = next_y
            current_x = current_x + 1

        # Otherwise, we don't have any other pixels that meet the threshold.  We have found the lowest pixel.
        else:
            found_lower_pixel = True

    return current_y

##################################################
##################################################

# Function to find the rightmost pixel of the beacon at the initial estimate
# input:    initial_pixel_loc       The initial location of the beacon
#           pixel_map               Map of pixels to examine
#           threshold               Minimum value to be considered a hit
# output:   current_x               The x coordinate of the rightmost pixel of the beacon
def find_right_pixel(initial_pixel_loc, pixel_map, threshold):

    # Pull off starting x and y locations
    current_x = initial_pixel_loc[0]
    current_y = initial_pixel_loc[1]
    next_x = current_x
    max_y, max_x = pixel_map.shape
    found_right_pixel = False

    # Continue trying to find a pixel farther up until we have found the uppermost pixel and current_y is not at the top
    # of the map
    while (not found_right_pixel) and current_x < max_x-1:

        # Save the next values because we now know they are on and farther up than the previous values.
        next_x = current_x + 1

        # First try the current x.  If that is greater than the threshold and there are enough neighbors on, set up for
        # the next pixel right.
        if (pixel_map[current_y, next_x] > threshold) \
                and neighbors_are_on((current_x, current_y), pixel_map, threshold):
            current_x = next_x

        # Otherwise, if there is a pixel up and to the right, it is greater than the threshold, and it has enough
        # neighbors on, set up for that to be the next pixel to check.  current_x + 1 must be in the map because loop
        # would have ended if it was at the top.
        elif (current_y + 1 < max_y) and (pixel_map[current_y + 1, next_x] > threshold) \
                and neighbors_are_on((next_x, current_y + 1), pixel_map, threshold):
            current_x = next_x
            current_y = current_y + 1

        # Otherwise, if there is a pixel down and to the right, it is greater than the threshold, and it has enough
        # neighbors on, set up for that to be the next pixel to check.  current_x + 1 must be in the map because loop
        # would have ended if it was at the top.
        elif (current_y - 1 >= 0) and (pixel_map[current_y - 1, next_x] > threshold) \
                and neighbors_are_on((next_x, current_y - 1), pixel_map, threshold):
            current_x = next_x
            current_y = current_y - 1

        # Otherwise, we don't have any other pixels that meet the threshold.  We have found the rightmost pixel.
        else:
            found_right_pixel = True

    return current_x

##################################################
##################################################

# Function to find the leftmost pixel of the beacon at the initial estimate
# input:    initial_pixel_loc       The initial location of the beacon
#           pixel_map               Map of pixels to examine
#           threshold               Minimum value to be considered a hit
# output:   current_x               The x coordinate of the rightmost pixel of the beacon
def find_left_pixel(initial_pixel_loc, pixel_map, threshold):

    # Pull off starting x and y locations
    current_x = initial_pixel_loc[0]
    current_y = initial_pixel_loc[1]
    max_y, max_x = pixel_map.shape
    found_left_pixel = False

    # Continue trying to find a pixel farther up until we have found the uppermost pixel and current_y is not at the top
    # of the map
    while (not found_left_pixel) and current_x > 0:
        # Save the next values because we now know they are on and farther up than the previous values.
        next_x = current_x - 1

        # First try the current x.  If that is greater than the threshold and there are enough neighbors on, set up for
        # the next pixel left.
        if (pixel_map[current_y, next_x] > threshold) \
                and neighbors_are_on((current_x, current_y), pixel_map, threshold):
            current_x = next_x

        # Otherwise, if there is a pixel up and to the left, it is greater than the threshold, and it has enough
        # neighbors on, set up for that to be the next pixel to check.  current_x - 1 must be in the map because loop
        # would have ended if it was at the top.
        elif (current_y + 1 < max_y) and (pixel_map[current_y + 1, next_x] > threshold) \
                and neighbors_are_on((next_x, current_y + 1), pixel_map, threshold):
            current_x = next_x
            current_y = current_y + 1

        # Otherwise, if there is a pixel down and to the left, it is greater than the threshold, and it has enough
        # neighbors on, set up for that to be the next pixel to check.  current_x - 1 must be in the map because loop
        # would have ended if it was at the top.
        elif (current_y - 1 >= 0) and (pixel_map[current_y - 1, next_x] > threshold) \
                and neighbors_are_on((next_x, current_y - 1), pixel_map, threshold):
            current_x = next_x
            current_y = current_y - 1

        # Otherwise, we don't have any other pixels that meet the threshold.  We have found the lefttmost pixel.
        else:
            found_left_pixel = True

    return current_x


##################################################
##################################################
def hough_circles(img, blur=5, canny_thresh=200, dp=1, center_dist=200, accum=18, min_rad=0, max_rad=0, high_noise=False, show_img=False):
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
    @param high_noise uses binary thresholding to before edge detection
    @return An array of circles that should have the best matches first if there are multiple
    """
    image_scaling = 10.0
    img = np.uint8(img)

    kernel = np.ones((3, 3), np.uint8)
    # erosion followed by dilation
    img = cv2.morphologyEx(img, cv2.MORPH_OPEN, kernel)

    if show_img:
        cv2.namedWindow('before', cv2.WINDOW_NORMAL)
        cv2.resizeWindow('before', 600, 600)
        cv2.imshow('before', img)
        # cv2.imwrite('beacon_orig.png', img)
        cv2.waitKey(0)
        cv2.destroyAllWindows()

    height, width = img.shape

    resized = False

    if height < 30 or width < 30:
        img = cv2.resize(img, (int(image_scaling*width), int(image_scaling*height)), interpolation=cv2.INTER_CUBIC)
        # the interpolation adds some smoothing to the image so we reduce the rest of the smoothing
        # and lower the threshold of the edge detection
        canny_thresh /= (image_scaling / 2)
        blur = blur * 2
        if not blur % 2:
            blur += 1
        resized = True
        height, width = img.shape

    orrig_img = cv2.cvtColor(img, cv2.COLOR_GRAY2BGR)
    orrig_img = np.copy(orrig_img)
    min_rad = int(min(height, width) / 2)

    if show_img and resized:
        cv2.namedWindow('before', cv2.WINDOW_NORMAL)
        cv2.resizeWindow('before', 600, 600)
        cv2.imshow('before', img)
        # cv2.imwrite('beacon_orig.png', img)
        cv2.waitKey(0)
        cv2.destroyAllWindows()

    img = cv2.medianBlur(img, blur)
    img = cv2.GaussianBlur(img, (blur, blur), 0)

    if show_img:
        cv2.namedWindow('blur', cv2.WINDOW_NORMAL)
        cv2.resizeWindow('blur', 600, 600)
        cv2.imshow('blur', img)
        cv2.waitKey(0)
        cv2.destroyAllWindows()

    # try and figure out the proper canny threshold
    if high_noise:
        canny_thresh, img = cv2.threshold(img, 0, 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU)
    else:
        canny_thresh, _ = cv2.threshold(img, 0, 255, cv2.THRESH_BINARY | cv2.THRESH_OTSU)

    # if show_img:
    if True:
        canny_img = cv2.Canny(img, canny_thresh, canny_thresh / 15)
        cv2.namedWindow('canny', cv2.WINDOW_NORMAL)
        cv2.resizeWindow('canny', 600, 600)
        cv2.imshow('canny', canny_img)
        cv2.waitKey(0)
        # cv2.imwrite('beacon_orig_canny.png', canny_img)
        cv2.destroyAllWindows()

    timeout = 30
    while timeout:
        timeout -= 1
        accum = 1 if accum <= 0 else accum
        circles = cv2.HoughCircles(img, cv2.HOUGH_GRADIENT, dp, center_dist,
                                   param1=canny_thresh, param2=accum, minRadius=min_rad,
                                   maxRadius=max_rad)

        if high_noise and circles is not None: # and circles.size/3 < 20:
            tmp = np.uint16(np.around(circles))
            average = np.uint16(np.mean(tmp[0, :], axis=0))
            average = average[:2]/image_scaling
            return average

        if show_img and circles is not None:
            tmp = np.uint16(np.around(circles))
            average = np.uint16(np.mean(tmp[0, :], axis=0))

            tmp_img = np.copy(orrig_img)
            # for i in circles[0, :]:
            #     cv2.circle(tmp_img, (i[0], i[1]), i[2], (0, 255, 0), 2)
            #     cv2.circle(tmp_img, (i[0], i[1]), 2, (0, 0, 255), 3)

            cv2.circle(tmp_img, (average[0], average[1]), average[2], (255, 0, 0), 2)

            cv2.namedWindow('Detected circles', cv2.WINDOW_NORMAL)
            cv2.resizeWindow('Detected circles', 600, 600)
            cv2.imshow('Detected circles', tmp_img)
            # cv2.imwrite('beacon_orig_circle.png', img_circle)
            cv2.waitKey(0)
            cv2.destroyAllWindows()

        if circles is None:
            # print "unable to find any circles"
            accum -= 1
            center_dist += 5
            continue

        num_circles = circles.size / 3
        # print num_circles
        if num_circles > 1:
            accum += 1
        elif num_circles is 1:
            break
        else:
            accum -= 1
            center_dist += 5

    if not timeout:
        return None

    for i in circles[0, :]:
        cv2.circle(orrig_img, (i[0], i[1]), i[2], (0, 255, 0), 2)
        cv2.circle(orrig_img, (i[0], i[1]), 2, (0, 0, 255), 3)

    if show_img:
        cv2.namedWindow('Detected circles', cv2.WINDOW_NORMAL)
        cv2.resizeWindow('Detected circles', 600, 600)
        cv2.imshow('Detected circles', orrig_img)
        # cv2.imwrite('beacon_orig_circle.png', img_circle)
        cv2.waitKey(0)
        cv2.destroyAllWindows()

    circles = np.squeeze(circles)

    do_plots = False
    if do_plots == True:
        print '\nInternal Hough Transform Output: ', circles
        fig, ax = plt.subplots(1)
        plt.imshow(img, interpolation='none', cmap='viridis')
        # plt.scatter(circles[0], circles[1], marker='x', s=150, linewidth=2, c='r')
        # circ = plt.Circle((circles[0], circles[1]), circles[2], color='w', fill=False)
        # ax.add_patch(circ)
        # plt.savefig('centerfinding_roi.png')
        plt.show()

    if resized is True:
        circles[:2] = circles[:2] / image_scaling

    return circles

##################################################
##################################################

# Output:   loc_centroid        (x, y) pixel/line coordinate of centroid locations
#           DN                  total signal count of ROI

def find_centroid_point_source(pixel_map, pixel_line_beacon_i, ROI_parameters):

    loc_centroid = np.empty([len(pixel_line_beacon_i)], dtype=tuple)
    DN = np.empty([len(pixel_line_beacon_i)], dtype=int)

    # crop original image to an ROI based on initial
    corner_ROI, image_ROI = generate_point_source_ROI(pixel_map, pixel_line_beacon_i, ROI_parameters)

    for i in range(0, len(image_ROI)):

        # determine average value of region of interest border, subtract from rest of pixel map
        image_ROI[i] = apply_ROI_border(image_ROI[i], ROI_parameters)

        # calculate centroid position and ROI brightness value <-- centroid location is pixel number not index value
        # (starts with 1)
        loc_centroid[i], DN[i] = find_centroid(image_ROI[i], corner_ROI[i], ROI_parameters)

    return loc_centroid, DN



def find_center_resolved_body(pixel_map, pixel_line_beacon_i, ROI_parameters):

    loc_center = np.empty([len(pixel_line_beacon_i)], dtype=tuple)
    DN = np.empty([len(pixel_line_beacon_i)], dtype=int)

    # crop original image to an ROI based on initial
    corner_ROI, image_ROI = generate_point_source_ROI(pixel_map, pixel_line_beacon_i, ROI_parameters)


    for i in range(0, len(image_ROI)):

        # determine average value of region of interest border, subtract from rest of pixel map
        image_ROI[i] = apply_ROI_border(image_ROI[i], ROI_parameters)


        # plt.savefig('saved_output/cropped_image_' + str(i) + '.png')
        # np.savez('saved_output/cropped_image_' + str(i)  + '.npz')

        maxROIvalue = np.amax(image_ROI[i])
        # print maxROIvalue

        # normalize ROI for center-finding function
        currentROI = (image_ROI[i] / maxROIvalue) * 255.

        center = hough_circles(currentROI, center_dist=50, canny_thresh=175, blur=5, accum=5, show_img=False)

        if center is not None:
            loc_center[i] = (center[0] + corner_ROI[i][1], center[1] + corner_ROI[i][0])
        else:
            loc_center = None

    return loc_center, DN

##################################################
##################################################

# Function to determine if the body at the entered pixel location has already been discovered
# input:    pixel_loc               The pixel location to check
#           corner_ROI              List of corner ROIs found so far
#           image_ROI               List of images for ROIs found so far
# output:   alreadyFound            True if the pixel location is inside one of the previously found ROIs and false
#                                   otherwise

def already_found_ROI(pixel_loc, corner_ROI, image_ROI):

    # Initially have not found it
    alreadyFound = False
    i = 0

    # Iterate through by checking the pixel location against each corner_ROI
    while (i < len(corner_ROI) and not alreadyFound):

        #print "\n\nLoop ", i
        #print "\n", pixel_loc
        #print corner_ROI[i]
        #print image_ROI[i].shape

        beaconLengthY, beaconLengthX = image_ROI[i].shape

        # OFFSET - Included offset of 2 pixels to eliminate close proximity values
        # Check x location of pixel to see if it is between x of corner_ROI and x of corner_ROI plus x distance of
        # image_ROI
        if (pixel_loc[1] >= corner_ROI[i][1] - 2) and (pixel_loc[1] < corner_ROI[i][1] + beaconLengthX + 2):

            # Check y location of pixel to see if it is between y of corner_ROI and y of corner_ROI plus y distance of
            # image_ROI
            if (pixel_loc[0] >= corner_ROI[i][0] - 2) and (pixel_loc[0] < corner_ROI[i][0] + beaconLengthY + 2):

                # This ROI has already been found.  Disregard it.
                #print "Found the thing!"
                alreadyFound = True

        # Increment i
        i = i + 1

    #print "\n\n"

    return alreadyFound


##################################################
##################################################

def checkBeaconStatus(posBeacons, posSC, DCM_BN, camFoV, radiusBeacons, angularSizeMin):
    """Check current status of beacons in target beacon lists
       :param pos_cb: list of target beacon position (1x3) numpy arrays in heliocentric coordinates [km]
       :param pos_obs: position of camera in heliocentric coordinates [km]
       :param DCM_NB: direction cosine matrix of observer attitude [body to heliocentric coord. frame]
       :param radiusBeacons: field of view tuple (horizontal, vertical) [degrees]
       :param angularSizeMin: minimum beacon size for it to be treated as a resolved body [deg]
       :return: list of beacon status as defined below
                        '0' if beacon is not in camera field of view
                        '1' if beacon is in camera field of view as a point source
                        '2' if beacon is in camera field of view as a resolved body
       """

    numBeacons = len(posBeacons)
    beaconStatus = []

    for ind in range(numBeacons):

        #################################################
        # Determine if beacon is in camera field of view

        # compute relative position of celestial body to observer in body coordinates
        pos_obs2cb_helio = (posBeacons[ind] - posSC)
        pos_obs2cb_body = np.matmul(DCM_BN, pos_obs2cb_helio)
        e_obs2cb_body = pos_obs2cb_body / np.linalg.norm(pos_obs2cb_body)

        # take field of view angles from centerline
        horiz_fovcenter = camFoV[0] / 2.
        vert_fovcenter = camFoV[1] / 2.

        # compute azimuth and elevation of celestial body
        az = math.degrees(math.atan2(e_obs2cb_body[1], e_obs2cb_body[0]))
        el = math.degrees(math.asin(e_obs2cb_body[2]))

        # check if limits in camera specs
        chkFoV = True
        if (abs(az) > horiz_fovcenter) or (abs(el) > vert_fovcenter):
            chkFoV = False
            currentStatus = 0

        if np.linalg.norm(pos_obs2cb_body) < radiusBeacons[ind]:
            chkFoV = False
            currentStatus = 0

        #################################################
        # Determine if beacon should be treated as a point-source or a resolved body

        if chkFoV == True:

            angularSizeBeacon = 2*math.degrees(math.atan(radiusBeacons[ind]/np.linalg.norm(pos_obs2cb_helio)))

            # print 'Check Beacon Status Angular Size: ', angularSizeBeacon

            if angularSizeBeacon < angularSizeMin:
                currentStatus = 1
            else:
                currentStatus = 2

        beaconStatus.append(currentStatus)

    return beaconStatus
