# Title   : find_centroid_functions.py
# Author  : Joe Park
# Date    : 08/22/2017
# Synopsis: Finds centroid of a point source from a pixel map image. DINO C-REx module.

import math         # common math functions
import numpy as np  # matrix algebra
import sqlite3 as sql
import dynamics as dyn
import matplotlib.pyplot as plt


##################################################
##################################################

# Function to convert s/c state and beacon position to initial estimate of beacon centroid pixel/line location
# (may not be necessary if pixel line location can be directly provided by navigation module)
# input:    pos_beacon              position of beacon in inertial heliocentric coordinates [km]
#                                   (n x 3) array
#           pos_sc                  position of s/c in inertial heliocentric coordinates [km]
#           attde_sc                BN direction cosine matrix of s/c attitude
#                                   (b1 = camera bore-sight, b2 = left, b3 = up)
#           cam_res                 Resolution of camera sensor (horizontal, vertical)
#           cam_focal_length        Camera focal length [mm]
#           cam_pixel_size          Camera pixel size (width, height) [mm]

def initial_beacons_estimate(pos_beacon, pos_sc, attde_sc, cam_res, cam_focal_length, cam_pixel_size):

    n_beacons = pos_beacon.shape[0]

    # image upper left corner distance
    x_upperleft = cam_res[0] / 2. * cam_pixel_size[0]
    y_upperleft = cam_res[1] / 2. * cam_pixel_size[1]

    # cycle through each beacon position and determine pixel / line estimate
    pixel = []
    line = []
    for ind_beacon in range(n_beacons):

        # compute unit direction from s/c to beacon in camera body coordinates
        print ind_beacon
        print pos_beacon[ind_beacon, :]
        current_pos_beacon = pos_beacon[ind_beacon, :]
        pos_sc2beacon = current_pos_beacon - pos_sc
        e_sc2beacon = pos_sc2beacon / np.linalg.norm(pos_sc2beacon)
        e_sc2beacon_cam = np.matmul(attde_sc, e_sc2beacon)

        # convert to pixel/line coordinates and shift origin to upper left corner
        x_focal = -(cam_focal_length/e_sc2beacon_cam[0]) * e_sc2beacon_cam[1] + x_upperleft
        y_focal = -(cam_focal_length/e_sc2beacon_cam[0]) * e_sc2beacon_cam[2] + y_upperleft

        # convert to pixel units
        current_pixel = x_focal / cam_pixel_size[0]
        current_line = y_focal / cam_pixel_size[1]

        print '\nS/C to Beacon Unit Vector (HCI): ', e_sc2beacon
        print 'S/C to Beacon Unit Vector (body coord): ',e_sc2beacon_cam
        print 'Pixel/Line Location: ', current_pixel, current_line

        # check to see if beacon is in field of view
        if current_pixel >= 0 and current_pixel <= cam_res[0] and \
            current_line >= 0 and current_line <= cam_res[1]:
            pixel.append(current_pixel)
            line.append(current_line)

    #if len(pixel) != 0:
    #    pixel_out = np.resize(pixel, (1, n_beacons))
    #    line_out = np.resize(line, (1, n_beacons))
    #else:
    #    pixel_out = np.empty((1,0))
    #    line_out = np.empty((1,0))

    return (pixel, line)
    #return (pixel_out, line_out)


##################################################
##################################################

# Function to convert s/c state and beacon position to initial estimate of beacon centroid pixel/line location
# (may not be necessary if pixel line location can be directly provided by navigation module)
# NOTE CURRENTLY BROKEN FOR HORIZONTAL FIELDS OF VIEW!!!!!!!!!!!!!!!
# Input:
#           attde_sc                BN direction cosine matrix of s/c attitude
#                                   (b1 = camera boresight, b2 = left, b3 = up)
#           cam_res                 Resolution of camera sensor (horizontal, vertical)
#           cam_focal_length        Camera focal length [mm]
#           cam_pixel_size          Camera pixel size (width, height) [mm]
#           fname_catalog           Filename for star catalog file to be used
# Output:
#           (pixel, line)           Tuple of 1D arrays of pixel and line coordinates

def initial_stars_estimate(attde_sc, cam_res, cam_focal_length, cam_pixel_size, fname_catalog):

    # Currently utilizes subset of Tycho2 catalog (average error for RA and Dec is 10E-8 degrees)

    # generate additional camera parameters
    cam_sensor_size = (cam_res[0] * cam_pixel_size[0], cam_res[1] * cam_pixel_size[1])  # [mm]
    cam_fov = (2 * math.degrees(math.atan2(cam_sensor_size[0] / 2., cam_focal_length)),
               2 * math.degrees(math.atan2(cam_sensor_size[1] / 2., cam_focal_length)))

    print 'Cam FoV (Initial Star Estimate Function): ', cam_fov

    # calculate unit vectors to four corners of camera field of view in camera body frame
    hplane = math.cos(math.radians(cam_fov[1]/2.))      # projection of unit vector on horizontal plane
    vplane = math.cos(math.radians(cam_fov[0]/2.))
    e2 = hplane * math.sin(math.radians(cam_fov[0]/2.))
    # e3 = vplane * math.sin(math.radians(cam_fov[1]/2.))
    e1 = hplane * math.cos(math.radians(cam_fov[0]/2.))
    e3 = math.sqrt(1 - (e1**2 + e2**2))

    upperleft = np.array([e1, e2, e3])
    upperright = np.array([e1, -e2, e3])
    lowerleft = np.array([e1, e2, -e3])
    lowerright = np.array([e1, -e2, -e3])

    # print '\nFoV Corner Unit Vectors (body coord.):'
    # print upperleft, np.linalg.norm(upperleft)
    # print upperright, np.linalg.norm(upperright)
    # print lowerleft, np.linalg.norm(lowerleft)
    # print lowerright, np.linalg.norm(lowerright)
    # print 'Upper Corner Angular Separation: ', math.degrees(math.acos(np.dot(upperleft, upperright)))
    # print 'Lower Corner Angular Separation: ', math.degrees(math.acos(np.dot(lowerleft, lowerright)))
    # print 'Right Side Angular Separation: ', math.degrees(math.acos(np.dot(upperright, lowerright)))
    # print 'Left Side Angular Separation: ', math.degrees(math.acos(np.dot(upperleft, lowerleft)))

    # convert unit vectors from body frame to heliocentric inertial frame
    ul_helio = np.matmul(attde_sc.T, upperleft)
    ur_helio = np.matmul(attde_sc.T, upperright)
    ll_helio = np.matmul(attde_sc.T, lowerleft)
    lr_helio = np.matmul(attde_sc.T, lowerright)

    # calculate inertial RA and Dec from inertial unit vectors
    ra_N_ul = math.degrees(math.atan2(ul_helio[1], ul_helio[0]))
    ra_N_ur = math.degrees(math.atan2(ur_helio[1], ur_helio[0]))
    ra_N_ll = math.degrees(math.atan2(ll_helio[1], ll_helio[0]))
    ra_N_lr = math.degrees(math.atan2(lr_helio[1], lr_helio[0]))

    # wrap right ascension so that it only has positive values
    if ra_N_ul < 0:
        ra_N_ul = ra_N_ul + 360.
    if ra_N_ll < 0:
        ra_N_ll = ra_N_ll + 360.
    if ra_N_ur < 0:
        ra_N_ur = ra_N_ur + 360.
    if ra_N_lr < 0:
        ra_N_lr = ra_N_lr + 360.

    dec_N_ul = math.degrees(math.asin(ul_helio[2] / np.linalg.norm(ul_helio)))
    dec_N_ur = math.degrees(math.asin(ur_helio[2] / np.linalg.norm(ur_helio)))
    dec_N_ll = math.degrees(math.asin(ll_helio[2] / np.linalg.norm(ll_helio)))
    dec_N_lr = math.degrees(math.asin(lr_helio[2] / np.linalg.norm(lr_helio)))

    radec_N_ul = (ra_N_ul, dec_N_ul)
    radec_N_ur = (ra_N_ur, dec_N_ur)
    radec_N_ll = (ra_N_ll, dec_N_ll)
    radec_N_lr = (ra_N_lr, dec_N_lr)

    # find stars in catalog within field of view
    radec_stars = find_stars_in_FoV((radec_N_ul, radec_N_ur, radec_N_ll, radec_N_lr), fname_catalog)

    # convert stars in catalog to pixel / line values
    pixel, line = radec_to_pixelline(
        radec_stars, attde_sc, cam_res, cam_focal_length, cam_pixel_size)

    n_stars = len(pixel)
    pixel_out = np.resize(pixel, (1, n_stars))
    line_out = np.resize(line, (1, n_stars))

    return (pixel_out, line_out, radec_stars)


##################################################
##################################################

# Function to convert right ascension and declination into pixel / line estimates
# Note this currently returns all entries in inertial normal rectangular field of view (does not account for rotation)
# Inputs:   radec_corners   list of tuple values for each corner right ascension and declination value
#                           [deg, deg] (inertial coord.)
#           fname_catalog   filename of star catalog sql database file
# Outputs:  radec_inertial  tuple of ra, dec, and id lists

def find_stars_in_FoV(radec_corners, fname_catalog):

    M_CUTOFF = 7.

    # open star catalog
    star_catalog = sql.connect(fname_catalog)
    s = star_catalog.cursor()
    s.execute("SELECT * FROM tycho_data")

    # find min and max ra dec values in field of view
    ra_min = min(radec_corners[0][0], radec_corners[1][0], radec_corners[2][0], radec_corners[3][0])
    ra_max = max(radec_corners[0][0], radec_corners[1][0], radec_corners[2][0], radec_corners[3][0])
    dec_min = min(radec_corners[0][1], radec_corners[1][1], radec_corners[2][1], radec_corners[3][1])
    dec_max = max(radec_corners[0][1], radec_corners[1][1], radec_corners[2][1], radec_corners[3][1])

    # check if field of view covers inertial poles (will determine if max declination needs to extend to cover poles)
    # (process assumes max angular field of view in either direction is 180 degrees)
    if ra_max - ra_min  > 180:

        # check which pole the field of view covers
        # +k covered
        if abs(dec_min) < abs(dec_max): dec_max = 90

        # -k covered
        elif abs(dec_min) > abs(dec_max): dec_min = -90

        else:
            print "ERROR: Star Catalog Field of View Min/Max RA and Dec Calculation"


    # print '\nInertial RA and Dec for star catalog lookup'
    # print radec_corners[0]
    # print radec_corners[1]
    # print radec_corners[2]
    # print radec_corners[3]

    # remove star catalog entries outside of rectangular field of view
    s.execute("SELECT * FROM tycho_data WHERE RA BETWEEN (?) AND (?)"
              " AND DEC BETWEEN (?) AND (?)"
              " AND BTmag IS NOT NULL "
              " AND BTmag <= (?)"
              " ORDER BY BTmag ASC",
              (ra_min, ra_max, dec_min, dec_max, M_CUTOFF))
    rows = s.fetchall()
    n_stars = len(rows)

    print '\nRA and Dec limits ', (ra_min, ra_max), (dec_min, dec_max)
    print 'Number of catalog entries: ', n_stars

    # generate numpy array of results
    ra = []
    dec = []
    id = []
    for ind in range(n_stars):
        ra.append(rows[ind][1])
        dec.append(rows[ind][2])
        id.append(rows[ind][0])
        if ind < 4:
            print rows[ind][0], rows[ind][1], rows[ind][2], rows[ind][3]

    radec_inertial = (ra, dec, id)

    return radec_inertial


##################################################
##################################################

# Function to convert right ascension and declination into pixel / line estimates
# Note:     Currently takes min/max RA and Dec values and calculates pixel / line estimates.
#           Only those within sensor dimensions are returned.
#
# Inputs:   radec_inertial      [n x 2] array of right ascension and declination of stars in field of view [deg, deg]
#           attde_sc            BN transformation matrix
#           cam_res             camera sensor resolution [horiz, vertical]
#           cam_focal_length    camera focal length [mm]
#           cam_pixel_size      individual pixel size [horiz mm, vertical mm]
# Outputs:  (pixel, line)       tuple of 1D numpy arrays for pixel and line coordinates

def radec_to_pixelline(radec_inertial, attde_sc, cam_res, cam_focal_length, cam_pixel_size):

    n_stars = len(radec_inertial[0])

    # image upper left corner distance
    x_upperleft = cam_res[0] / 2. * cam_pixel_size[0]
    y_upperleft = cam_res[1] / 2. * cam_pixel_size[1]

    # cycle through each star position and determine pixel / line estimate
    pixel = []
    line = []
    for ind in range(n_stars):

        # compute unit direction from s/c to star in camera body coordinates
        ra = radec_inertial[0][ind]
        dec = radec_inertial[1][ind]

        ehat = dyn.radec_to_unitvector((ra,dec))
        e_star_helio = np.array([ehat[0], ehat[1], ehat[2]])

        e_star_cam = np.matmul(attde_sc, e_star_helio)

        #print 'Unit vector to star (HCI): ', e_star_helio
        #print 'Unit vector to star (camera coord.): ', e_star_cam
        e1_cam_helio = np.matmul(attde_sc.T, np.array([1,0,0]))
        # print 'e1 camera unit vector (HCI): ', e1_cam_helio
        # print 'Angular separation: ', math.degrees(math.acos(np.dot(e_star_helio, e1_cam_helio)))

        # convert to pixel/line coordinates and shift origin to upper left corner
        x_focal = -(cam_focal_length / e_star_cam[0]) * e_star_cam[1] + x_upperleft
        y_focal = -(cam_focal_length / e_star_cam[0]) * e_star_cam[2] + y_upperleft

        current_pixel = x_focal / cam_pixel_size[0]
        current_line = y_focal / cam_pixel_size[1]

        #print 'Current RA / Dec: ', ra, dec
        #print 'Current pixel/line coordinate: ', current_pixel, current_line

        # check to make sure pixel and line coordinates are within sensor limits
        if current_pixel <= cam_res[0] and current_pixel >= 0 \
                and current_line <= cam_res[1] and current_line >= 0:
            pixel.append(current_pixel)
            line.append(current_line)
            #print 'STAR IN FIELD OF VIEW FOUND'

    pixel_out = np.resize(pixel, (1, n_stars))
    line_out = np.resize(line, (1, n_stars))

    return (pixel, line)

