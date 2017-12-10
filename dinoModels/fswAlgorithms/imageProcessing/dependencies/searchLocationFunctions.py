# Title   : find_centroid_functions.py
# Author  : Joe Park
# Date    : 08/22/2017
# Synopsis: Finds centroid of a point source from a pixel map image. DINO C-REx module.

import math         # common math functions
import numpy as np  # matrix algebra
import sqlite3 as sql
import dynamicFunctions as dyn
import matplotlib.pyplot as plt

import sys, os, inspect
# sys.path.append('../../../../../external/')
import objectIDFunctions as oID


##################################################
##################################################

def initial_beacons_estimate(N_r_beacons, N_r_sc, BN_dcm_sc, cameraParameters):
    """
    Generates initial pixel, line coordinate estimates for beacons in the camera field of view.
    @param N_r_beacons: list of (1x3) numpy arrays of beacons positions in HCI coordinate frame [km]
    @param N_r_sc: (1x3) numpy array of spacecraft position in HCI coordinate frame [km]
    @param BN_dcm_sc: (3x3) numpy array of the BN direction cosine matrix (HCI to camera body coord. frame)
    @param cameraParameters: python dict with the following entries
                            ['resolution'] tuple of horizontal x vertical camera sensor resolution
                            ['focal length'] camera sensor effective focal length [m]
                            ['pixel size'] tuple of horizontal x vertical camera sensor pixel size [m]
    @return pixelLine: list of (pixel, line) coordinate initial estimates tuples
    """

    cam_res = cameraParameters['resolution']
    cam_focal_length = cameraParameters['focal length']
    cam_pixel_size = cameraParameters['pixel size']

    n_beacons = len(N_r_beacons)

    # image upper left corner distance
    x_upperleft = cam_res[0] / 2. * cam_pixel_size[0]
    y_upperleft = cam_res[1] / 2. * cam_pixel_size[1]

    # cycle through each beacon position and determine pixel / line estimate
    pixel_line = []
    for ind_beacon in range(n_beacons):

        # compute unit direction from s/c to beacon in camera body coordinates
        N_r_currentBeacon = N_r_beacons[ind_beacon]
        N_r_sc2beacon = N_r_currentBeacon - N_r_sc
        N_e_sc2beacon = N_r_sc2beacon / np.linalg.norm(N_r_sc2beacon)
        B_e_sc2beacon = np.matmul(BN_dcm_sc, N_e_sc2beacon)

        # convert to pixel/line coordinates and shift origin to upper left corner
        x_focal = -(cam_focal_length/B_e_sc2beacon[0]) * B_e_sc2beacon[1] + x_upperleft
        y_focal = -(cam_focal_length/B_e_sc2beacon[0]) * B_e_sc2beacon[2] + y_upperleft

        # convert to pixel units
        current_pixel = x_focal / cam_pixel_size[0]
        current_line = y_focal / cam_pixel_size[1]

        # check to see if beacon is in field of view
        if current_pixel >= 0 and current_pixel <= cam_res[0] and \
            current_line >= 0 and current_line <= cam_res[1]:
            pixel_line.append((current_pixel, current_line))

    return pixel_line


##################################################
##################################################

def initial_stars_estimate(BN_dcm_sc, cameraParameters, fname_catalog):
    """
    Generates initial pixel, line coordinate estimates for stars in the camera field of view using a reference catalog
    @param BN_dcm_sc: (3x3) numpy array of the BN direction cosine matrix (HCI to camera body coord. frame)
    @param cameraParameters: python dict with the following entries
                            ['resolution'] tuple of horizontal x vertical camera sensor resolution
                            ['focal length'] camera sensor effective focal length [m]
                            ['pixel size'] tuple of horizontal x vertical camera sensor pixel size [m]
                            ['field of view'] tuple of horizontal x vertical camera field of view [deg]
    @param fname_catalog: reference catalog filename in '.db' SQL database format with the following fields
                            'id' unique ID
                            'BTmag' visual magnitude
                            'RA' HCI right ascension [deg]
                            'DEC' HCI declination [deg]
    @return pixelLine: list of (pixel, line) coordinate initial estimates tuples
    @return referenceID: list of reference catalog IDs of stars in field of view
    """

    # Currently utilizes subset of Tycho2 catalog (average error for RA and Dec is 10E-8 degrees)

    cam_res = cameraParameters['resolution']
    cam_focal_length = cameraParameters['focal length']
    cam_pixel_size = cameraParameters['pixel size']
    cam_fov = cameraParameters['field of view']

    # print 'Cam FoV (Initial Star Estimate Function): ', cam_fov

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

    print '\nVectors to Corners in Body Coord.'
    print upperleft
    print upperright
    print lowerleft
    print lowerright

    print '\nBN DCM Matrix'
    print BN_dcm_sc
    # convert unit vectors from body frame to heliocentric inertial frame
    ul_helio = np.matmul(BN_dcm_sc.T, upperleft)
    ur_helio = np.matmul(BN_dcm_sc.T, upperright)
    ll_helio = np.matmul(BN_dcm_sc.T, lowerleft)
    lr_helio = np.matmul(BN_dcm_sc.T, lowerright)

    print '\nUL, UR, LL, LR Corner in HCI'
    print ul_helio
    print ur_helio
    print ll_helio
    print lr_helio

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
    pixel_line = radec_to_pixelline(
        radec_stars, BN_dcm_sc, cameraParameters)

    catalog_id = radec_stars[2]

    return pixel_line, catalog_id


##################################################
##################################################

def find_stars_in_FoV(radec_corners, fname_catalog):
    """
    Returns all stars in field of view in a reference star catalog based on right ascension and declination bounds.
    @param radec_corners: list of tuple values for each right ascension and declination corner point (ra, dec) [deg]
                            (upper left, upper right, lower left, lower right)
    @param fname_catalog: reference catalog filename in '.db' SQL database format with the following fields
                            'id' unique ID
                            'BTmag' visual magnitude
                            'RA' HCI right ascension [deg]
                            'DEC' HCI declination [deg]
    @return ra: list of right ascension matches [deg]
    @return dec: list of declination matches [deg]
    @return id: list of reference catalog id matches
    """

    # open star catalog
    # print fname_catalog
    # os.getcwd()
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

    # check if the delta between the corner RA values exceed 135 degrees
    # (assume field of view cannot exceed 180 degrees)
    deltaRA1 = abs(radec_corners[0][0] - radec_corners[3][0])
    deltaRA2 = abs(radec_corners[1][0] - radec_corners[2][0])
    deltaRA = ra_max-ra_min
    if deltaRA1 > 180:
        deltaRA1 = 360 - deltaRA

    if deltaRA2 > 180:
        deltaRA2 = 360 - deltaRA2

    # if deltaRA > 90:
    if deltaRA1 > 135 or deltaRA2 > 135:

        # check which pole the field of view covers
        # +k covered
        if abs(dec_min) < abs(dec_max): dec_max = 90

        # -k covered
        elif abs(dec_min) > abs(dec_max): dec_min = -90

        else:
            print "ERROR: Star Catalog Field of View Min/Max RA and Dec Calculation"


    print '\nInertial RA and Dec for star catalog lookup'
    print radec_corners[0]
    print radec_corners[1]
    print radec_corners[2]
    print radec_corners[3]
    print ra_min, ra_max
    print dec_min, dec_max

    # remove star catalog entries outside of rectangular field of view
    s.execute("SELECT * FROM tycho_data WHERE RA BETWEEN (?) AND (?)"
              " AND DEC BETWEEN (?) AND (?)"
              " AND VTmag IS NOT NULL "
              " ORDER BY VTmag ASC",
              (ra_min, ra_max, dec_min, dec_max))
    rows = s.fetchall()
    n_stars = len(rows)

    # print '\nRA and Dec limits ', (ra_min, ra_max), (dec_min, dec_max)
    # print 'Number of catalog entries: ', n_stars

    # generate numpy array of results
    ra = []
    dec = []
    id = []

    for ind in range(n_stars):
        ra.append(rows[ind][1])
        dec.append(rows[ind][2])
        id.append(rows[ind][0])

    return (ra, dec, id)


##################################################
##################################################

def radec_to_pixelline(N_radec, BN_dcm_cam, cameraParameters):
    """
    Converts right ascension and declination values to pixel line coordinates.
    @param N_radec: list of tuple values for each HCI right ascension and declination coordinate [deg]
    @param BN_dcm_cam: reference catalog filename in '.db' SQL database format with the following fields
    @param BN_dcm_sc: (3x3) numpy array of the BN direction cosine matrix (HCI to camera body coord. frame)
    @param cameraParameters: python dict with the following entries
                            ['resolution'] tuple of horizontal x vertical camera sensor resolution
                            ['focal length'] camera sensor effective focal length [m]
                            ['pixel size'] tuple of horizontal x vertical camera sensor pixel size [m]
                            ['field of view'] tuple of horizontal x vertical camera field of view [deg]
    @return pixel_line: list of (pixel, line) coordinate tuples
    """

    cam_res = cameraParameters['resolution']
    cam_focal_length = cameraParameters['focal length']
    cam_pixel_size = cameraParameters['pixel size']

    n_stars = len(N_radec[0])

    # image upper left corner distance
    x_upperleft = cam_res[0] / 2. * cam_pixel_size[0]
    y_upperleft = cam_res[1] / 2. * cam_pixel_size[1]

    # cycle through each star position and determine pixel / line estimate
    pixel_line = []
    for ind in range(n_stars):

        # compute unit direction from s/c to star in camera body coordinates
        ra = N_radec[0][ind]
        dec = N_radec[1][ind]

        ehat = dyn.radec_to_unitvector((ra,dec))
        e_star_helio = np.array([ehat[0], ehat[1], ehat[2]])

        e_star_cam = np.matmul(BN_dcm_cam, e_star_helio)

        # print 'Unit vector to star (HCI): ', e_star_helio
        # print 'Unit vector to star (camera coord.): ', e_star_cam
        # e1_cam_helio = np.matmul(BN_dcm_cam.T, np.array([1,0,0]))
        # print 'e1 camera unit vector (HCI): ', e1_cam_helio
        # print 'Angular separation: ', math.degrees(math.acos(np.dot(e_star_helio, e1_cam_helio)))

        # convert to pixel/line coordinates and shift origin to upper left corner
        x_focal = -(cam_focal_length / e_star_cam[0]) * e_star_cam[1] + x_upperleft
        y_focal = -(cam_focal_length / e_star_cam[0]) * e_star_cam[2] + y_upperleft

        current_pixel = x_focal / cam_pixel_size[0]
        current_line = y_focal / cam_pixel_size[1]

        # check to make sure pixel and line coordinates are within sensor limits
        if current_pixel <= cam_res[0] and current_pixel >= 0 \
                and current_line <= cam_res[1] and current_line >= 0:
            pixel_line.append((current_pixel, current_line))
            #print 'STAR IN FIELD OF VIEW FOUND'


    # for ind in range(len(pixel_line)):
    #     print 'Check: ', N_radec[0][ind], N_radec[1][ind]
    #     print 'Verify: ', oID.pixelline_to_radec(pixel_line[ind],
    #                                              BN_dcm_cam,
    #                                              cameraParameters)

    return pixel_line

