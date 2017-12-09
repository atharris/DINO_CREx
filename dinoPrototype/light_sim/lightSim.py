#	Title   :   lightSim.py
#	Author  :   Joe Park
#	Date    :   03/19/2017
#	Synopsis:   Runs scenarios for lighting simulation of DINO C-REx
#               For a given scenario, simulates solar illumination of a celestial body
#
#   Parameters      :   Celestial Body (cb) Parameters [label, albedo, radius, position]
#                       Spacecraft (s/c) Parameters [attitude, position]
#                       Camera Parameters [horizontal field of view, vertical field of view]
#                       Spherical Mapping Parameters [latitude resolution, longitude resolution]
#   Dependencies    :   CELESTIAL_BODY_PARAMETERS.csv
#                       /functions/lightSimFunctions.py
#                       /functions/lightSimPlots.py
#                       packages - csv, math, matplotlib, numpy


import numpy as np
import sys
sys.path.insert(0, 'functions')
import lightSimFunctions as lSim
import csv


###################################################
###################################################

# Constants
AU = 149597870.7            # km

###################################################
###################################################

# Define and create a database of necessary parameters for celestial bodies of interest (label, albedo, radius [km])
# Values pulled from https://nssdc.gsfc.nasa.gov/planetary/planetfact.html (accessed 2/20/2017).

cb_names = ('Mercury', 'Venus', 'Earth', 'Mars')

CB_params = [['Mercury', .142, 2439.7],
             ['Venus', .689, 6051.8],
             ['Earth', .434, 6378.137],
             ['Mars', .170, 3396.2],
             ['Jupiter', .538, 71492],
             ['Earth Moon', .12, 1738.1],
             ['Phobos', .07, 13],
             ['Deimos', .08, 7.8],
             ['Ceres', .09, 487.3]]

with open("CELESTIAL_BODY_PARAMETERS.csv", "w") as output:
    writer = csv.writer(output)
    writer.writerows(CB_params)


###################################################
###################################################

# angular field of view for camera
dtheta_horiz = 45.           # [deg]
dtheta_vert = 45.            # [deg]

# fidelity of spherical mapping used in lighting simulation
lat_res = 120               # [n pts]
long_res = 120              # [n pts]


####################################################
####################################################

doSingleCB = True

if doSingleCB:

    # s/c attitude, DCM for camera body to heliocentric coordinates (BN matrix)
    # camera body frame: {+i right, +j forward, +k vertical}
    DCM_camera = np.eye(3)

    # 45 degree rotation about z-axis
    DCM_camera = np.array([[0.70710678, 0.70710678, -0.],
                           [-0.70710678, 0.70710678, 0.],
                           [0., 0., 1.]])

    # 90 degree rotation about z-axis
    #DCM_camera = np.array([[0, 1, 0],
    #                  [-1, 0, 0],
    #                  [0, 0, 1]])

    pos_camera = np.array([AU, 0, 0])           # [m]
    pos_cbs = np.array([[2.2*AU, AU, .2*AU]])       # Mercury
    do_pt_source = [True]

    # conduct simulation
    lsim_output = lSim.lightSim(DCM_camera, pos_camera, pos_cbs, (dtheta_horiz, dtheta_vert),
                                lat_res, long_res, do_pt_source)

    if len(lsim_output) > 0:

        output1 = lsim_output[cb_names[0]]

        print cb_names[0], ' Light Sim Output'
        print 'RA: ', output1['facetRA']
        print 'Dec: ', output1['facetDec']
        print 'Pt Source Body Position: ', output1['bodypos']
        print 'Facet Area Min and Max: ', min(output1['facetArea']), max(output1['facetArea'])
        print 'Net Albedo Min and Max: ', min(output1['netAlbedo']), max(output1['netAlbedo'])

    else:

        print 'No visible objects in field of view.'






