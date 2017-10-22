#!/usr/local/bin/python
'''

 <<Description>>


 <<Summary>>

'''

__author__ = 'Marc Balducci'
__version__ = '$Revision$'[11:-2]
__date__ = '$Date$'[7:26]

################################################################################
#                     I M P O R T     L I B R A R I E S
################################################################################

import numpy as np
import pdb
from numpy import cos as cos
from numpy import sin as sin


################################################################################
#                  I N T E R N A L    F U N C T I O N S:
################################################################################

def I_to_TV(input):
	# pull out the angles from the input (RA, dec, twist) in that order
    #pdb.set_trace()
    angles = input
    I_to_TV_matrix = np.dot( rot3( angles[0] ), \
                     np.dot( rot2( angles[1] ), rot1( angles[2] ) ) )

    return I_to_TV_matrix

def rot1(input):
    # rotation matrix about axis 1
    angle = input
    rot1_matrix = np.array( [ [ 1., 0., 0. ], [ 0., cos( angle ), -sin( angle ) ],\
                                              [ 0., sin( angle ),  cos( angle ) ] ] )
    return rot1_matrix

def rot2(input):
    # rotation matrix about axis 2
    angle = input
    rot2_matrix = np.array( [ [ cos( angle ), 0., -sin( angle ) ], [ 0., 1., 0. ],\
                              [ sin( angle ), 0.,  cos( angle ) ] ] )
    return rot2_matrix

def rot3(input):
    # rotation matrix about axis 3
    angle = input
    rot3_matrix = np.array( [ [  cos( angle ), sin( angle ), 0. ],\
                              [ -sin( angle ), cos( angle ), 0. ],\
                                                     [ 0., 0., 1. ] ] )
    return rot3_matrix

################################################################################
#                  E X P O R T     F U N C T I O N S:
################################################################################

def fncH(input):
    # pull out the inputs for the generation of estimated observation data
    state      = input[0]
    SPICE_data = input[1]
    angles     = input[2]
    extras     = input[-1]

    # number of beacon observations
    n_beacons  = len(extras['obs_beacons'])

    # focal length
    FoL = extras['FoL']

    # camera and P&L parameters
    resolution = extras['resolution']
    p0         = resolution[0] / 2.
    l0         = resolution[1] / 2.
    Kx         = 1. / extras['pixel_width']
    Ky         = 1. / extras['pixel_height']
    Dx         = extras['pixel_direction']
    Dy         = extras['line_direction']

    # count the number of QoIs
    n_state = state.shape[1]

    # initiate the H matrix
    H = np.zeros((2*n_beacons, n_state))

    # loop through beacons
    for ii in xrange(n_beacons):
        beacon_state = SPICE_data[ii,:]
        # calculate the difference between the positions of the beacon and state
        r_diff = beacon_state[0:3] - state[ii, 0:3]

        # calculate norm of pointing vector
        rho = np.linalg.norm( r_diff )

        # create inertial unit pointing vector (A_hat_I)
        A_hat_I = r_diff / rho

        # partials of A_hat_I with respect to position components
        dA_IdX = np.array([ r_diff[0]**2 / rho**3 - 1. / rho,\
                            r_diff[0] * r_diff[1] / rho**3,\
                            r_diff[0] * r_diff[2] / rho**3 ])

        dA_IdY = np.array([ r_diff[1] * r_diff[0] / rho**3,\
                            r_diff[1]**2 / rho**3 - 1. / rho,\
                            r_diff[1] * r_diff[2] / rho**3 ])

        dA_IdZ = np.array([ r_diff[2] * r_diff[0] / rho**3,\
                            r_diff[2] * r_diff[1] / rho**3,\
                            r_diff[2]**2 / rho**3 - 1. / rho ])

        # partials of A_hat_TV with respect to position components
        # Compute DCM from the Inertia frame to camera frame:
        DCM_TVI = np.dot(extras['DCM_TVB'], extras['DCM_BI'])
        
        dA_TVdX = np.dot( DCM_TVI, dA_IdX )
        dA_TVdY = np.dot( DCM_TVI, dA_IdY )
        dA_TVdZ = np.dot( DCM_TVI, dA_IdZ )

        # calculate camera frame unit pointing vector (A_hat_TV)
        A_hat_TV = np.dot( DCM_TVI, A_hat_I )

        # partials of millimeter frame with respect to state
        dMMdX = FoL * np.array( [ -1. / A_hat_TV[2]**2 * dA_TVdX[2] * A_hat_TV[0] +\
                                               dA_TVdX[0] / A_hat_TV[2],\
                        -1. / A_hat_TV[2]**2 * dA_TVdX[2] * A_hat_TV[1] +\
                                               dA_TVdX[1] / A_hat_TV[2]  ] )

        dMMdY = FoL * np.array( [ -1. / A_hat_TV[2]**2 * dA_TVdY[2] * A_hat_TV[0] +\
                                               dA_TVdY[0] / A_hat_TV[2],\
                        -1. / A_hat_TV[2]**2 * dA_TVdY[2] * A_hat_TV[1] +\
                                               dA_TVdY[1] / A_hat_TV[2]  ] )

        dMMdZ = FoL * np.array( [ -1. / A_hat_TV[2]**2 * dA_TVdZ[2] * A_hat_TV[0] +\
                                               dA_TVdZ[0] / A_hat_TV[2],\
                        -1. / A_hat_TV[2]**2 * dA_TVdZ[2] * A_hat_TV[1] +\
                                               dA_TVdZ[1] / A_hat_TV[2]  ] )  
  
        # partials of pixel with respect to state
        H[2*ii, :]     = Kx * Dx * np.array([ dMMdX[0], dMMdY[0], dMMdZ[0], 0, 0, 0])
        
        # partials of line with respect to state
        H[2*ii + 1, :] = Ky * Dy * np.array([ dMMdX[1], dMMdY[1], dMMdZ[1], 0, 0, 0])

    return H

def fncG(input):
    # pull out the inputs for the generation of estimated observation data
    state      = input[0]
    SPICE_data = input[1]
    extras     = input[-1]

    # number of beacon observations
    n_beacons  = len(extras['beacons'])

    # focal length
    FoL = extras['FoL']

    # camera and P&L parameters
    resolution = extras['resolution']
    p0         = resolution[0] / 2.
    l0         = resolution[1] / 2.
    Kx         = 1. / extras['pixel_width']
    Ky         = 1. / extras['pixel_height']
    Dx         = extras['pixel_direction']
    Dy         = extras['line_direction']

    # create array to be twice the size of the number of beacons. this is because there are
    # two data types: pixel and line
    G = np.zeros([n_beacons,2])
    for ii in range(n_beacons):

        beacon_state = SPICE_data[ii,:]
        # calculate the difference between the positions of the beacon and state
        r_diff = beacon_state[0:3] - state[ii, 0:3]

        # create inertial unit pointing vector (A_hat)
        Ahat_I = r_diff / np.linalg.norm( r_diff )

        # Compute DCM from the Inertia frame to camera frame:
        DCM_TVI = np.dot(extras['DCM_TVB'], extras['DCM_BI'])

        # rotate inertial pointing vector to camera (TV) frame
        Ahat_TV = np.dot( DCM_TVI, Ahat_I )

        # convert TV frame to millimeter units
        x_mm = FoL / Ahat_TV[2] * Ahat_TV[0]
        y_mm = FoL / Ahat_TV[2] * Ahat_TV[1]

        # convert millimeter to P&L
        pixel = Kx * Dx * x_mm + p0
        line  = Ky * Dy * y_mm + l0

        # add pixel and line to output array
        G[ii,0] = pixel
        G[ii,1] = line
    return G


###############################################################################
#                       M A I N     F U N C T I O N:
################################################################################

def main():
    return


if __name__ == "__main__":
    main()
