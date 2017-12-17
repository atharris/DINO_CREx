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


## \defgroup pixel_and_line pixelLineBatch - pixel and line measurements generation
##   @{
## The module for the creation of pixel and line data.
#
# Overview {#overview}
# ====
#
# Purpose
# -----
# This script contains two exportable functions as well as four internal functions. The primary goal of the module is to contain a reference observation function and observation-state mapping matrix function that take inputs and outputs as dictated by this software suite's batch functions. 
#
# Contents
# -----
# The following exportable functions are contained in this module:
#
# - `fncH`
# - `fncG`
#
# `fncH` is the function that computes the method for computing the observation-state mapping matrix, while `fncG` computes the observed quantities for the respective inputs. This language is chosen for the latter due to the fact that although the intended use of `fncG.py` is to provide reference measurements for the batch filter, it is possible that it could be used to generate a data set representing "true" measurements. It would then be up to the operator to provide appropriate inputs for either case.
#
# The following internal functions are contained in this module:
#
# - `I_to_TC`
# - `rot1`
# - `rot2`
# - `rot3`
#
# All four functions are linear transformations, with `I_to_TV` calling the remaining three rotation functions. Specifically, `I_to_TV` creates an inertial to camera coordinate transformation matrix, while the other functions compute rotation matrices for a 3-D coordinate system. 
#
# As with other modules in the state estimation nav filter, there is a reliance on the `extras` dictionary to pass through parameters to various functions. It is noted that the `extras` dictionary should never be an output from a function.
#
# The Code
# =====
#
# `fncH`
# -----
# `fncH` is a function that produces the observation-state mapping matrix (H) used within a batch filter. The following is a table of inputs and associated brief descriptions:
#
# Name      | Description                                           | Size/Type
# -----     | -------------------                                   | -----
# state | input state of the quantities of interest for each observation time | (N,d) numpy array  
# beaconStates | propagated beacon states | (N,6) numpy array
# angles | reference attitude angles for each observation time| (N,3) numpy array
# extras    | dictionary of various parameters                      | dictionary
#
# The same is repeated for outputs:
#
# Name      | Description         | Size/Type
# -----     | ------------------- | -----
# H | observation-state mapping matrices for each time step | (2*N,6) numpy array
#
# The size of the output is determined from the fact that there are two rows of partials at each time. This is due to the fact that a pixel and line filter has two observation types. 
#
# The significant lines of code in this function are within a loop that counts the number of observations. At each iteration, the relative position vector between the selected beacon and spacecraft is calculated
# ~~~~~~~~~~~~~~~~{.py}
#        # calculate the difference between the positions of the beacon and state
#        positionDiff = selectedBeacon[0:3] - state[ii, 0:3]
# ~~~~~~~~~~~~~~~~
#
# This `positionDiff` is then used to calculate a unit pointing vector and its derivatives.
# ~~~~~~~~~~~~~~~~{.py}
#        # create inertial unit pointing vector (A_hat_I)
#        A_hat_I = positionDiff / rho
#
#        # partials of A_hat_I with respect to position components
#        dA_IdX = np.array([ positionDiff[0]**2 / rho**3 - 1. / rho,\
#                            positionDiff[0] * positionDiff[1] / rho**3,\
#                            positionDiff[0] * positionDiff[2] / rho**3 ])
#
#        dA_IdY = np.array([ positionDiff[1] * positionDiff[0] / rho**3,\
#                            positionDiff[1]**2 / rho**3 - 1. / rho,\
#                            positionDiff[1] * positionDiff[2] / rho**3 ])
#
#        dA_IdZ = np.array([ positionDiff[2] * positionDiff[0] / rho**3,\
#                            positionDiff[2] * positionDiff[1] / rho**3,\
#                            positionDiff[2]**2 / rho**3 - 1. / rho ])
# ~~~~~~~~~~~~~~~~
#
# These values are transformed into the millimeter coordinates
# ~~~~~~~~~~~~~~~~{.py}
#        # partials of millimeter frame with respect to state
#        dMMdX = FoL * np.array( [ -1. / A_hat_TV[2]**2 * dA_TVdX[2] * A_hat_TV[0] +\
#                                               dA_TVdX[0] / A_hat_TV[2],\
#                        -1. / A_hat_TV[2]**2 * dA_TVdX[2] * A_hat_TV[1] +\
#                                               dA_TVdX[1] / A_hat_TV[2]  ] )
#
#        dMMdY = FoL * np.array( [ -1. / A_hat_TV[2]**2 * dA_TVdY[2] * A_hat_TV[0] +\
#                                               dA_TVdY[0] / A_hat_TV[2],\
#                        -1. / A_hat_TV[2]**2 * dA_TVdY[2] * A_hat_TV[1] +\
#                                               dA_TVdY[1] / A_hat_TV[2]  ] )
#
#        dMMdZ = FoL * np.array( [ -1. / A_hat_TV[2]**2 * dA_TVdZ[2] * A_hat_TV[0] +\
#                                               dA_TVdZ[0] / A_hat_TV[2],\
#                        -1. / A_hat_TV[2]**2 * dA_TVdZ[2] * A_hat_TV[1] +\
#                                               dA_TVdZ[1] / A_hat_TV[2]  ] ) 
# ~~~~~~~~~~~~~~~~
# The resulting values are a then a constant multiplication away from being the pointing vector derivatives in the camera frame.
#
# `fncG` 
# -----
# `fncG` is a function that computes measurement data for each observation time. It relies on inputs similar to `fncH`. In the case of this module, the measurement data is pixel and line. Specifically:
#
# Name      | Description                                           | Size/Type
# -----     | -------------------                                   | -----
# state | input state of the quantities of interest for each observation time | (N,d) numpy array  
# beaconStates | propagated beacon states | (N,6) numpy array
# angles | reference attitude angles for each observation time | (N,3) numpy array
# extras    | dictionary of various parameters | dictionary
#
# The same is repeated for outputs:
#
# Name      | Description         | Size/Type
# -----     | ------------------- | -----
# G | measurement data for each time step | (N,2) numpy array
#
# Due to the nature of the relationship between the H matrix and measurement generation, the code of `fncG` largely mirrors that of `fncH`. There is a significant when considering the fact that `fncH` effectively computes the derivatives of `fncG`. The measurement generation begins with computing the relative position `positionDiff`. For this function, however, the derivative of the pointing vector is never calculated. Therefore, the code consists of the calculation of the unit vector
# ~~~~~~~~~~~~~~~~{.py}
#        # create inertial unit pointing vector (A_hat)
#        Ahat_I = positionDiff / np.linalg.norm( positionDiff )
# ~~~~~~~~~~~~~~~~
#
# Followed by a rotation from inertial to camera and a conversion to millimeter coordinates
# ~~~~~~~~~~~~~~~~{.py}
#        # rotate inertial pointing vector to camera (TV) frame
#        Ahat_TV = np.dot( DCM_TVI, Ahat_I )
#
#        # convert TV frame to millimeter units
#        x_mm = FoL / Ahat_TV[2] * Ahat_TV[0]
#        y_mm = FoL / Ahat_TV[2] * Ahat_TV[1]
# ~~~~~~~~~~~~~~~~
#
# Millimeter is then converted to pixel and line by adding the resolution offset and multiplying by pixel width (for pixel units) and pixel height (for line units).
# ~~~~~~~~~~~~~~~~{.py}
#        # convert millimeter to P&L
#        pixel = Kx * Dx * x_mm + p0
#        line  = Ky * Dy * y_mm + l0
# ~~~~~~~~~~~~~~~~
# The pixel and line measurements are then stored in the G array for output
#
# `rot#` 
# -----
# The rotation matrices included in this module are typical 3-D coordinate rotation matrices in the 1-2-3 axis format. Therefore, the inputs and outputs are functionally similar for all three functions. The inputs are the appropriate attitude angle for the stated rotation direction:
#
# Name      | Description                                           | Size/Type
# -----     | -------------------                                   | -----
# angle | rotation angle, in radians, to be applied to the matrix formula | double
#
# The outputs are all numpy arrays calculated using the relevant rotation matrix equations:
#
# Name      | Description         | Size/Type
# -----     | ------------------- | -----
# matrixRot# | resulting rotation matrix | (3,3) numpy array
#
# These matrices are called in the intertia to camera coordinate transformation matrix function
#
# `I_to_TV` 
# -----
# The `I_to_TV` function constructs a transformation matrix for interial to the camera coordinate frame. This differs from the `rot#` functions as the input is an array of all three angles, which are then used as inputs to each rotation function individually:
#
# Name      | Description                                           | Size/Type
# -----     | -------------------                                   | -----
# angles | rotation angles, in radians | (3,) numpy array
#
# The output is an array that functions as a transformation matrix:
#
# Name      | Description         | Size/Type
# -----     | ------------------- | -----
# I_to_TV_matrix | inertial to camera transformation matrix | (3,3) numpy array
#
# This matrix is composed of a 1-2-3 rotation matrix application, as seen in the code:
# ~~~~~~~~~~~~~~~~{.py}
#     I_to_TV_matrix = np.dot( rot3( angles[0] ), \
#                     np.dot( rot2( angles[1] ), rot1( angles[2] ) ) )
# ~~~~~~~~~~~~~~~~
#
## @}

################################################################################
#                  I N T E R N A L    F U N C T I O N S:
################################################################################

def I_to_TV(input):

	# pull out the angles from the input (RA, dec, twist) in that order
    angles = input
    I_to_TV_matrix = np.dot( rot3( angles[0] ), \
                     np.dot( rot2( angles[1] ), rot1( angles[2] ) ) )

    return I_to_TV_matrix

def rot1(input):
    # rotation matrix about axis 1
    angle = input
    matrixRot1 = np.array( [ [ 1., 0., 0. ], [ 0., cos( angle ), -sin( angle ) ],\
                                              [ 0., sin( angle ),  cos( angle ) ] ] )
    return matrixRot1

def rot2(input):
    # rotation matrix about axis 2
    angle = input
    matrixRot2 = np.array( [ [ cos( angle ), 0., -sin( angle ) ], [ 0., 1., 0. ],\
                              [ sin( angle ), 0.,  cos( angle ) ] ] )
    return matrixRot2

def rot3(input):
    # rotation matrix about axis 3
    angle = input
    matrixRot3 = np.array( [ [  cos( angle ), sin( angle ), 0. ],\
                              [ -sin( angle ), cos( angle ), 0. ],\
                                                     [ 0., 0., 1. ] ] )
    return matrixRot3

################################################################################
#                  E X P O R T     F U N C T I O N S:
################################################################################

def fncH(input):
    # pull out the inputs for the generation of estimated observation data
    state        = input[0]
    beaconStates = input[1]
    angles       = input[2]
    extras       = input[-1]

    # number of beacon observations
    nObservations  = len(extras['obs_beacons'])

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
    stateDimension = state.shape[1]

    # initiate the H matrix
    H = np.zeros((2*nObservations, stateDimension))

    # loop through beacons
    for ii in xrange(nObservations):
        selectedBeacon = beaconStates[ii,:]
        # calculate the difference between the positions of the beacon and state
        positionDiff = selectedBeacon[0:3] - state[ii, 0:3]

        # calculate norm of pointing vector
        rho = np.linalg.norm( positionDiff )

        # create inertial unit pointing vector (A_hat_I)
        A_hat_I = positionDiff / rho

        # partials of A_hat_I with respect to position components
        dA_IdX = np.array([ positionDiff[0]**2 / rho**3 - 1. / rho,\
                            positionDiff[0] * positionDiff[1] / rho**3,\
                            positionDiff[0] * positionDiff[2] / rho**3 ])

        dA_IdY = np.array([ positionDiff[1] * positionDiff[0] / rho**3,\
                            positionDiff[1]**2 / rho**3 - 1. / rho,\
                            positionDiff[1] * positionDiff[2] / rho**3 ])

        dA_IdZ = np.array([ positionDiff[2] * positionDiff[0] / rho**3,\
                            positionDiff[2] * positionDiff[1] / rho**3,\
                            positionDiff[2]**2 / rho**3 - 1. / rho ])

        # partials of A_hat_TV with respect to position components
        # Compute DCM from the Inertia frame to camera frame:
        if len( angles ) == 0:
          DCM_TVI = np.dot(extras['DCM_TVB'], extras['DCM_BI'])
        else:
          DCM_TVI = I_to_TV( angles ) 

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
        H[2*ii, 0:3]     = Kx * Dx * np.array([ dMMdX[0], dMMdY[0], dMMdZ[0] ])
        
        # partials of line with respect to state
        H[2*ii + 1, 0:3] = Ky * Dy * np.array([ dMMdX[1], dMMdY[1], dMMdZ[1] ])

    return H

def fncG(input):
    # pull out the inputs for the generation of estimated observation data
    state      = input[0]
    beaconStates = input[1]
    angles     = input[2]
    extras     = input[-1]

    # number of beacon observations
    nObservations  = len(extras['obs_beacons'])

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
    G = np.zeros([nObservations,2])

    for ii in range(nObservations):

        selectedBeacon = beaconStates[ii,:]
        # calculate the difference between the positions of the beacon and state
        positionDiff = selectedBeacon[0:3] - state[ii, 0:3]

        # create inertial unit pointing vector (A_hat)
        Ahat_I = positionDiff / np.linalg.norm( positionDiff )
        
        # Compute DCM from the Inertia frame to camera frame:
        if len( angles ) == 0:
          DCM_TVI = np.dot(extras['DCM_TVB'], extras['DCM_BI'])
        else:
          DCM_TVI = I_to_TV( angles ) 

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
