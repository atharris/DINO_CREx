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
##############################################################################
# Log path in order to get pyswice from BSK

import sys, os, inspect

filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))
path2 = os.path.dirname(os.path.abspath(filename))
bskName = 'Basilisk'
dinoName = 'DINO_CREx'
splitPath = path.split(dinoName)
dinoSpicePath = splitPath[0] + dinoName + '/DINObatch/SPICE/'
bskSpicePath = splitPath[0] + bskName + '/External/EphemerisData/'
bskPath = splitPath[0] + bskName + '/'
sys.path.append(bskPath + 'modules')
sys.path.append(bskPath + 'PythonModules')
sys.path.append(dinoSpicePath)


import numpy as np
try:
    import pyswice
except ImportError:
    from Basilisk import pyswice

import pdb

## \defgroup EOMs_vanilla posVel - vanilla EOMs of position and velocity
##   @{
## The module for the EOMs and A matrix for position and velocity.
#
# Overview {#overview}
# ====
#
# Purpose
# -----
# This script contains two exportable functions. This module hosts the equations of motion (EOMs) to be called by a desired propagator, as well as the code necessary to formulate a relevant A matrix. This matrix is essential to propagating deviations under linear assumptions and is utilized in the EOMs for the propagation of the STM. 
#
# Contents
# -----
# The following exportable functions are contained in this module:
#
# - `matrixA`
# - `EOM`
#
# The `matrixA` function is used for the calculation of the A matrix, which is composed of various derivatives of the quantities of interest. This matrix is further explained in this documentation. The `EOM` function ultimately calls the `matrixA` function, and calculates the time derivatives of the quantities of interest for the filter propagator. `EOM` contains a number of accelerations such as that from gravity or SRP. 
#
# The Code
# =====
#
# `matrixA`
# -----
# The A matrix is the first order derivative of the state derived from linear assumptions. Because of this, its formulation relies on knowing the time derivatives of the quantities of interest, and the derivatives of these with respect to the quantities of interest themselves. The mathematical formulation of this matrix is found in Eq. 4.2.6 of Tapley, Schutz and Born, and it is illustrated in Example 4.2.1. A table of required input is provided below:
#
# Name      | Description                                           | Size/Type
# -----     | -------------------                                   | -----
# spacecraftPosition | position of spacecraft. original input of full state | (d,) numpy array  
# nSecondaries | number of secondary gravitational bodies | int
# muPrimary | gravitational parameter of primary gravitational body | float
# muSecondaries | list of secondary body gravitional parameters | list
# kSRP | solar radiation pressure constant | float
# cR | spacecraft coefficient of reflectivity | float
# sunPosition | position of the sun | (3,) numpy array
# position_secondaries_primary | relative positions of secondary bodies with respect to primary | (3, nSecondaries) numpy array
#
# The output of this function is the A matrix:
#
# Name      | Description         | Size/Type
# -----     | ------------------- | -----
# A | matrix derived from first order linear assumptions | (d,d) numpy array
#
# This function runs off of the assumption that the filter quantities of interest are limited to position and velocity. Therefore, the derivations needed involve taking derivatives of acceleration with respect to position, as the EOMs are invariant with respect to velocity, and the velocity is invariant with respect to position. A figure is provided to illustrate the compositon of the A matrix in this case.
#
# \image html matrixA.svg
#
# Here, we consider each component of the acceleration derivatives. It is noted that each derivative is solved for in a vector format, so the accelerations are able to be summed up with little fanfare at the end. 
#
# First, the acceleration due to the primary body gravity is taken with respect to position:
# ~~~~~~~~~~~~~~~~{.py}
#    # the spacecraftPosition derivate associated with the primary gravitational force
#    dFdR_p = -muPrimary * ( np.identity(3) / np.linalg.norm(spacecraftPosition) ** 3 -\
#    3 * np.dot(spacecraftPosition, spacecraftPosition.T) / \
#    np.linalg.norm(spacecraftPosition) ** 5 )
# ~~~~~~~~~~~~~~~~
#
# This is followed up by looping through the list of secondary bodies and calculating the position derivatives of the contributed gravitational accelerations:
# ~~~~~~~~~~~~~~~~{.py}
#    # loop through the secondary bodies
#    for ii in xrange(nSecondaries) :
#        positionSecondary = np.expand_dims( position_secondaries_primary[:, ii], axis = 1 )
#        dFdR_s += -muSecondaries[ii] * (
#        np.identity(3) / np.linalg.norm(spacecraftPosition - positionSecondary) ** 3 -\
#        3 * np.dot(spacecraftPosition - positionSecondary,(spacecraftPosition -\
#        positionSecondary).T)/np.linalg.norm(spacecraftPosition - positionSecondary) ** 5)
# ~~~~~~~~~~~~~~~~
#
# Lastly, we address the derivatives of the SRP acceleration and sum up all accleration derivatives:
# ~~~~~~~~~~~~~~~~{.py}
#    # the spacecraftPosition derivative associated with the SRP force
#    dFdR_SRP = cR * kSRP * (np.identity(3) / \
#       np.linalg.norm(spacecraftPosition - sunPosition) ** 3 -
#       3 * np.dot(spacecraftPosition - sunPosition, (spacecraftPosition - sunPosition).T) /\
#                            np.linalg.norm(spacecraftPosition - sunPosition) ** 5)
#    # total spacecraftPosition derivatives of forces
#    dFdR = dFdR_p + dFdR_s + dFdR_SRP
# ~~~~~~~~~~~~~~~~
#
# The A matrix is then populated with these values, as well as 1s for velocity derivatives, and is returned to the function that called it.
#
#
# `EOM` 
# -----
# `EOM` is a function that computes the equations of motion when given appropriate inputs including the state, timestep and various force parameters. The list of inputs is long, but all are necessary for the function to operate properly
#
# Name      | Description                                           | Size/Type
# -----     | -------------------                                   | -----
# state | input state of the quantities of interest and STM for the considered time step| (d(1+d),) numpy array  
# et | time step for propagator | float
# primaryIndex| index of the primary body in the bodies list | int
# secondaryIndices| indices of the seoncdary bodies in the bodies list | list of int
# nSecondaries| number of secondary bodies | int
# muPrimary| gravitational parameter of primary body | float
# muSecondaries| gravitational parameters of secondary bodies | list of floats
# kSRP| solar radiation pressure coefficient | float
# cR| spacecraft coefficient of reflectivity | float
# abcorr| aberration correction for SPICE (not considered as of 11/17) | float
# refFrame| reference frame for SPICE | str
# bodies| list of gravitational bodies for SPICE | list of str
# stateDimension | dimension of the filter quantities of interest | int
#
# The same is repeated for outputs:
#
# Name      | Description         | Size/Type
# -----     | ------------------- | -----
# dState | quantity of interest and STM derivatives to be integrated | list of float
#
# As with many EOM functions, this code is largely devoted to calculating and arranging the desired accelerations in the scenario. This function does contain interesting behavior when considering the STM and the gravitational bodies. The STM is resized for input so that the `state` variable is not a matrix, but a vector. The position of gravitational bodies in the considered reference frame (`refFrame`) a found via the SPICE tool. We first see this when calculating the gravitational force contributed by the secondary bodies. Before this, however, we consider the gravitational acceleration of the primary body:
# ~~~~~~~~~~~~~~~~{.py}
#    # gravitational force from primary body
#    fPrimary = -muPrimary * state[0:3] / np.linalg.norm(state[0:3]) ** 3
# ~~~~~~~~~~~~~~~~
#
# Here, we note that there is an assumption that the input `state` is the position of the spacecraft with respect to the sun. This is true for all data that is generated via the `data_generation` function. 
#
# Following this, we have a loop over the number of secondary bodies. The accelerations due to these bodies are linearly summed. The position of the secondary body with respect to the primary body is found first:
# ~~~~~~~~~~~~~~~~{.py}
#        # determine distance from secondary to primary body
#        positionArray = np.zeros(3)
#        stateSpice = pyswice.new_doubleArray(6)
#        lt = pyswice.new_doubleArray(1)
#        pyswice.spkezr_c(bodies[secondaryIndices[ii]], et, refFrame,
#                        abcorr, bodies[primaryIndex], stateSpice, lt)
#        for i in range(3):
#            positionArray[i] = pyswice.doubleArray_getitem(stateSpice, i)
#        position_secondaries_primary[:, ii] = positionArray
# ~~~~~~~~~~~~~~~~
#
# We note that the `pyswice.spkezr_c()` line contains the indexing for the secondary body `bodies[secondaryIndices[ii]]`, as well as centering the position vector on the primary body by using the corresponding index `primaryIndex`. Using this relative position, we are then able to calculate the acceleration contribution.
# ~~~~~~~~~~~~~~~~{.py}
#        # calculate the "third body" force
#        f3rdBodies += -muSecondaries[ii] * \
#           ( ( state[0:3] - position_secondaries_primary[:, ii] ) / \
#             np.linalg.norm(state[0:3] - position_secondaries_primary[:, ii] ) ** 3 + \
#             position_secondaries_primary[:, ii] / \
#             np.linalg.norm( position_secondaries_primary[:, ii] )**3 )
# ~~~~~~~~~~~~~~~~
#
# The final acceleration in this function is that of solar radiation pressure (SRP). This is a function of the distance from the sun to the spacecraft, i.e.,
# ~~~~~~~~~~~~~~~~{.py}
#    # SRP force
#    fSRP = cR * kSRP * state[0:3] / np.linalg.norm(state[0:3]) ** 3
# ~~~~~~~~~~~~~~~~
# 
# The acceleration vectors are then able to be summed. Various force sources can be added to this basic template. This is seen in the \ref pos_vel_acc `posVelAcc.py` function. 
#
# Finally, the derivative of the STM is calculate by computing the A matrix and multiplying it with the most recent STM. This is found in Eq. 4.2.10 and 4.2.24 of Tapley, Schutz and Born, and is accomplished via the following lines of code
# ~~~~~~~~~~~~~~~~{.py}
#    # A matrix calculation
#    A = matrixA(args)
#
#    # compute the derivative of the STM
#    dPhi = np.dot(A, phi)
# ~~~~~~~~~~~~~~~~
#
# The STM is then reshaped to a vector and appended to the derivative of the state. 
## @}


################################################################################
#                  E X P O R T E D     F U N C T I O N S
################################################################################

# -------------------------------------------------------------------------------

def matrixA(input):
    stateDimension     = len(input[0])
    spacecraftPosition = np.expand_dims(input[0][0:3],axis=1)
    nSecondaries       = input[1]
    muPrimary          = input[2]
    muSecondaries      = input[3]
    kSRP               = input[4]
    cR                 = input[5]
    position_secondaries_primary = input[-1]

    # set the size for the A matrix and premptively populate with zeros
    A = np.zeros((stateDimension, stateDimension))

    # the spacecraftPosition derivate associated with the primary gravitational force
    dFdR_p = -muPrimary * ( np.identity(3) / np.linalg.norm(spacecraftPosition) ** 3 -\
    3 * np.dot(spacecraftPosition, spacecraftPosition.T) / \
    np.linalg.norm(spacecraftPosition) ** 5 )

    # the spacecraftPosition derivatives associated 
    # with gravitational force from secondary bodies
    dFdR_s = np.zeros( ( 3, 3 ) )

    # loop through the secondary bodies
    for ii in xrange(nSecondaries) :
        positionSecondary = np.expand_dims( position_secondaries_primary[:, ii], axis = 1 )
        dFdR_s += -muSecondaries[ii] * (
        np.identity(3) / np.linalg.norm(spacecraftPosition - positionSecondary) ** 3 -\
        3 * np.dot(spacecraftPosition - positionSecondary,(spacecraftPosition -\
        positionSecondary).T)/np.linalg.norm(spacecraftPosition - positionSecondary) ** 5)

    # the spacecraftPosition derivative associated with the SRP force
    dFdR_SRP = cR * kSRP * (np.identity(3) / \
       np.linalg.norm(spacecraftPosition) ** 3 -
       3 * np.dot(spacecraftPosition, spacecraftPosition.T) /\
                            np.linalg.norm(spacecraftPosition) ** 5)

    # total spacecraftPosition derivatives of forces
    dFdR = dFdR_p + dFdR_s + dFdR_SRP

    # populate the A matrix, where each row (D) is a time derivative of a QoI
    # and each column (E) is a derivative taken with respect to each QoI
    # [D,E]
    #       X | Y | Z | dX | dY | dZ |
    # dX
    # --
    # dY
    # --
    # dZ
    # --
    # d(dX)
    # --
    # d(dY)
    # --
    # d(dZ)
    # --

    A[0, 3] = 1
    A[1, 4] = 1
    A[2, 5] = 1
    A[3:6, 0:3] = dFdR

    return A


def EOM(state, et, primaryIndex, secondaryIndices, nSecondaries, muPrimary, muSecondaries,
        kSRP, cR, abcorr, refFrame, bodies, stateDimension):

    # pull out the STM
    phi = np.array(state[stateDimension:], copy=True)
    phi = np.reshape(phi, (stateDimension, stateDimension) )

    # gravitational force from primary body
    fPrimary = -muPrimary * state[0:3] / np.linalg.norm(state[0:3]) ** 3

    # gravitational force from secondary bodies
    f3rdBodies = 0

    # set the size of the spacecraftPositionitionSecondaries between 
    # secondary bodies and primary body
    position_secondaries_primary = np.zeros((3, nSecondaries))

    # loop through the secondary bodies
    for ii in range(nSecondaries):
        # determine distance from secondary to primary body
        positionArray = np.zeros(3)
        stateSpice = pyswice.new_doubleArray(6)
        lt = pyswice.new_doubleArray(1)
        pyswice.spkezr_c(bodies[secondaryIndices[ii]], et, refFrame,
                        abcorr, bodies[primaryIndex], stateSpice, lt)
        for i in range(3):
            positionArray[i] = pyswice.doubleArray_getitem(stateSpice, i)
        position_secondaries_primary[:, ii] = positionArray

        # calculate the "third body" force
        f3rdBodies += -muSecondaries[ii] * \
           ( ( state[0:3] - position_secondaries_primary[:, ii] ) / \
             np.linalg.norm(state[0:3] - position_secondaries_primary[:, ii] ) ** 3 + \
             position_secondaries_primary[:, ii] / \
             np.linalg.norm( position_secondaries_primary[:, ii] )**3 )


    # SRP force
    fSRP = cR * kSRP * state[0:3] / \
            np.linalg.norm(state[0:3]) ** 3

    # total force (acceleration) vector
    f = fPrimary + f3rdBodies + fSRP

    # args for the A matrix function
    args = (state[0:stateDimension], nSecondaries, muPrimary, muSecondaries, kSRP, cR,\
            position_secondaries_primary)

    # A matrix calculation
    A = matrixA(args)

    # compute the derivative of the STM
    dPhi = np.dot(A, phi)
    dPhi = np.reshape(dPhi,stateDimension * stateDimension)

    # acceleration vector to be returned to the integrator
    dState = [state[3], state[4], state[5], f[0], f[1], f[2]]
    dState += list(dPhi)

    return dState


################################################################################
#                       M A I N     F U N C T I O N:
################################################################################

def main():
    return


if __name__ == "__main__":
    main()
