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
    bskSpicePath = splitPath[0] + bskName + '/supportData/EphemerisData/'
import pdb

## \defgroup EOMs_acc posVelAcc - EOMs for unmodeled acceleration
##   @{
## The module for the EOMs and A matrix for position, velocity, and constant accelerations.
#
# Overview {#overview}
# ====
#
# Purpose
# -----
# This script contains two exportable functions. This module hosts the equations of motion (EOMs) to be called by a desired propagator, as well as the code necessary to formulate a relevant A matrix. This matrix is essential to propagating deviations under linear assumptions and is utilized in the EOMs for the propagation of the STM. 
#
# This code is largely similar to that of \ref EOMs_vanilla "`posVel.py`". There are exceptions, however. All of these center around the inclusion of constant accelerations in the x, y and z directions. This affects lines of code that will be identified in this document. For some theory and a verbose explanation of the EOMs, it is suggested that the reader look at the documentation of \ref EOMs_vanilla "`posVel.py`".
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
# This function runs off of the assumption that the filter quantities of interest are expandsed from only position and velocity to include unmodeled, constant accelerations. A figure is provided to illustrate the compositon of the A matrix in this case.
#
# \image html matrixA_acc.svg
#
# Due to the constant nature of the estimatable accelerations, the A matrix is largely unchanged. It is expanded, however, to be a (9x9) matrix rather than (6x6). An identity matrix occupies a (3x3) section of the matrix devoted to the derivatives of acceleration with respect to the constant, unmodeled accelerations. The addition of these 1s indicates the constant nature of the unmodeled accelerations.
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
# The code of `EOM` is almost identical to that of the `EOM` function found in \ref EOMs_vanilla "`posVel.py`". There are two exceptions to this. The first is the addition of the last three values of the quantities of interest to the force summation
# ~~~~~~~~~~~~~~~~{.py}
#     # total force (acceleration) vector
#     f = fPrimary + f3rdBodies + fSRP + state[6:9]
# ~~~~~~~~~~~~~~~~
#
# The addition of this code accounts for the constant accelerations in the EOMs. The other alteration is the inclusion of extra zeros appended to the derivative vector. This is seen in the line
# ~~~~~~~~~~~~~~~~{.py}
#     # acceleration vector to be returned to the integrator
#     dState      = [0] * stateDimension 
# ~~~~~~~~~~~~~~~~
# Since `stateDimension` = 9, the last three entries of `dState` will be zeros, and therefore indicative of three constant values. These values are the constant, unmodeled accelerations. 
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
    dFdR_p = -muPrimary * (
    np.identity(3) / np.linalg.norm(spacecraftPosition) ** 3 -
    3 * np.dot(spacecraftPosition, spacecraftPosition.T) / np.linalg.norm(spacecraftPosition) ** 5 
    )

    # the spacecraftPosition derivatives associated with gravitational force from secondary bodies
    dFdR_s = np.zeros( ( 3, 3 ) )

    # loop through the secondary bodies
    for ii in xrange(nSecondaries) :
        positionSecondary = np.expand_dims( position_secondaries_primary[:, ii], axis = 1 )
        dFdR_s += -muSecondaries[ii] * (
        np.identity(3) / np.linalg.norm(spacecraftPosition - positionSecondary) ** 3 -
        3 * np.dot(spacecraftPosition - positionSecondary,(spacecraftPosition - positionSecondary).T)/
        np.linalg.norm(spacecraftPosition - positionSecondary) ** 5
        )

    # the spacecraftPosition derivative associated with the SRP force
    dFdR_SRP = cR * kSRP * (np.identity(3) / np.linalg.norm(spacecraftPosition) ** 3 -
                            3 * np.dot(spacecraftPosition, spacecraftPosition.T) /
                            np.linalg.norm(spacecraftPosition) ** 5)

    # total spacecraftPosition derivatives of forces
    dFdR = dFdR_p + dFdR_s + dFdR_SRP

    # populate the A matrix, where each row (D) is a time derivative of a QoI.
    # and each column (E) is the QoI that the A matrix derivative is taken wrt.
    # unmodeled acceleration terms (ax, ay, az) are constant wrt time.
    # [D,E]
    #       X | Y | Z | dX | dY | dZ | ax | ay | az |
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
    # d(ax)
    # --
    # d(ay)
    # --
    # d(az)
    # --

    A[0, 3] = 1
    A[1, 4] = 1
    A[2, 5] = 1
    A[3, 6] = 1
    A[4, 7] = 1
    A[5, 8] = 1
    A[3:6, 0:3] = dFdR

    return A


def EOM(state, et, primary_index, secondary_indices, nSecondaries, muPrimary, muSecondaries,
        kSRP, cR, abcorr, ref_frame, bodies, stateDimension):

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
        pyswice.spkezr_c(bodies[secondary_indices[ii]], et, ref_frame,
                        abcorr, bodies[primary_index], stateSpice, lt)
        for i in range(3):
            positionArray[i] = pyswice.doubleArray_getitem(stateSpice, i)
        position_secondaries_primary[:, ii] = positionArray

        # calculate the "third body" force
        f3rdBodies += -muSecondaries[ii] * \
           ( ( state[0:3] - position_secondaries_primary[:, ii] ) / \
             np.linalg.norm(state[0:3] - position_secondaries_primary[:, ii] ) ** 3 + \
             position_secondaries_primary[:, ii] / np.linalg.norm( position_secondaries_primary[:, ii] )**3 )

    # SRP force
    fSRP = cR * kSRP * state[0:3] / np.linalg.norm(state[0:3]) ** 3

    # total force (acceleration) vector
    f = fPrimary + f3rdBodies + fSRP + state[6:9]

    # args for the A matrix function
    args = (state[0:stateDimension], nSecondaries, muPrimary, muSecondaries, kSRP, cR,
            position_secondaries_primary)

    # A matrix calculation
    A = matrixA(args)

    # calculate the derivative of the STM
    dPhi = np.dot(A, phi)
    dPhi = np.reshape(dPhi,stateDimension * stateDimension)

    # acceleration vector to be returned to the integrator
    dState      = [0] * stateDimension 
    dState[0:6] = [state[3], state[4], state[5], f[0], f[1], f[2]]
    dState     += list(dPhi)

    return dState


################################################################################
#                       M A I N     F U N C T I O N:
################################################################################

def main():
    return


if __name__ == "__main__":
    main()
