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
    sunPosition        = np.expand_dims(input[6],axis=1)
    position_secondaries_primary = input[7]

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
       np.linalg.norm(spacecraftPosition - sunPosition) ** 3 -
       3 * np.dot(spacecraftPosition - sunPosition, (spacecraftPosition - sunPosition).T) /\
                            np.linalg.norm(spacecraftPosition - sunPosition) ** 5)

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
             position_secondaries_primary[:, ii] / \
             np.linalg.norm( position_secondaries_primary[:, ii] )**3 )

    # spacecraftPositionition of sun with respect to primary body
    sunPositionArray = np.zeros(3)
    stateSpice = pyswice.new_doubleArray(6)
    lt = pyswice.new_doubleArray(1)
    pyswice.spkezr_c(bodies[secondary_indices[ii]], et, ref_frame,
                     abcorr, bodies[primary_index], stateSpice, lt)
    for i in range(3):
        sunPositionArray[i] = pyswice.doubleArray_getitem(stateSpice, i)
    sunPosition = sunPositionArray

    # SRP force
    f_SRP = cR * kSRP * (state[0:3] - sunPosition) / \
            np.linalg.norm(state[0:3] - sunPosition) ** 3

    # total force (acceleration) vector
    f = fPrimary + f3rdBodies + f_SRP

    # args for the A matrix function
    args = (state[0:stateDimension], nSecondaries, muPrimary, muSecondaries, kSRP, cR,\
            sunPosition, position_secondaries_primary)

    # A matrix calculation
    A = matrixA(args)

    # calculate the derivative of the STM
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
