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

import os
import pickle
import sys
import time
import optparse
import socket
import pdb
import scipy.integrate as integ
import scipy.io as io
import spiceypy as SP
import numpy as np


################################################################################
#                  E X P O R T E D     F U N C T I O N S
################################################################################

# -------------------------------------------------------------------------------

def matrixA(input):
    n_state        = len(input[0])
    r_spacecraft   = np.expand_dims(input[0][0:3],axis=1)
    n_secondaries  = input[1]
    mu_primary     = input[2]
    mu_secondaries = input[3]
    kSRP           = input[4]
    cR             = input[5]
    r_sun          = np.expand_dims(input[6],axis=1)
    r_secondaries_primary = input[7]

    # set the size for the A matrix and premptively populate with zeros
    A = np.zeros((n_state, n_state))

    # the r_spacecraft derivate associated with the primary gravitational force
    dFdR_p = -mu_primary * (
    np.identity(3) / np.linalg.norm(r_spacecraft) ** 3 -
    3 * np.dot(r_spacecraft, r_spacecraft.T) / np.linalg.norm(r_spacecraft) ** 5 
    )

    # the r_spacecraft derivatives associated with gravitational force from secondary bodies
    dFdR_s = np.zeros( ( 3, 3 ) )

    # loop through the secondary bodies
    for ii in xrange(n_secondaries) :
        r_secondary = np.expand_dims( r_secondaries_primary[:, ii], axis = 1 )
        dFdR_s += -mu_secondaries[ii] * (
        np.identity(3) / np.linalg.norm(r_spacecraft - r_secondary) ** 3 -
        3 * np.dot(r_spacecraft - r_secondary,(r_spacecraft - r_secondary).T)/
        np.linalg.norm(r_spacecraft - r_secondary) ** 5
        )

    # the r_spacecraft derivative associated with the SRP force
    dFdR_SRP = cR * kSRP * (np.identity(3) / np.linalg.norm(r_spacecraft - r_sun) ** 3 -
                            3 * np.dot(r_spacecraft - r_sun, (r_spacecraft - r_sun).T) /
                            np.linalg.norm(r_spacecraft - r_sun) ** 5)

    # total r_spacecraft derivatives of forces
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


def EOM(state, et, primary_index, secondary_indices, n_secondaries, mu_primary, mu_secondaries,
        kSRP, cR, abcorr, ref_frame, bodies, n_state):

    # pull out the STM
    phi = np.array(state[n_state:], copy=True)
    phi = np.reshape(phi, (n_state, n_state) )

    # gravitational force from primary body
    f_primary = -mu_primary * state[0:3] / np.linalg.norm(state[0:3]) ** 3

    # gravitational force from secondary bodies
    f_3rd_bodies = 0

    # set the size of the r_spacecraftition_secondaries between 
    # secondary bodies and primary body
    r_secondaries_primary = np.zeros((3, n_secondaries))

    # loop through the secondary bodies
    for ii in xrange(n_secondaries):
        # determine distance from secondary to primary body
        r_secondaries_primary[:, ii] = \
             SP.spkezr(bodies[secondary_indices[ii]], et, ref_frame, 
                       abcorr, bodies[primary_index])[0][0:3]
        # calculate the "third body" force
        f_3rd_bodies += -mu_secondaries[ii] * \
           ( ( state[0:3] - r_secondaries_primary[:, ii] ) / \
             np.linalg.norm(state[0:3] - r_secondaries_primary[:, ii] ) ** 3 + \
             r_secondaries_primary[:, ii] / np.linalg.norm( r_secondaries_primary[:, ii] )**3 )

    # r_spacecraftition of sun with respect to primary body
    r_sun = SP.spkezr(bodies[primary_index], et, ref_frame, abcorr, bodies[0])[0][0:3]

    # SRP force
    f_SRP = cR * kSRP * (state[0:3] - r_sun) / np.linalg.norm(state[0:3] - r_sun) ** 3

    # total force (acceleration) vector
    f = f_primary + f_3rd_bodies + f_SRP

    # args for the A matrix function
    args = (state[0:n_state], n_secondaries, mu_primary, mu_secondaries, kSRP, cR, r_sun,
            r_secondaries_primary)

    # A matrix calculation
    A = matrixA(args)

    # calculate the derivative of the STM
    dPhi = np.dot(A, phi)
    dPhi = np.reshape(dPhi,n_state * n_state)

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
