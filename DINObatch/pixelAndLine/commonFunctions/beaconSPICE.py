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

from rngRngRt import fncH
from rngRngRt import fncG

import numpy as np
from numpy.random import normal as rndNrm


################################################################################
#                  E X P O R T E D     F U N C T I O N S:
################################################################################

# -------------------------------------------------------------------------------
def norm(input):
    norm = np.sqrt(sum(np.square(input)))
    return norm


def getObs(input):
    # pull out the inputs for the generation of observation data
    SPICE_data = input[0]
    observation_uncertainty = input[1]
    extras = input[-1]
    n_beacons = extras['n_beacons']

    # for simplification, give the location of the SC with respect to the sun it's own variable
    r_spacecraft = SPICE_data['spacecraft']

    # number of samples/observations
    n_samples = r_spacecraft.shape[1]

    # create noise for the estimated observations
    observation_noise = rndNrm(0, 1, (2, n_samples, n_beacons))

    # create array to be twice the size of the number of beacons. this is because there are
    # two data types: range and range rate
    obs = np.zeros((n_samples, 2 * n_beacons))

    for ii in xrange(n_beacons):
        # pull out relevant key
        key = extras['beacons'][ii]

        # calculate the difference between the positions of the sun-to-sc and the sun-to-object
        # as well as the difference in velocities. These are used for multiple calculations
        r_diff = r_spacecraft[0:3, :] - SPICE_data[key][0:3, :]
        v_diff = r_spacecraft[3:6, :] - SPICE_data[key][3:6, :]

        # calculate the range
        rng = np.array(np.sqrt(np.sum(np.square(r_diff), axis=0)))
        obs[:, 0 + 2 * ii] = rng

        # calculate the range rate
        rng_rate = np.divide(np.sum(np.multiply(r_diff, v_diff), axis=0), rng)
        obs[:, 1 + 2 * ii] = rng_rate

        # store these in the observation array
        obs[:, 0 + ii * 2:2 + ii * 2] += \
            np.dot(np.linalg.cholesky(observation_uncertainty), observation_noise[:, :, ii]).T

    return obs


################################################################################
#                    E X P O R T E D     C L A S S E S:
################################################################################

# -------------------------------------------------------------------------------
#
# -------------------------------------------------------------------------------

################################################################################
#             U N I T     T E S T     C A S E     F U N C T I O N:
################################################################################

# -------------------------------------------------------------------------------
#
# -------------------------------------------------------------------------------

################################################################################
#                       M A I N     F U N C T I O N:
################################################################################

def main():
    return


if __name__ == "__main__":
    main()
