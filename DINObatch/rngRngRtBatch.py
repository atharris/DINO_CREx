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
# import spiceypy as SP

import numpy as np

import math
################################################################################
#                  E X P O R T     F U N C T I O N S:
################################################################################

def fncH(input):
    # pull out the inputs for the H matrix
    state = input[0]
    SPICE_data = input[1]
    extras = input[-1]
    n_beacons  = len(extras['obs_beacons'])

    # count the number of QoIs
    n_state = state.shape[1]

    # initiate the H matrix
    H = np.zeros((2*n_beacons, n_state))

    # loop through beacons
    for ii in xrange(n_beacons):
        beacon_state = SPICE_data[ii,:]
        # calculate the difference between the positions and velocities
        r_diff = state[ii,0:3] - beacon_state[0:3]
        v_diff = state[ii,3:6] - beacon_state[3:6]

        # calculate the range
        rng = np.linalg.norm(r_diff)

        pv_summed = np.dot(r_diff, v_diff)

        H[2*ii, :] = np.array([r_diff[0]/rng, r_diff[1]/rng, r_diff[2]/rng,
                               0, 0, 0])
        H[2*ii + 1, :] = np.array([v_diff[0] / rng - r_diff[0] * pv_summed / (rng ** 3.),
                                   v_diff[1] / rng - r_diff[1] * pv_summed / (rng ** 3.),
                                   v_diff[2] / rng - r_diff[2] * pv_summed / (rng ** 3.),
                                   r_diff[0]/rng, r_diff[1]/rng, r_diff[2]/rng])

    return H

def fncG(input):
    # pull out the inputs for the generation of estimated pbservation data
    state = input[0]
    SPICE_data = input[1]
    extras = input[-1]
    n_beacons  = len(extras['obs_beacons'])

    # create array to be twice the size of the number of beacons. this is because there are
    # two data types: range and range rate
    G = np.zeros((n_beacons,2))

    for ii in xrange(n_beacons):
        beacon_state = SPICE_data[ii,:]
       
        # calculate the difference between the positions and velocites
        r_diff = state[ii, 0:3] - beacon_state[0:3]
        v_diff = state[ii, 3:6] - beacon_state[3:6]

        # calculate the range
        rng = np.linalg.norm(r_diff)
        G[ii,0] = rng
        # calculate the range rate
        rng_rate = np.divide(np.sum(np.multiply(r_diff, v_diff)), rng)
        G[ii,1] = rng_rate

    return G


###############################################################################
#                       M A I N     F U N C T I O N:
################################################################################

def main():
    return


if __name__ == "__main__":
    main()
