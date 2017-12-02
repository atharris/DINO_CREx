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

import math
################################################################################
#                  E X P O R T     F U N C T I O N S:
################################################################################

def fncH(input):
    # pull out the inputs for the H matrix
    state = input[0].flatten()
    SPICE_data = input[1]
    extras = input[-1]
    n_beacons  = len(extras['obs_beacons'])

    # count the number of QoIs
    n_state = state.shape[0]

    # initiate the H matrix
    H = np.zeros((2*n_beacons, n_state))

    # loop through beacons
    for ii in xrange(n_beacons):
        # pull out relevant key, add 1 to account for beacons[0] = SUN
        key = extras['obs_beacons'][ii]
        beacon_state = SPICE_data[key].flatten()
        # calculate the difference between the positions of the sun-to-sc and the sun-to-object
        # as well as the difference in velocities. These are used for multiple calculations
        r_diff = state[0:3] - beacon_state[0:3]
        v_diff = state[3:6] - beacon_state[3:6]

        # calculate the range and then tile it such that it is repeated 3 times
        rng = math.sqrt(np.sum(np.square(r_diff)))

        pv_summed = np.dot(r_diff, v_diff)

        H[2*ii, :] = np.array([r_diff[0]/rng, r_diff[1]/rng, r_diff[2]/rng,
                               0, 0, 0])
        H[2*ii + 1, :] = np.array([v_diff[0] / rng - r_diff[0] * pv_summed / (rng ** 3.),
                                   v_diff[1] / rng - r_diff[1] * pv_summed / (rng ** 3.),
                                   v_diff[2] / rng - r_diff[2] * pv_summed / (rng ** 3.),
                                   r_diff[0]/rng, r_diff[1]/rng, r_diff[2]/rng])

    return H

def fncG(input):
    # pull out the inputs for the generation of estimated observation data
    state = input[0]
    SPICE_data = input[1]
    extras = input[-1]
    n_beacons  = len(extras['obs_beacons'])

    # create array to be twice the size of the number of beacons. this is because there are
    # two data types: range and range rate
    G = np.zeros((state.shape[0], 2 * n_beacons))

    for ii in xrange(n_beacons):
        # pull out relevant key, add 1 to account for beacons[0] = SUN
        key = extras['obs_beacons'][ii]

        # calculate the difference between the positions of the sun-to-sc and the sun-to-object
        # as well as the difference in velocities. These are used for multiple calculations
        r_diff = state[:, 0:3] - SPICE_data[key].T[:, 0:3]
        v_diff = state[:, 3:6] - SPICE_data[key].T[:, 3:6]

        # calculate the range
        rng = np.array(np.sqrt(np.sum(np.square(r_diff), axis=1)))
        G[:, 0 + 2 * ii] = rng
        # calculate the range rate
        rng_rate = np.divide(np.sum(np.multiply(r_diff, v_diff), axis=1), rng)
        G[:, 1 + 2 * ii] = rng_rate

    return G


###############################################################################
#                       M A I N     F U N C T I O N:
################################################################################

def main():
    return


if __name__ == "__main__":
    main()
