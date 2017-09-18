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

from numpy.random import normal as rndNrm


################################################################################
#                  E X P O R T E D     F U N C T I O N S:
################################################################################

# -------------------------------------------------------------------------------

def getObs(input):
    # pull out the inputs for the generation of observation data
    SPICE_data = input[0]
    observation_uncertainty = input[1]
    extras = input[-1]
    n_beacons = extras['n_beacons']

    # number of observations in a mini observation set per beacon
    n_obs = extras['n_obs']

    # for simplification, give the location of the SC with respect to the sun it's own variable
    r_spacecraft = SPICE_data['spacecraft']

    # number of samples/observations
    n_samples = r_spacecraft.shape[1]

    # create noise for the estimated observations
    np.random.seed(42)
    observation_noise = rndNrm(0, 1, (2, n_samples, n_beacons))



    # the size of an observation "bin". A "bin" is a super set of 
    # observations for each beacon. Each beacon is observed
    # until it has been observed n_obs times. 
    # After this, the next beacon is observed and so on. This
    # pattern repeats until the last beacon in the list of 
    # possible objects has been reached. After this, the pattern repeats once again
    bin_size = n_obs * n_beacons

    # The largest number of whole bins that it will take to get as close as possible to the
    # amount of provided data samples
    n_bins = n_samples / bin_size

    # initialize the observation dictionary
    obs = {}

    # create a list for the beacon names associated with each measurement
    obs['beacons'] = list()
    obs['data']    = np.zeros((n_samples, 2))
    obs['truth']   = np.zeros((n_samples, 2))
    obs['SPICE']   = np.zeros((n_samples, 6))
    # loop through all the bins. make sure to go max bin + 1 so that the rest of the 
    # arrays get filled even after the last whole bin
    for bb in xrange( int(n_bins + 1) ) :    

      # while in a bin, loop through beacon keys
      for ii in xrange(n_beacons):
         # if the indices for the beacon mini bin are not outside the last whole bin, add the data
         if bb * bin_size + n_obs + n_obs * ii < n_samples :
            # create beacon mini bin indices for readability
            start_idx = bb * bin_size + n_obs * ii
            end_idx   = bb * bin_size + n_obs + n_obs * ii

            # pull out relevant key
            key = extras['beacons'][ii]
 
            # store the key associated with the data
            obs['beacons'] += [key] * n_obs
  
            # calculate the difference between the positions of 
            # the sun-to-sc and the sun-to-object
            # as well as the difference in velocities. 
            # These are used for multiple calculations
            r_beacon = np.copy(SPICE_data[key][0:3, start_idx : end_idx ])
            v_beacon = np.copy(SPICE_data[key][3:6, start_idx : end_idx ])

            # store the beacon data in an array to be used for the calculation
            # of H matrices and estimated observations.
            obs['SPICE'][start_idx : end_idx, 0:3 ] = r_beacon.T
            obs['SPICE'][start_idx : end_idx, 3:6 ] = v_beacon.T

            r_diff = r_spacecraft[0:3, start_idx : end_idx ] - \
                     r_beacon #+ np.dot(np.linalg.cholesky(observation_uncertainty[0:3,0:3]), \
                     #observation_noise[0:3, start_idx : end_idx, ii])

            v_diff = r_spacecraft[3:6, start_idx : end_idx ] - \
                     v_beacon #+ np.dot(np.linalg.cholesky(observation_uncertainty[3:6,3:6]), \
                     #observation_noise[3:6, start_idx : end_idx, ii])

            # calculate the range 

            rng = np.array(np.sqrt(np.sum(np.square(r_diff), axis=0)))
            obs['data'][start_idx : end_idx, 0] = rng

            # calculate the range rate
            rng_rate = np.divide(np.sum(np.multiply(r_diff, v_diff), axis=0), rng)
            obs['data'][start_idx : end_idx, 1] = rng_rate

            # store these in the observation array
            obs['truth'][start_idx : end_idx, 0:2] = np.copy( obs['data'][start_idx : end_idx, 0:2] )
            obs['data'][start_idx : end_idx, 0:2] += \
                np.dot(np.linalg.cholesky(observation_uncertainty), \
                observation_noise[:, start_idx : end_idx, ii]).T

         # if the indices are outside a whole bin, do a special data fill
         else :
            # begin where the begin should
            start_idx = bb * bin_size + n_obs * ii
            # but end at the last data sample
            end_idx   = n_samples

            # pull out relevant key
            key = extras['beacons'][ii]

            # copy the relevant keys in
            obs['beacons'] += [key] * (end_idx - start_idx)

            # calculate the difference between the positions of 
            # the sun-to-sc and the sun-to-object
            # as well as the difference in velocities. 
            # These are used for multiple calculations
            r_beacon = np.copy(SPICE_data[key][0:3, start_idx : end_idx ])
            v_beacon = np.copy(SPICE_data[key][3:6, start_idx : end_idx ])

            # store the beacon data in an array to be used for the calculation
            # of H matrices and estimated observations.

            obs['SPICE'][start_idx : end_idx, 0:3 ] = r_beacon.T
            obs['SPICE'][start_idx : end_idx, 3:6 ] = v_beacon.T

            r_diff = r_spacecraft[0:3, start_idx : end_idx ] - \
                     r_beacon #+ np.dot(np.linalg.cholesky(observation_uncertainty[0:3,0:3]), \
                     #observation_noise[0:3, start_idx : end_idx, ii])

            v_diff = r_spacecraft[3:6, start_idx : end_idx ] - \
                     v_beacon #+ np.dot(np.linalg.cholesky(observation_uncertainty[3:6,3:6]), \
                     #observation_noise[3:6, start_idx : end_idx, ii])

            # calculate the range  
            rng = np.array(np.sqrt(np.sum(np.square(r_diff), axis=0)))
            obs['data'][start_idx : end_idx, 0] = rng

            # calculate the range rate
            rng_rate = np.divide(np.sum(np.multiply(r_diff, v_diff), axis=0), rng)
            obs['data'][start_idx : end_idx, 1] = rng_rate

            # store these in the observation array
            obs['truth'][start_idx : end_idx, 0:2] = np.copy( obs['data'][start_idx : end_idx, 0:2] )
            obs['data'][start_idx : end_idx, 0:2] += \
                np.dot(np.linalg.cholesky(observation_uncertainty), \
                observation_noise[:, start_idx : end_idx, ii]).T

            # after this last data gets pulled, break the bin loop
            break

    return obs

def putObs( input ) :

   data   = input[0]
   extras = input[-1]

   n_samples = extra_data['Y']['data'].shape[0]

   n_obs     = extra_data['Y']['data']

   obs_data = {}

   for bb in xrange( extras['n_beacons'] ) :
      key = extras['beacons'][bb]
      obs_data[key] = {}
      

   for ii in xrange( n_samples ) :
      print wow
      
   return obs_data

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
