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
from numpy.random import normal as rndNrm
from pixelLineBatch import fncG
import pdb

## \defgroup beacon_bin beaconBinSPICE - measurement creation in house
##   @{
## The module for creating "in house" simulation measurement data for the batch filter.
#
# Overview {#overview}
# ====
#
# Purpose
# -----
# This script contains the exportable function `getObs` used for the purpose of creating measurement data (observations) for the batch filter when called from the "in house" perspective. For other functions in this configuration, please look for "in house" via the search bar
#
# Contents
# -----
# The following exportable function is contained in this module:
#
# - `getObs`
#
# The `getObs` function is used for the creation of observations that take the form of measured data at desired times. This data is analyzed in a batch filter called outside of a BSK simulation. This standalone functionality is referred to as the "in house" configuration of the batch. It is used when testing or other needs require the use of the Navigation Software portion of DINO-CREx. 
#
# The Code
# =====
#
# `getObs`
# -----
# Using inputs from \ref data_creation, `getObs` provides an ability to create observed values needed for the "in house" batch filter run of DINO-CREx. A table of required input is provided below::
#
# Name      | Description                                           | Size/Type
# -----     | -------------------                                   | -----
# ephemerisData | Dictonary containing the state of each beacon at every observation time | dictionary
# observationUncertainty | Uncertainty to be appleid to the measurements | list of str
# angles | Array of spacecraft attitude at each observation time | (N,3) Numpy array
# extras | Dictionary of various parameters | dictionary
#
# The output of this function is a dictionary of measurements
#
# Name      | Description         | Size/Type
# -----     | ------------------- | -----
# obs | dictonairy of true and measured observations | dictionary
#
# The first significant lines of this function handle the noise to be added to the true measurements. 
# ~~~~~~~~~~~~~~~~{.py}
#    # create noise for the estimated observations
#    if 'seed' in extras:
#      np.random.seed(extras['seed'])
#
#    observationNoise = np.zeros([nObservations, 2])
#    observationNoise[:,0] = rndNrm(0., observationUncertainty[0,0] , nObservations)
#    observationNoise[:,1] = rndNrm(0., observationUncertainty[1,1] , nObservations)
# ~~~~~~~~~~~~~~~~
#
# It is important to note that the measurement uncertainty is applied via a Gaussian noise that has the standard deviation associated with the `observationUncertainty` variable. the variable of `extras['seed']` allows an operator to produce repeatable results, if desired. 
#  
# A concept utilized in this function is a "bin". These bins are composed of an observation for each beacon in succession. For example, if Mars, Earth and the Moon are beacons a bin would be ['Mars', 'Earth', 'Moon']. Each bin is then repeated until the total number of observations is satisfied. The bin becomes important if the number of total observations is not divisible by the length of these repeatable bins. This consideration is seen in the dictated loop statements:
# ~~~~~~~~~~~~~~~~{.py}
#    # loop through all the bins
#    for bb in xrange( nBins ) :
#      # while in a bin, loop through beacon keys
#      for ii in xrange(nUniqueBeacons):
#            index = bb * nUniqueBeacons + ii
#
#            # pull out relevant key
#            beaconKey = extras['unique_beacon_IDs'][ii]
# ~~~~~~~~~~~~~~~~
# followed by
# ~~~~~~~~~~~~~~~~{.py}
#    # fill out the partial bin
#    for ii in xrange( partialBinLength ):
#      index = nObsInFullBins + ii
#      # pull out relevant key
#      beaconKey = extras['unique_beacon_IDs'][ii]
# ~~~~~~~~~~~~~~~~
#
# These loops collect the appropriate beacon state and the spacecraft state for each measurement time. The module then creates observation data by using these states as inputs to `fncG()` found in \ref pixel_and_line, i.e.,
# ~~~~~~~~~~~~~~~~{.py}
#    G_ref_inputs = (referenceState, obs['SPICE'], angles, extras)
#    # calculate the estimated observables and organize into an array
#    obs['truth'] = np.copy(fncG(G_ref_inputs))
#    obs['measurements'] = np.copy(fncG(G_ref_inputs)) + observationNoise
## ~~~~~~~~~~~~~~~~
## @}

################################################################################
#                  E X P O R T E D     F U N C T I O N S:
################################################################################

# -------------------------------------------------------------------------------

def getObs(input):
    # pull out the inputs for the generation of observation data
    ephemerisData          = input[0]
    observationUncertainty = input[1]
    angles                 = input[2]
    extras                 = input[-1]
    nUniqueBeacons         = extras['n_unique_beacons']

    # for simplification, give the location of the SC with respect to the sun it's own variable
    spacecraftPosition = ephemerisData['spacecraft']

    # number of samples/observations
    nObservations = spacecraftPosition.shape[1]

    # create noise for the estimated observations
    if 'seed' in extras:
      np.random.seed(extras['seed'])

    observationNoise = np.zeros([nObservations, 2])
    observationNoise[:,0] = rndNrm(0., observationUncertainty[0,0] , nObservations)
    observationNoise[:,1] = rndNrm(0., observationUncertainty[1,1] , nObservations)

    # calculate the number of observation bins. use int() to round down
    nBins            = int(nObservations / nUniqueBeacons)

    # length of full bins
    nObsInFullBins = nBins * nUniqueBeacons

    # number of leftover observations
    partialBinLength = int(nObservations - nObsInFullBins)

    # initialize the observation dictionary
    obs = {}

    # create a list for the beacon names associated with each measurement
    obs['beacons']         = list()
    obs['measurements']    = np.zeros((nObservations, 2))
    obs['truth']           = np.zeros((nObservations, 2))
    obs['SPICE']           = np.zeros((nObservations, 6))

    referenceState = np.zeros([nObservations, 6])
    # loop through all the bins
    for bb in xrange( nBins ) :
      # while in a bin, loop through beacon keys
      for ii in xrange(nUniqueBeacons):
            index = bb * nUniqueBeacons + ii

            # pull out relevant key
            beaconKey = extras['unique_beacon_IDs'][ii]
 
            # store the key associated with the data
            obs['beacons'] += [beaconKey]
  
            # pull out beacon position and velocity
            beaconPosition = np.copy(ephemerisData[beaconKey][0:3, index ])
            beaconVelocity = np.copy(ephemerisData[beaconKey][3:6, index ])

            # store the beacon data in an array to be used for the calculation
            # of H matrices and estimated observations.
            obs['SPICE'][index, 0:3 ] = beaconPosition.T
            obs['SPICE'][index, 3:6 ] = beaconVelocity.T
            
            # pull out spacecraft state
            rSC = np.copy(spacecraftPosition[:, index])
            referenceState[index ,:] = rSC.T

    # fill out the partial bin
    for ii in xrange( partialBinLength ):
      index = nObsInFullBins + ii
      # pull out relevant key
      beaconKey = extras['unique_beacon_IDs'][ii]

      # store the key associated with the data
      obs['beacons'] += [beaconKey]

      # pull out beacon position and velocity
      beaconPosition = np.copy(ephemerisData[beaconKey][0:3, index ])
      beaconVelocity = np.copy(ephemerisData[beaconKey][3:6, index ])

      # store the beacon data in an array to be used for the calculation
      # of H matrices and estimated observations.
      obs['SPICE'][index, 0:3 ] = beaconPosition.T
      obs['SPICE'][index, 3:6 ] = beaconVelocity.T
      pdb.set_trace()
      # pull out spacecraft state
      rSC = np.copy(spacecraftPosition[:, index])
      referenceState[index ,:] = rSC.T



    extras['obs_beacons'] = list(obs['beacons'])
    G_ref_inputs = (referenceState, obs['SPICE'], angles, extras)
    # calculate the estimated observables and organize into an array
    obs['truth'] = np.copy(fncG(G_ref_inputs))
    obs['measurements'] = np.copy(fncG(G_ref_inputs)) + observationNoise
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
