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

################################################################################
#                  E X P O R T E D     F U N C T I O N S:
################################################################################

# -------------------------------------------------------------------------------

def getObs(input):
    # pull out the inputs for the generation of observation data
    spiceData              = input[0]
    observationUncertainty = input[1]
    angles = input[2]
    extras                 = input[-1]
    nUniqueBeacons         = extras['n_unique_beacons']
    # number of observations in a mini observation set per beacon
    repeatObservations     = extras['repeat_obs']

    # for simplification, give the location of the SC with respect to the sun it's own variable
    spacecraftPosition = spiceData['spacecraft']

    # number of samples/observations
    nObservations = spacecraftPosition.shape[1]

    # create noise for the estimated observations
    # np.random.seed(42)
    if 'seed' in extras:
      np.random.seed(extras['seed'])

    observationNoise = np.zeros([nObservations, 2])
    observationNoise[:,0] = rndNrm(0., observationUncertainty[0,0] , nObservations)
    observationNoise[:,1] = rndNrm(0., observationUncertainty[1,1] , nObservations)

    # the size of an observation "bin". A "bin" is a super set of 
    # observations for each beacon. Each beacon is observed
    # until it has been observed repeatObservations times. 
    # After this, the next beacon is observed and so on. This
    # pattern repeats until the last beacon in the list of 
    # possible objects has been reached. After this, the pattern repeats once again
    binSize = repeatObservations * nUniqueBeacons
    # The largest number of whole bins that it will take to get as close as possible to the
    # amount of provided data samples
    nBins   = nObservations / binSize

    # initialize the observation dictionary
    obs = {}

    # create a list for the beacon names associated with each measurement
    obs['beacons'] = list()
    obs['measurements']    = np.zeros((nObservations, 2))
    obs['truth']   = np.zeros((nObservations, 2))
    obs['SPICE']   = np.zeros((nObservations, 6))

    referenceState = np.zeros([nObservations, 6])
    # loop through all the bins. make sure to go max bin + 1 so that the rest of the 
    # arrays get filled even after the last whole bin
    for bb in xrange( int(nBins ) ) :
      # while in a bin, loop through beacon keys
      for ii in xrange(nUniqueBeacons):
         # if the indices for the beacon mini bin are not outside the last whole bin, add the data
         indice =  bb * binSize + repeatObservations * ii
         if indice < nObservations :
            # create beacon mini bin indices for readability
            startIndex = bb * binSize + repeatObservations * ii
            endIndex   = bb * binSize + repeatObservations + repeatObservations * ii

            # pull out relevant key
            key = extras['unique_beacon_IDs'][ii]
 
            # store the key associated with the data
            obs['beacons'] += [key] * repeatObservations
  
            # calculate the difference between the positions of 
            # the sun-to-sc and the sun-to-object
            # as well as the difference in velocities. 
            # These are used for multiple calculations
            beaconPosition = np.copy(spiceData[key][0:3, startIndex : endIndex ])
            beaconVelocity = np.copy(spiceData[key][3:6, startIndex : endIndex ])

            # store the beacon data in an array to be used for the calculation
            # of H matrices and estimated observations.
            obs['SPICE'][startIndex : endIndex, 0:3 ] = beaconPosition.T
            obs['SPICE'][startIndex : endIndex, 3:6 ] = beaconVelocity.T

            rSC = np.copy(spacecraftPosition[:, startIndex : endIndex])
            referenceState[indice ,:] = rSC.T

         # if the indices are outside a whole bin, do a special data fill
         else :
            # begin where the begin should
            startIndex = bb * binSize + repeatObservations * ii
            # but end at the last data sample
            endIndex   = nObservations

            # pull out relevant key
            key = extras['unique_beacon_IDs'][ii]

            # copy the relevant keys in
            obs['beacons'] += [key] * (endIndex - startIndex)

            # calculate the difference between the positions of 
            # the sun-to-sc and the sun-to-object
            # as well as the difference in velocities. 
            # These are used for multiple calculations
            beaconPosition = np.copy(spiceData[key][0:3, startIndex : endIndex ])
            beaconVelocity = np.copy(spiceData[key][3:6, startIndex : endIndex ])

            # store the beacon data in an array to be used for the calculation
            # of H matrices and estimated observations.

            obs['SPICE'][startIndex : endIndex, 0:3 ] = beaconPosition.T
            obs['SPICE'][startIndex : endIndex, 3:6 ] = beaconVelocity.T

            rSC = np.copy(spacecraftPosition[:, startIndex : endIndex])
            # print bb * binSize  + repeatObservations * ii
            referenceState[indice,:] = rSC.T

            break
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
