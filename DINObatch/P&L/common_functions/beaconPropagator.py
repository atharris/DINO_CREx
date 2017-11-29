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

try:
    import pyswice
except ImportError:
    from Basilisk import pyswice
    bskSpicePath = splitPath[0] + bskName + '/supportData/EphemerisData/'

import numpy as np

import pdb

################################################################################
#                  E X P O R T E D     F U N C T I O N S:
################################################################################

# -------------------------------------------------------------------------------

def beaconStates(input):
    # pull out the inputs for the generation of observation data
    beacon_list       = input[0]
    observation_times = input[1]
    extras            = input[-1]

    # compute and save total number of observations
    n_observations    = len( beacon_list )

    # Load Ephemeris Files from extras dictionary
    pyswice.furnsh_c(bskSpicePath  + extras['basic_bsp'])
    pyswice.furnsh_c(dinoSpicePath + extras['mission_bsp'])
    pyswice.furnsh_c(dinoSpicePath + extras['tls'])

    # instantiate beacon_states
    beacon_states = []

    # for each observation time, compute the state of the dictated beacon
    for ii in xrange( n_observations ):
        stateArray = np.zeros(6)
        state = pyswice.new_doubleArray(6)
        lt = pyswice.new_doubleArray(1)
        pyswice.spkezr_c( beacon_list[ii], observation_times[ii], extras['ref_frame'],\
                          'None', 'SUN', state, lt)
        for i in range(6):
            stateArray[i] = pyswice.doubleArray_getitem(state, i)
        beacon_states.append(stateArray)
    
    # convert final data set into numpy array
    beacon_states = np.array(beacon_states)

    return beacon_states

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
