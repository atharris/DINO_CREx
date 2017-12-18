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
dinoCommonPath = splitPath[0] + dinoName + '/DINObatch/pixelAndLine/commonFunctions/'
bskSpicePath = splitPath[0] + bskName + '/External/EphemerisData/'
bskPath = splitPath[0] + bskName + '/'
sys.path.append(bskPath + 'modules')
sys.path.append(bskPath + 'PythonModules')
sys.path.append(dinoSpicePath)
sys.path.append(dinoCommonPath)

try:
    import pyswice
except ImportError:
    from Basilisk import pyswice
    bskSpicePath = splitPath[0] + bskName + '/supportData/EphemerisData/'

import numpy as np
from batchFilter    import run_batch as runBatchVanilla
#from batchFilterAcc import run_batch as runBatchAcc

from batchFilter    import runRef

from plotIntegratedFilter import plotFunction as PF

import pickle

import pdb


## \defgroup init_batch_function initBatch - in-line callable batch filter
##   @{
## The module for running a DINO-CREx integrated batch filter.
#
# Overview {#overview}
# ====
#
# Purpose
# -----
# This module is a runable script for the batch filter when called from the "in house" perspective. For other functions in this configuration, please look for "in house" via the search bar. As a result of this, the `init.py` script is intended to run entirely independent of the rest of DINO-CREx and may be used for simplified scenarios or testing, among other applications. 
#
# Contents
# -----
# The `init.py` script contains a secondary function and a main function:
#
# - `main`
# - `writingText`
#
# It is here that we note `init.py` has largely served as a means to test the functionality of the other batch filter modules. It is provided for such reasons and is not considered an integral part of DINO-CREx. In fact, the DINO-CREx software package should never run a call of `init.py`, as it should be running the function found in \ref init_batch_function 
#
# The Code
# =====
#
# `main`
# -----
# Because `init.py` is a script meant to be run from the terminal, there are no inputs. The operator may choose to tweak various mission settings or spacecraft parameters, which will be covered shortly. One other consequence of this architecture is that there are also no outputs. The script does write to various .pkl files and create plots. However, these are not outputs in the Python function sense. 
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
#                  S E C O N D A R Y     F U N C T I O N S:
################################################################################

# -------------------------------------------------------------------------------

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

def initBatchFnc( stateValues, timeSpan, filterObservations, angles, extras ):

    ##################################################################################

    print '------------------'
    print 'Filter Image Span : ' ,(timeSpan[-1] - timeSpan[0])/(60*60*24), 'days'
    print '------------------'

    # pull out relevant state values from the associated input dictionary
    IC0         = stateValues['IC']
    phi0        = stateValues['phi0']
    covBar      = stateValues['covBar']
    stateDevBar = stateValues['stateDevBar']
    initialTime = stateValues['initial time']
    observationUncertainty = filterObservations['observation uncertainty']

    print 'IC', IC0

    # propagate the a priori state to the first observation time
    # THIS MATTER SHOULD BE INVESTIGATED FURTHER. 

    # prep pyswice for the extraction of initial data
    pyswice.furnsh_c(bskSpicePath  + extras['basic_bsp'])
    pyswice.furnsh_c(dinoSpicePath + extras['tls'])
    pyswice.furnsh_c(dinoSpicePath + extras['mission_bsp'])

    # time to first observation
    timeToObs       = np.append( np.array(initialTime), timeSpan[0] )
    propagatorInput = ( IC0[0:6], phi0[0:6,0:6], timeToObs, extras ) 

    # execute propagation
    IC_posVel  = runRef( propagatorInput )

    IC = np.append( IC_posVel[-1,:6], IC0[6:] )

    # initiate a filter output dictionary
    filterOutputs = {}

    ##################################################################################
    #
    # Run the Filter
    #
    ##################################################################################

   # run the filter and output the referenceStates (including STMs), est states and extra data
    for itr in xrange(extras['iterations']):

        if itr > 0:
            IC     = estimatedState[0, :]
            stateDevBar -= extraData['stateDevHatArray'][0, :]

        if itr==0:
            extras['oldPost'] = np.zeros([len(timeSpan), 2])

        # the arguments for the filter: the IC, the first STM, the time span, the observables
        # data dictionary, a priori uncertainty, and the measurables' uncertainty,
        # as well as any extras
        filterInputs = (IC, phi0, timeSpan, filterObservations,\
                         covBar, observationUncertainty, stateDevBar, angles, extras)
        # run filter function
        if extras['acc_est'] == 'OFF':
          referenceState, estimatedState, extraData = runBatchVanilla(filterInputs)
        elif extras['acc_est'] == 'ON':
          referenceState, estimatedState, extraData = runBatchAcc(filterInputs)
        else :
          print 'Acceleration estimation was not specified. Running Vanilla...'
          referenceState, estimatedState, extraData = runBatchVanilla(filterInputs)

        extras['oldPost'] = extraData['postfit residuals']

        # Check for anomaly:

        [anomaly_bool , anomaly_num] = extraData['anomaly_detected']
        if anomaly_bool == True:
            print '**********************************************************'
            print 'Anomaly Detected - Estimates are not to be trusted'
            print '**********************************************************'
            print anomaly_num, 'Residuals out of bounds'
        #    return

        # save all outputs into the dictionary with a name associated with the iteration
        filterOutputs[str(itr)] = {}
        filterOutputs[str(itr)]['referenceState'] = referenceState
        filterOutputs[str(itr)]['estimatedState'] = estimatedState
        filterOutputs[str(itr)]['extraData']      = extraData

        ##################################################################################
        #
        # \ BLOCK A page 196
        #
        ##################################################################################

        # Iteration Directory
        dirIt = '../Batch_Iteration' + str(itr+1)

        # Make directory for the iterations
        if not os.path.exists(dirIt):
            os.makedirs(dirIt)

        if extras['nav plots'] == 'ON':
          plotData = extraData

          plotData['postfit delta']   = extraData['postfit changes']
          plotData['states']          = estimatedState
          plotData['beacon_list']     = filterObservations['beaconIDs']
          plotData['timeSpan']        = timeSpan
          plotData['dirIt']           = dirIt
          plotData['obs_uncertainty'] = observationUncertainty
          plotData['referenceState']  = referenceState
          plotData['extras']          = extras
          plotData['acc_est']         = extras['acc_est']
          PF( plotData )

        #  Write the output to the pickle file
        fileTag = 'nominal'
        file = dirIt+'/'+fileTag+'_data.pkl'
        pklFile = open( file, 'wb')
        pickle.dump( plotData, pklFile, -1 )
        pklFile.flush()

        pklFile.close()

    return filterOutputs

if __name__ == "__main__":
    main()