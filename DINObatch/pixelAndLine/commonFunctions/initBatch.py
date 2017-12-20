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
from batchFilterAcc import run_batch as runBatchAcc

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
# This module is similar to those of \ref init_vanilla or \ref init_acc , except it is intended to be run as function call from a larger script/module such as a DINO sim.  
#
# Contents
# -----
# The `initBatch.py` contains a single exportable function:
#
# - `initBatchFnc`
#
# This function serves to unpack inputs, determine which batch filter to run (currently \ref batch_vanilla or \ref batch_acc ) and loop through iterations before packaging the outputs to be passed to the module from where the function was called.
#
# The Code
# =====
#
# `initBatchFnc`
# -----
# Before we discuss the code contained within this module, we list the inputs needed:
#
# Name      | Description                                           | Size/Type
# -----     | -------------------                                   | -----
# stateValues | dictionary of arrays related to the quantities of interest | dictionary  
# timeSpan  | times at which observations occur | (N,) numpy array
# filterObersvations | dictionary of arrays related to measured data  | dictionary
# angles    | angles at which the observations are taken at each time in `timeSpan`  | numpy array (N,3)
# extras    | dictionary of various parameters  | dictionary
#
# The same is repeated for outputs:
#
# Name      | Description         | Size/Type
# -----     | ------------------- | -----
# filterOutputs | various results of the filter. One key for each iteration | dictionary
#
# Although `initBatch.py` is largely similar in function to the "init" class scripts of \ref init_vanilla and \ref init_acc , there are some differences. Notably in the first lines, some code is needed to bridge any gaps between the provided state and that of the first observation, i.e.,
# ~~~~~~~~~~~~~~~~{.py}
#    # time to first observation
#    timeToObs       = np.append( np.array(initialTime), timeSpan[0] )
#    propagatorInput = ( IC0[0:6], phi0[0:6,0:6], timeToObs, extras ) 
#
#    # execute propagation
#    IC_posVel  = runRef( propagatorInput )
#
#    IC = np.append( IC_posVel[-1,:6], IC0[6:] )
# ~~~~~~~~~~~~~~~~
# This is needed because the batch filter algorithms provided for DINO-CREx process an observation for each time given. In the case of the inputs, the first time must have an observation. Therefore, if the provided state initial conditition occurs before the first observation, a bit of propagation is needed to align the first state with the first observation. 
# 
# After this initial propagation, `initBatchFnc()` behaves almost identically with other batch initialization scripts. One other notable change, the outputs are packaged into a single dictionary `filterOutputs`:
# ~~~~~~~~~~~~~~~~{.py}
#        # save all outputs into the dictionary with a name associated with the iteration
#        filterOutputs[str(itr)] = {}
#        filterOutputs[str(itr)]['referenceState'] = referenceState
#        filterOutputs[str(itr)]['estimatedState'] = estimatedState
#        filterOutputs[str(itr)]['extraData']      = extraData
# ~~~~~~~~~~~~~~~~
# This dictionary is then passed to where `initBatch.py` was called from. Note, the dictionaty of `filterOutputs` is organized by iteration, e.g., ``filterOutputs['0']`` if the first iteration.
#
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
