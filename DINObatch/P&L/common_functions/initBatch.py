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
dinoCommonPath = splitPath[0] + dinoName + '/DINObatch/P&L/common_functions/'
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
from batchFilter import run_batch
import data_generation as dg
from plotFilter import plotFunction as PF
from beaconBinSPICE import getObs
import pickle

import pdb

################################################################################
#                  S E C O N D A R Y     F U N C T I O N S:
################################################################################

# -------------------------------------------------------------------------------

def writingText(itr, referenceState, estimatedState, trueEphemeris, extraData, initialPositionError , initialVelocityError):
    # calculate the difference between the perturbed reference and true trajectories: reference state errors
    err = referenceState[:, 0:6] - trueEphemeris['spacecraft'].T

    # compare the estimated and true trajectories: estimated state errors
    stateDevHat = estimatedState[:, 0:6] - trueEphemeris['spacecraft'].T

    resultString = ''

    resultString += '---------------------------------------------------' + '\n'
    resultString += 'Iteration number '+ str(itr) + '\n'
    resultString += '---------------------------------------------------'+ '\n'
    resultString += '\n'
    resultString += 'Estimated x_hat_0 = ' + str(extraData['stateDevHat'])+ '\n'
    resultString += 'Actual Error = ' + str(initialPositionError) + str(initialVelocityError) + '\n'
    resultString += '\n'

    resultString += 'Ref X Pos err = ' + str(err[-1, 0]) + '\n'
    resultString += 'Ref Y Pos err = ' + str(err[-1, 1]) + '\n'
    resultString += 'Ref Z Pos err = ' + str(err[-1, 2]) + '\n'
    resultString += 'Ref X Vel err = ' + str(err[-1, 3]) + '\n'
    resultString += 'Ref Y Vel err = ' + str(err[-1, 4]) + '\n'
    resultString += 'Ref Z Vel err = ' + str(err[-1, 5]) + '\n'
    resultString += '\n'
    resultString += 'Est X Pos err = ' + str(stateDevHat[-1, 0]) + '\n'
    resultString += 'Est Y Pos err = ' + str(stateDevHat[-1, 1]) + '\n'
    resultString += 'Est Z Pos err = ' + str(stateDevHat[-1, 2]) + '\n'
    resultString += 'Est X Vel err = ' + str(stateDevHat[-1, 3]) + '\n'
    resultString += 'Est Y Vel err = ' + str(stateDevHat[-1, 4]) + '\n'
    resultString += 'Est Z Vel err = ' + str(stateDevHat[-1, 5]) + '\n'
    resultString += '\n'

    print resultString

    text_file = open('Batch_Iteration' + str(itr) + "/Batch" + str(itr) + ".txt", "w")
    text_file.write(resultString)
    text_file.close()

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

def initBatchFnc():

    ##################################################################################

    print '------------------'
    print 'Filter Image Span : ' ,(timeSpan[-1] - timeSpan[0])/(60*60*24), 'days'
    print '------------------'

    # copy the initial conditions as the first sun to SC referenceStates from the SPICE file
    IC = np.copy(trueEphemeris['spacecraft'][:, 0])

    print 'IC', IC

    # initiate a filter output dictionary
    filterOutputs = {}

    ##################################################################################
    #
    # Get the noisy observations
    #
    ##################################################################################

    # create dictionary for observation data to be inputs in filter. This is a more limited
    # dictionary than dataObservations and serves as the most "real" input
    filterObservations = {}
    filterObservations['measurements'] = dataObservations['measurements']
    filterObservations['beaconIDs']    = dataObservations['beacons']

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
        referenceState, estimatedState, extraData = run_batch(filterInputs)
        extras['oldPost'] = extraData['postfit residuals']

        # Check for anomaly:

        [anomaly_bool , anomaly_num] = extraData['anomaly_detected']
        if anomaly_bool == True:
            print '**********************************************************'
            print 'Anomaly Detected - Estimates are not to be trusted'
            print '**********************************************************'
            print anomaly_num, 'Residuals out of bounds'
            return

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
        dirIt = 'Batch_Iteration' + str(itr+1)

        # Make directory for the iterations
        if not os.path.exists(dirIt):
            os.makedirs(dirIt)

        # File to write data
        writingText( itr+1, referenceState, estimatedState, trueEphemeris, extraData,\
                     initialPositionError , initialVelocityError)

        # calculate the difference between the perturbed reference and 
        # true trajectories: reference state errors
        stateError = referenceState[:, 0:6] - trueEphemeris['spacecraft'].T

        # compare the estimated and true trajectories: estimated state errors
        stateErrorHat = estimatedState[:, 0:6] - trueEphemeris['spacecraft'].T

        plotData = extraData

        plotData['postfit delta']   = extraData['postfit changes']
        plotData['states']          = estimatedState
        plotData['truth']           = dataObservations['truth']
        plotData['beacon_list']     = dataObservations['beacons']
        plotData['timeSpan']        = timeSpan
        plotData['dirIt']           = dirIt
        plotData['err']             = stateError
        plotData['stateErrorHat']   = stateErrorHat
        plotData['obs_uncertainty'] = observationUncertainty
        plotData['referenceState']  = referenceState
        plotData['trueEphemeris']   = trueEphemeris
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

    [anomaly_bool , anomaly_num] = extraData['anomaly_detected']
    if anomaly_bool == True:
        print '**********************************************************'
        print 'Anomaly Detected - Estimates are not to be trusted'
        print '**********************************************************'
        print anomaly_num, 'Residuals out of bounds'

    return


if __name__ == "__main__":
    main()
