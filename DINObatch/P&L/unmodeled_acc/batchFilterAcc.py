#!/usr/local/bin/python
'''

 <<Description>>


 <<Summary>>

'''

__author__ = 'Marc Balducci'
__version__ = '$Revision$'[11:-2]
__date__ = '$Date$'[7:26]

################################################################################
#                     I M primaryBodyIndex O R T     L I B R A R I E secondaryBodyIndices
################################################################################

import scipy.integrate as integ
import numpy as np

from pixelLineBatch import fncH
from pixelLineBatch import fncG

from posVelAcc import EOM
from posVelAcc import matrixA

from beaconPropagator import beaconStates

from numpy.linalg import inv as aInv

import pdb

## \defgroup batch_acc Batch Filter - unmodeled acceleration
##   @{
## The module for the unmodeled acceleration estimation batch filter.
#
# Overview {#overview}
# ====
#
# Purpose
# -----
# This script contains two functions that may be exported for the purposes of reference trajectory propagation and the running of a batch filter. 
#
# It is noted here and in latter sections of this document, this module is largely identical to that of the vanilla batch filter \ref batch_vanilla "`batchFilter.py`". The exception is the equations of motions (EOMs) that are imported. Rather than utilizing \ref EOMS_vanilla "`posVel.py`", the EOM function containing support for a constant acceleration in the x, y and z coordinates is imported as \ref EOMs_acc "`posVelAcc.py`". 
#
# A significant amount of material originiates from Statistical Orbit Determination (2004) a work by Tapley, Schutz, and Born. Therefore, equation numbers from the text will be included with the theory and presentation. The purpose of this documentation is to provide a succinct foundation of theory as well as a walkthrough of the overall code structure.
#
# Contents
# -----
# The following functions are contained in this module:
#
# - `runRef.py`
# - `run_batch.py`
#
# The former (`runRef.py`) is an exportable function that parses inputs and calls an ODE solver for specified equations of motion (EOMs). The outputs are an integrated state and state transition matrix (STM). 
#
# The latter function (`run_batch.py`) is also an exportable function. It organizes inputs, calls the reference propagator runRef.py, runs beaconStates to create reference beacon positions, creates reference observations and associated H matrix by calling fncG and fncH, respectively, and computes deviations using a batch algorithm.
#
# Neither function is a stand alone script. As with other modules in the state estimation nav filter, there is a reliance on the `extras` dictionary to pass through parameters to various functions. It is noted that the `extras` dictionary should never be an output from a function.
#
# The Code
# =====
#
# `runRef.py` 
# -----
# The purpose of this function is relatively straightforward: calculate a reference trajectory when given an initial condition, a set of times, and a dictionary of extras. The following is a table of inputs and associated brief descriptions:
#
# Name      | Description                                           | Size/Type
# -----     | -------------------                                   | -----
# IC0       | initial conditions of state                           | (d,N) numpy array  
# phi0      | initial condition of STM                              | (d,d) numpy array
# timeSpan  | an array of times to integrate to                     | (N,) numpy array
# extras    | dictionary of various parameters                      | dictionary
#
# The same is repeated for outputs:
#
# Name      | Description         | Size/Type
# -----     | ------------------- | -----
# state     | propagated state    | (N,d(1+d)) numpy array
#
# Thusly, we have established the purpose of this function. It is noted that the output variable `state` contains the quantities of interest for the scenario (position, velocity, etc) as well as the propagates STM that has been resized to (,dXd). 
#
# `run_batch.py`
# -----
# This function consists of the meat of the batch filter algorithm. As previously described, it contains and calls a variety of calculations. To begin, we once again list inputs: 
#
# Name      | Description                                           | Size/Type
# -----     | -------------------                                   | -----
# IC        | initial conditions of state                           | (d,) numpy array  
# phi       | initial condition of STM                              | (d,d) numpy array
# timeSpan  | an array of times to integrate to                     | (N,) numpy array
# filterObservations | dictionary of observation related data     | dictionary
# P_bar     | initial covariance                                    | (d,d) numpy array
# observationUncertainty | a priori uncertainty of measurements. diagonal array | (m,m) numpy array
# x_bar     | a priori state deviation                              | (d,) numpy array  
# extras    | dictionary of various parameters                      | dictionary
#
# The same is repeated for outputs:
#
# Name      | Description         | Size/Type
# -----     | ------------------- | -----
# referenceState | propagated reference state | (N,d(1+d)) numpy array
# estimatedState | propagated reference state + estimated state deviations | (N,d) numpy array
# extraData | dictionary of various outputs, e.g., deviations | dictionary
#
# As previously stated, the code contained within this module is identical to that of \ref batch_vanilla "`batchFilter.py`" with the exception of the imported EOM functions. Rather than utilizing \ref EOMS_vanilla "`posVel.py`", the EOM function containing support for a constant acceleration in the x, y and z coordinates is imported as \ref EOMs_acc "`posVelAcc.py`". 
## @}

################################################################################
#                  E X P O R T E D     F U N C T I O N S:
################################################################################

#-------------------------------------------------------------------------------

def runRef( input ) :

  IC0          = input[0]
  phi0         = input[1]
  timeSpan       = input[2]
  extras       = input[-1]

  # size of estimation state
  stateDimension           = IC0.shape[0]

  # get the primary body index
  primaryBodyIndex     = extras['primary']
  # secondary body indices
  secondaryBodyIndices = extras['secondary']
  # number of secondary bodies
  nSecondaries     = len( secondaryBodyIndices )
  
  # list of gravitational parameters (GP)
  mu                = extras['mu']
  # GP of the primary body
  muPrimary        = mu[primaryBodyIndex]
  # GP of secondaries
  muSecondaries    = [ mu[ii] for ii in secondaryBodyIndices ] 
  
  # SRP index
  kSRP   = extras['SRP']
  # coefficient of reflectivity
  cR     = extras['cR']

  # body I.D. list
  bodies = extras['bodies']
  
  # abberation correction for spkzr
  abcorr = extras['abcorr']
  # reference frame
  referenceFrame    = extras['ref_frame']

  # organize args
  args   = ( primaryBodyIndex, secondaryBodyIndices, nSecondaries, muPrimary, muSecondaries, kSRP, cR, abcorr, referenceFrame, bodies, stateDimension )
  
  # reshape the STM to fit in 1D
  phi = np.reshape(phi0, stateDimension**2 )

  # append STM to IC
  IC     = np.append( IC0, phi )

  # propagate the IC and STM 
  state  = integ.odeint( EOM, IC, timeSpan, args )

  return state 

def run_batch( input ) :

  IC           = input[0]
  phi          = input[1]
  timeSpan       = input[2]
  filterObservations   = input[3]
  P_bar        = input[4]
  observationUncertainty = input[5]
  x_bar        = input[6]
  angles       = input[7]
  extras       = input[-1]

  # number of estimated states
  stateDimension   = IC.shape[0]

  # number of samples/observations
  nObservations = timeSpan.shape[0]

  # observed measurement data
  observations = filterObservations['measurements']

  # initiate stateDevHat
  stateDevHat = np.zeros( (stateDimension,) )

  IC0  = np.copy( IC )


  ##################################################################################
  #
  # Integrate Reference Trajectory page 196
  #
  ##################################################################################
 
  # if there exists a priori information. if not, use default
  if np.sum(np.sum( P_bar )) == 0 :
    infoMatrix = np.zeros( (stateDimension,stateDimension) )
    normalMatrix = np.zeros( (stateDimension,1) )
 
  else :
    infoMatrix = aInv( P_bar )
    normalMatrix = np.dot( aInv( P_bar ), 
                    np.expand_dims(x_bar,axis=1) )

  phi0       = np.identity( stateDimension )

  # input to the propagator takes the referenceState and STM at t0, as well as the list of times
  propagatorInput = ( IC0, phi0, timeSpan, extras ) 
  
  # execute propagation
  referenceState  = runRef( propagatorInput )
  
    #######################################################################
    # The measurements are currently of size 2 for each time step. 
    # This is hard coded in a few places with the placeholder int 2. 
    #######################################################################

  ##################################################################################
  #
  # Calculate Observation Block - aka accumulate current observation at all times page 196
  #
  ##################################################################################

  # collect the list of beacons (in observational order) into an extras list
  extras['obs_beacons'] = list(filterObservations['beaconIDs'])

  # generate the positions of beacons using dictated function beaconPositions
  beaconPropInputs = ( filterObservations['beaconIDs'], timeSpan, extras )
  beaconStateArray = beaconStates( beaconPropInputs )

  # create observation weight matrix (W)
  W = aInv( observationUncertainty )

  # inputs for referenceObservations (G) calculation
  referenceObservationInputs = ( referenceState[:,0:stateDimension], beaconStateArray, angles, extras )

  # calculate the estimated observables and organize into an array
  referenceObservations = fncG( referenceObservationInputs )


  # using the inputs of G, calculate the H matrix
  mappingMatrixInputs  = ( referenceState[:,0:stateDimension], beaconStateArray, angles, extras )
  mappingMatrix        = fncH( mappingMatrixInputs )


  # calculate the deviation of the observables ( Y - G )
  observationDeviations = observations - referenceObservations

  # initiate an array to hold the filtered covariances
  covArray = np.zeros( (nObservations, stateDimension, stateDimension) )
  # cycle through the observations/reference states and build up the filter data
  for ii in xrange(nObservations) :
      # pull out the STM
      phi_t_t0  = np.reshape( referenceState[ii,stateDimension:],\
                             (stateDimension,stateDimension) )
      # matrix multiply the H matrix at time tii with that of the contemporary STM
      # RANGE AND RANGE RATE BATCH - hard coded 2 for observation type size
      H    = np.dot( mappingMatrix[0+2*ii:2+2*ii,:], phi_t_t0 )
      # add the new H^T H result to the information matrix
      infoMatrix   += np.dot( H.T, np.dot( W, H ) )
      # covArray[ii,:,:]   = aInv( infoMatrix )
      # add the H^T Y result to the observation information matrix
      yii = np.zeros([len( observationDeviations[ii,:]),1])
      yii[:,0] = observationDeviations[ii,:]
      normalMatrix += np.dot( H.T, np.dot( W, yii))

  ##################################################################################
  #
  # \ Calculate Observation Block - aka accumulate current observation at all times page 196
  #
  ##################################################################################
    
  # perform least squares on the infoMatrix and observation matrix to compute the residuals
  stateDevHat = np.reshape(np.linalg.lstsq( infoMatrix, normalMatrix )[0], [stateDimension])

  # initiate a filtered referenceState
  estimatedState = np.zeros( (referenceState.shape[0],stateDimension) )

  # initiate an array to hold the filtered covariances
  # covArray = np.zeros( (nObservations, stateDimension, stateDimension) )

  # the first filtered covariance is the inverse of the covariance matrix
  P = aInv( infoMatrix )

  # initiate an array for the state deviation vectors
  stateDevHatArray = np.zeros( (nObservations, stateDimension) )
  stateDevBarArray = np.zeros( (nObservations, stateDimension) )

  for ii in xrange( referenceState.shape[0] ) :
    # pull the STM that is able to transform from the current time to t0
    phi_t_t0  = np.reshape( referenceState[ii,stateDimension:],\
                            (stateDimension,stateDimension) )
    # linearly transform deviation and at it to the referenceState and save
    stateDevHatArray[ii,:] = np.dot( phi_t_t0, stateDevHat ).T
    stateDevBarArray[ii,:] = np.dot( phi_t_t0, x_bar ).T
    covArray[ii,:,:] = np.dot(np.dot( phi_t_t0, P ),  phi_t_t0.T)
    # add the deviation to the reference state 
    estimatedState[ii,:]   = referenceState[ii,0:stateDimension] + stateDevHatArray[ii,:]
    # store the transformed covariance matrix
    # covArray[ii,:,:]   = np.dot( np.dot( phi_t_t0, P ), phi_t_t0.T )

  # inputs for estimatedObservations (G) calculation
  estimatedObservationInputs = ( estimatedState[:,0:stateDimension], beaconStateArray, angles, extras )

  # calculate the estimated observables and organize into an array
  estimatedObservations = fncG( estimatedObservationInputs )

  # compute the postfits using the updated observables and the measured values
  postfits      = np.zeros([np.shape(stateDevHatArray)[0], \
                            np.shape(observationDeviations)[1]])
  postfitsDelta = np.zeros([np.shape(stateDevHatArray)[0], \
                            np.shape(observationDeviations)[1]])

  for ii in range(np.shape(stateDevHatArray)[0]):
    postfitsDelta[ii,:] = extras['oldPost'][ii,:] - observationDeviations[ii,:] + np.dot(mappingMatrix[0+2*ii:2+2*ii,:], stateDevHatArray[ii,0:np.shape(mappingMatrix)[1]])
    postfits[ii,:] = observationDeviations[ii,:] - np.dot(mappingMatrix[0+2*ii:2+2*ii,:], stateDevHatArray[ii,:])

  # Anomaly detection
  for ii in range(np.shape(stateDevHatArray)[0]):
    for jj in range(np.shape(postfits)[1]):
      if np.abs(postfits[ii,jj]) - 3*observationUncertainty[jj, jj] > 0:
        extras['anomaly_num']+=1
        print 'Anomalous measurement detected at time ' , ii , 'on measurement type ', jj

  if extras['anomaly_num'] > extras['anomaly_threshold']:
    extras['anomaly'] = True

  prefits = np.zeros([np.shape(stateDevHatArray)[0], np.shape(observationDeviations)[1]])
  for ii in range(1,np.shape(stateDevHatArray)[0]):
    prefits[ii,:] = observationDeviations[ii,:] - np.dot(mappingMatrix[0+2*(ii):2+2*(ii),:], stateDevBarArray[ii,:])

  # store various arrays in a data dictionary
  extraData                       = {}
  extraData['Y']                  = observations
  extraData['covArray']           = covArray
  extraData['stateDevHatArray']   = stateDevHatArray
  extraData['prefit residuals']   = prefits
  extraData['postfit residuals']  = postfits
  extraData['postfit changes']    = postfitsDelta
  extraData['stateDevHat']        = extras['x_hat_0']
  extraData['anomaly_detected']   = [extras['anomaly'], extras['anomaly_num']]


  return referenceState, estimatedState, extraData


################################################################################
#                       M A I N     F U N C T I O N:
################################################################################

def main():

  return

if __name__ == "__main__":
  main()
