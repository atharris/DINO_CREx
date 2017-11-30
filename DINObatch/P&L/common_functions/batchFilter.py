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

from posVel import EOM
from posVel import matrixA

from beaconPropagator import beaconStates

from numpy.linalg import inv as aInv

import pdb

################################################################################
#                  E X P O R T E D     F U N C T I O N S:
################################################################################

#-------------------------------------------------------------------------------
def norm( input ) :
  norm = np.sqrt( sum( np.square( input ) ) )
  return norm

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
  extras       = input[-1]

  # number of estimated states
  stateDimension   = IC.shape[0]

  # number of samples/observations
  nObservations = timeSpan.shape[0]

  # observed measurement data
  observations = filterObservations['measurements']

  # initiate stateErrorHat
  stateErrorHat = np.zeros( (stateDimension,) )

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
  referenceObservationInputs = ( referenceState[:,0:stateDimension], beaconStateArray, extras )

  # calculate the estimated observables and organize into an array
  referenceObservations = fncG( referenceObservationInputs )


  # using the inputs of G, calculate the H matrix
  mappingMatrixInputs  = ( referenceState[:,0:stateDimension], beaconStateArray, extras )
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
  stateErrorHat = np.reshape(np.linalg.lstsq( infoMatrix, normalMatrix )[0], [stateDimension])

  # initiate a filtered referenceState
  estimatedState = np.zeros( (referenceState.shape[0],stateDimension) )

  # initiate an array to hold the filtered covariances
  # covArray = np.zeros( (nObservations, stateDimension, stateDimension) )

  # the first filtered covariance is the inverse of the covariance matrix
  P = aInv( infoMatrix )

  # initiate an array for the state deviation vectors
  stateErrorHatArray = np.zeros( (nObservations, stateDimension) )
  stateErrorBarArray = np.zeros( (nObservations, stateDimension) )

  for ii in xrange( referenceState.shape[0] ) :
    # pull the STM that is able to transform from the current time to t0
    phi_t_t0  = np.reshape( referenceState[ii,stateDimension:],\
                            (stateDimension,stateDimension) )
    # linearly transform deviation and at it to the referenceState and save
    stateErrorHatArray[ii,:] = np.dot( phi_t_t0, stateErrorHat ).T
    stateErrorBarArray[ii,:] = np.dot( phi_t_t0, x_bar ).T
    covArray[ii,:,:] = np.dot(np.dot( phi_t_t0, P ),  phi_t_t0.T)
    # add the deviation to the reference state 
    estimatedState[ii,:]   = referenceState[ii,0:stateDimension] + stateErrorHatArray[ii,:]
    # store the transformed covariance matrix
    # covArray[ii,:,:]   = np.dot( np.dot( phi_t_t0, P ), phi_t_t0.T )

  # inputs for estimatedObservations (G) calculation
  estimatedObservationInputs = ( estimatedState[:,0:stateDimension], beaconStateArray, extras )

  # calculate the estimated observables and organize into an array
  estimatedObservations = fncG( estimatedObservationInputs )

  # compute the postfits using the updated observables and the measured values
  postfits      = np.zeros([np.shape(stateErrorHatArray)[0], \
                            np.shape(observationDeviations)[1]])
  postfitsDelta = np.zeros([np.shape(stateErrorHatArray)[0], \
                            np.shape(observationDeviations)[1]])

  # Anomaly detection
  for ii in range(np.shape(x_hat_array)[0]):
    for jj in range(np.shape(postfits)[1]):
      if np.abs(postfits[ii,jj]) - 3*observation_uncertainty[jj, jj] > 0:
        extras['anomaly_num']+=1
        print 'Anomalous measurement detected at time ' , ii , 'on measurement type ', jj

  if extras['anomaly_num'] > extras['anomaly_threshold']:
    extras['anomaly'] = True

  prefits = np.zeros([np.shape(x_hat_array)[0], np.shape(y)[1]])
  for ii in range(1,np.shape(x_hat_array)[0]):
    prefits[ii,:] = y[ii,:] - np.dot(H_tilde[0+2*(ii):2+2*(ii),:], x_bar_array[ii,:])



  prefits = np.zeros([np.shape(stateErrorHatArray)[0], np.shape(observationDeviations)[1]])
  for ii in range(1,np.shape(stateErrorHatArray)[0]):
    prefits[ii,:] = observationDeviations[ii,:] - np.dot(mappingMatrix[0+2*(ii):2+2*(ii),:],\
                    stateErrorBarArray[ii,:])

  # store various arrays in a data dictionary
  extra_data                      = {}
  extra_data['Y']                 = Y_obs
  extra_data['P_array']           = P_array
  extra_data['x_hat_array']       = x_hat_array
  extra_data['prefit residuals']  = prefits
  extra_data['postfit residuals'] = postfits
  extra_data['postfit changes'] = postfitsDelta
  extras['x_hat_0']               += x_hat
  extra_data['x_hat_0']           = extras['x_hat_0']
  extra_data['anomaly_detected']  = [extras['anomaly'], extras['anomaly_num']]


  return referenceState, estimatedState, extraData


################################################################################
#                       M A I N     F U N C T I O N:
################################################################################

def main():

  return

if __name__ == "__main__":
  main()
