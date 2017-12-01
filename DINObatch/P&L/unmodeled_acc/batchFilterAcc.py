#!/usr/local/bin/python
'''

 <<Description>>


 <<Summary>>

'''

__author__ = 'Marc Balducci'
__version__ = '$Revision$'[11:-2]
__date__ = '$Date$'[7:26]

################################################################################
#                     I M primary_index O R T     L I B R A R I E secondary_indices
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
  t_span       = input[2]
  extras       = input[-1]

  # size of estimation state
  n_state           = IC0.shape[0]

  # get the primary body index
  primary_index     = extras['primary']
  # secondary body indices
  secondary_indices = extras['secondary']
  # number of secondary bodies
  n_secondaries     = len( secondary_indices )
  
  # list of gravitational parameters (GP)
  mu                = extras['mu']
  # GP of the primary body
  mu_primary        = mu[primary_index]
  # GP of secondaries
  mu_secondaries    = [ mu[ii] for ii in secondary_indices ] 
  
  # SRP index
  kSRP   = extras['SRP']
  # coefficient of reflectivity
  cR     = extras['cR']

  # body I.D. list
  bodies = extras['bodies']
  
  # abberation correction for spkzr
  abcorr = extras['abcorr']
  # reference frame
  ref_frame    = extras['ref_frame']

  # organize args
  args   = ( primary_index, secondary_indices, n_secondaries, mu_primary, mu_secondaries, kSRP, cR, abcorr, ref_frame, bodies, n_state )
  
  # reshape the STM to fit in 1D
  phi = np.reshape(phi0, n_state**2 )

  # append STM to IC
  IC     = np.append( IC0, phi )

  # propagate the IC and STM 
  state  = integ.odeint( EOM, IC, t_span, args )

  return state 

def run_batch( input ) :

  IC           = input[0]
  phi          = input[1]
  t_span       = input[2]
  obs_filter   = input[3]
  P_bar        = input[4]
  observation_uncertainty = input[5]
  x_bar        = input[6]
  extras       = input[-1]

  # number of estimated states
  n_state   = IC.shape[0]

  # number of samples/observations
  n_samples = t_span.shape[0]

  # observed measurement data
  Y_obs = obs_filter['measurements']

  # initiate x_hat
  x_hat = np.zeros( (n_state,) )

  IC0  = np.copy( IC )


  ##################################################################################
  #
  # Integrate Reference Trajectory page 196
  #
  ##################################################################################
 
  # if there exists a priori information. if not, use default
  if np.sum(np.sum( P_bar )) == 0 :
    info_matrix = np.zeros( (n_state,n_state) )
    normal_matrix = np.zeros( (n_state,1) )
 
  else :
    info_matrix = aInv( P_bar )
    normal_matrix = np.dot( aInv( P_bar ), 
                    np.expand_dims(x_bar,axis=1) )

  phi0       = np.identity( n_state )

  # input to the propagator takes the ref_state and STM at t0, as well as the list of times
  prop_input = ( IC0, phi0, t_span, extras ) 
  
  # execute propagation
  state      = runRef( prop_input )
  ref_state  = np.copy( state )
  
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
  extras['obs_beacons'] = list(obs_filter['beaconIDs'])

  # generate the positions of beacons using dictated function beaconPositions
  beacon_inputs = ( obs_filter['beaconIDs'], t_span, extras )
  beacon_states = beaconStates( beacon_inputs )

  # create observation weight matrix (W)
  W = aInv( observation_uncertainty )

  # inputs for Y_refs (G) calculation
  G_ref_inputs = ( ref_state[:,0:n_state], beacon_states, extras )

  # calculate the estimated observables and organize into an array
  Y_refs = fncG( G_ref_inputs )


  # using the inputs of G, calculate the H matrix
  H_inputs = ( ref_state[:,0:n_state], beacon_states, extras )
  H_tilde   = fncH( H_inputs )


  # calculate the deviation of the observables ( Y - G )
  y    = Y_obs - Y_refs

  # initiate an array to hold the filtered covariances
  P_array = np.zeros( (n_samples, n_state, n_state) )
  # cycle through the observations/reference states and build up the filter data
  for ii in xrange(n_samples) :
      # pull out the STM
      phi_t_t0  = np.reshape( ref_state[ii,n_state:], (n_state,n_state) )
      # matrix multiply the H matrix at time tii with that of the contemporary STM
      # RANGE AND RANGE RATE BATCH - hard coded 2 for observation type size
      H    = np.dot( H_tilde[0+2*ii:2+2*ii,:], phi_t_t0 )
      # add the new H^T H result to the information matrix
      info_matrix   += np.dot( H.T, np.dot( W, H ) )
      # P_array[ii,:,:]   = aInv( info_matrix )
      # add the H^T Y result to the observation information matrix
      yii = np.zeros([len( y[ii,:]),1])
      yii[:,0] = y[ii,:]
      normal_matrix += np.dot( H.T, np.dot( W, yii))

  ##################################################################################
  #
  # \ Calculate Observation Block - aka accumulate current observation at all times page 196
  #
  ##################################################################################
    
  # perform least squares on the info_matrix and observation matrix to compute the residuals
  x_hat = np.dot(np.linalg.inv(info_matrix), normal_matrix).reshape([len(IC)])

  # initiate a filtered ref_state
  est_state = np.zeros( (ref_state.shape[0],n_state) )

  # initiate an array to hold the filtered covariances
  # P_array = np.zeros( (n_samples, n_state, n_state) )

  # the first filtered covariance is the inverse of the covariance matrix
  P = aInv( info_matrix )

  # initiate an array for the state deviation vectors
  x_hat_array = np.zeros( (n_samples, n_state) )
  x_bar_array = np.zeros( (n_samples, n_state) )

  for ii in xrange( ref_state.shape[0] ) :
    # pull the STM that is able to transform from the current time to t0
    phi_t_t0  = np.reshape( ref_state[ii,n_state:], (n_state,n_state) )
    # linearly transform deviation and at it to the ref_state and save
    x_hat_array[ii,:] = np.dot( phi_t_t0, x_hat ).T
    x_bar_array[ii,:] = np.dot( phi_t_t0, x_bar ).T
    P_array[ii,:,:] = np.dot(np.dot( phi_t_t0, P ),  phi_t_t0.T)
    # add the deviation to the reference state 
    est_state[ii,:]   = ref_state[ii,0:n_state] + x_hat_array[ii,:]
    # store the transformed covariance matrix
    # P_array[ii,:,:]   = np.dot( np.dot( phi_t_t0, P ), phi_t_t0.T )

  # inputs for Y_est (G) calculation
  G_est_inputs = ( est_state[:,0:n_state], beacon_states, extras )

  # calculate the estimated observables and organize into an array
  Y_est = fncG( G_est_inputs )

  # compute the postfits using the updated observables and the measured values
  postfits = np.zeros([np.shape(x_hat_array)[0], np.shape(y)[1]])
  postfitsDelta = np.zeros([np.shape(x_hat_array)[0], np.shape(y)[1]])
  for ii in range(np.shape(x_hat_array)[0]):
    postfitsDelta[ii,:] = extras['oldPost'][ii,:] - y[ii,:] + np.dot(H_tilde[0+2*ii:2+2*ii,:], x_hat_array[ii,0:np.shape(H_tilde)[1]])
    postfits[ii,:] = y[ii,:] - np.dot(H_tilde[0+2*ii:2+2*ii,:], x_hat_array[ii,0:np.shape(H_tilde)[1]])

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

  # store various arrays in a data dictionary
  extra_data                      = {}
  extra_data['Y']                 = Y_obs
  extra_data['P_array']           = P_array
  extra_data['x_hat_array']       = x_hat_array
  extra_data['prefit residuals']  = prefits
  extra_data['postfit residuals'] = postfits
  extra_data['postfit changes'] = postfitsDelta
  # is this legacy code?
  extras['x_hat_0']               += x_hat
  extra_data['x_hat_0']           = extras['x_hat_0']
  extra_data['anomaly_detected']  = [extras['anomaly'], extras['anomaly_num']]



  return ref_state, est_state, extra_data


################################################################################
#                       M A I N     F U N C T I O N:
################################################################################

def main():

  return

if __name__ == "__main__":
  main()