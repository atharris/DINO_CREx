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
import matplotlib.pyplot as plt

from rngRngRt import fncH
from rngRngRt import fncG

from posVelCr import EOM
from posVelCr import matrixA

import numpy as np

import TurboProp

from numpy.linalg import inv as aInv
from beaconBinSPICE import getObs


################################################################################
#                  E X P O R T E D     F U N C T I O N S:
################################################################################

# -------------------------------------------------------------------------------
def runRef(input):
    IC0 = input[0]
    phi0 = input[1]
    t_span = input[2]
    extras = input[-1]

    # size of estimation state
    n_state = IC0.shape[0]

    # get the primary body index
    primary_index = extras['primary']
    # secondary body indices
    secondary_indices = extras['secondary']
    # number of secondary bodies
    n_secondaries = len(secondary_indices)

    # list of gravitational parameters (GP)
    mu = extras['mu']
    # GP of the primary body
    mu_primary = mu[primary_index]
    # GP of secondaries
    mu_secondaries = [mu[ii] for ii in secondary_indices]

    # SRP index
    kSRP = extras['SRP']

    # coefficient of reflectivity
    cR = extras['cR']

    # body I.D. list
    bodies = extras['bodies']

    # abberation correction for spkzr
    abcorr = extras['abcorr']
    # reference frame
    ref_frame = extras['ref_frame']

    # organize args
    args = (
    primary_index, secondary_indices, n_secondaries, mu_primary, 
    mu_secondaries, kSRP, cR, abcorr, ref_frame, bodies, n_state)

    # resize the STM to fit in 1D
    phi_append = np.resize(phi0, n_state ** 2)

    # append STM to IC
    IC = np.append(IC0, phi_append)

    # propagate the IC and STM
    state = integ.odeint(EOM, IC, t_span, args)

    return state


def run_CKF(input):
    IC0        = input[0]
    phi0       = input[1]
    time       = input[2]
    SPICE_data = input[3]
    barP0      = input[4]
    observation_uncertainty = input[5]
    extras     = input[-1]

    # number of estimated states
    n_state = IC0.shape[0]

    # number of samples/observations
    n_samples = time.shape[0]

    # create very large state array for the reference as well as estimated states
    # this includes the STM values reshaped into 1D arrays
    ref_state = np.zeros((n_samples, n_state ** 2 + n_state))
    est_state = np.zeros((n_samples, n_state ** 2 + n_state))

    # add the IC to both the est and reference state. NOTE: the first observation
    # will not be computed as it is considered t0
    ref_state[0, :] = np.append(IC0, np.resize(phi0, (1, n_state ** 2)))
    est_state[0, :] = np.append(IC0, np.resize(phi0, (1, n_state ** 2)))

    # initiate x_hat and an array to store these values
    x_hat = np.zeros((n_state, 1))
    x_hat_array = np.zeros((n_samples, n_state))

    # initiate the covariance and covariance array to store these values
    P_array = np.zeros((n_samples, n_state, n_state))
    P_array[0, :, :] = np.copy(barP0)

    # get the observations
    #
    # observation inputs
    obs_inputs = (SPICE_data, observation_uncertainty, extras)
    # get observations and the associated beacon ket
    Y_obs = getObs(obs_inputs)

    # create array for observation differenes to be stored in the filter
    y_array  = np.zeros(Y_obs['data'].shape)
    postfits = np.zeros(Y_obs['data'].shape)

    # inverse of weighting matrix
    # R = aInv( observation_uncertainty )

    # loop through the observations. the -1 is because the first observation is at t0, but the
    # propagator starts at t0 to t1
    for tt in xrange(n_samples - 1):
        # SPICE data dict for input of data into estimated observation calculations
        SPICE_data_GH = {}
        # initial conditions of t0, where t0 = t_(tt-1)
        IC = np.copy( ref_state[ tt, 0:n_state ] )
        # initial condition of STM is identity
        phi = np.identity( n_state )

        # t_span is from current time step to next
        t_span = np.copy( time[ tt:(tt + 2) ] )
        print tt

        # for each beacon, retrieve the state at the filter considered time ( t + 1 )
        for ii in xrange( len( Y_obs['beacons'][tt+1] ) ):
            extras['obs_beacons'] = list(Y_obs['beacons'][tt+1])
            key = extras['obs_beacons'][ii] 
            SPICE_data_GH[key] = np.expand_dims(SPICE_data[key][:, tt + 1], axis=1)

        # input to the propagator takes the state and STM at t0, as well as the list of times
        prop_input = (IC, phi, t_span, extras)

        # execute reference propagation
        state = runRef(prop_input)

        # save the state into a very large array for later
        ref_state[tt + 1, :] = np.copy(state[-1, :])
        est_state[tt + 1, :] = np.copy(state[-1, :])

        # pull out the STM
        phi_t_t0 = np.reshape(state[-1, n_state:], (n_state, n_state))

        # linearly propagate x_hat
        x_bar = np.dot(phi_t_t0, x_hat)

        # linearly propagate covariance
        P_bar = np.dot(np.dot(phi_t_t0, P_array[tt,:,:]), phi_t_t0.T)
        
        if 'SNC' in extras :
           delta_t = t_span[-1] - t_span[0]
           gamma = np.identity(n_state)*delta_t**2
           # P_bar += np.dot( gamma * extras['SNC'], gamma.T )
           SNC = np.identity( n_state ) 
           SNC[0:3,0:3] = np.identity( 3 ) * 10
           SNC[3:6,3:6] = np.identity( 3 ) * .001

        # inputs for Y_refs (G) calculation include the ref state at the lastest time and
        # the contemporary SPICE data
        G_ref_inputs = (np.expand_dims(state[-1, 0:n_state], axis=0), SPICE_data_GH, extras)

        # calculate the estimated observables and organize into an array
        Y_refs = fncG(G_ref_inputs)

        # using the inputs of G, calculate the H matrix
        H_inputs = (np.expand_dims(state[-1, 0:n_state], axis=0), SPICE_data_GH, extras)
        H_tilde = fncH(H_inputs)

        # calculate the deviation of the observables ( Y - G )
        y = Y_obs['data'][tt+1, :] - Y_refs
        y_array[tt+1, :] = np.copy( y )

        KR = aInv(np.dot(np.dot(H_tilde, P_bar), H_tilde.T))
        KL = np.dot(P_bar, H_tilde.T)
        K  = np.dot(KL, KR)

        # perform least squares on the info and observation matrix to compute the residuals

        x_hat                  = x_bar + np.dot(K, (y.T - np.dot(H_tilde, x_bar)))
        x_hat_array[tt + 1, :] = x_hat.T

        postfit           = (y.T - np.dot(H_tilde, x_hat))
        postfits[tt+1, :] = postfit.T

        P_array[tt+1, :, :] = np.dot( np.dot( np.identity(n_state)-np.dot(K,H_tilde), P_bar ), 
                               ( np.identity(n_state)-np.dot(K,H_tilde) ).T ) + np.dot(K,K.T)
        est_state[tt + 1, 0:n_state] += np.squeeze(x_hat)

    extra_data = {}
    extra_data['P_array'] = P_array
    extra_data['x_hat_array'] = x_hat_array
    extra_data['prefit residuals'] = y_array
    extra_data['postfit residuals'] = postfits
    extra_data['Y'] = Y_obs

    return ref_state, est_state, extra_data


def run_EKF(input, tt_switch):
    IC0 = input[0]
    phi0 = input[1]
    time = input[2]
    SPICE_data = input[3]
    barP0 = input[4]
    observation_uncertainty = input[5]
    extras = input[-1]

    # number of estimated states
    n_state = IC0.shape[0]

    # number of samples/observations
    n_samples = time.shape[0]

    # create very large state array for the reference as well as estimated states
    # this includes the STM values reshaped into 1D arrays
    ref_state = np.zeros((n_samples, n_state ** 2 + n_state))
    est_state = np.zeros((n_samples, n_state ** 2 + n_state))

    # copy the STM for the IC
    phi = np.resize(phi0, (1, n_state ** 2))

    # append STM to state
    IC = np.append(IC0, phi).flatten()

    # add the IC to both the est and reference state. NOTE: the first observation
    # will not be computed as it is considered t0
    ref_state[0, :] = IC
    est_state[0, :] = IC

    # initiate x_hat and an array to store these values
    x_hat = np.zeros((n_state, 1))
    x_hat_array = np.zeros((n_samples, n_state))

    # initiate the covariance and covariance array to store these values
    P = np.copy(barP0)
    # P_array = np.zeros((n_samples, n_state ** 2))
    P_array = np.zeros((n_samples, n_state, n_state))

    # P_array[0, :] = np.resize(P_bar, (1, n_state ** 2))
    P_array[0, :, :] = P
    # get the observations
    #
    # observation inputs
    obs_inputs = (SPICE_data, observation_uncertainty, extras)
    # get observations and associated beacon keys
    Y_obs = getObs(obs_inputs)
    print Y_obs['data'].shape

    # create array for observation differenes to be stored in the filter
    y_array = np.zeros(Y_obs['data'].shape)
    postfits = np.zeros(Y_obs['data'].shape)

    # loop through the observations. the -1 is because the first observation is at t0, but the
    # propagator starts at t0 to t1
    for tt in xrange(n_samples - 1) :
        # SPICE data dict for input of data into estimated observation calculations
        SPICE_data_GH = {}
        # initial conditions of t0, where t0 = t_(tt-1)
        if tt > tt_switch:
            pdb.set_trace()
            IC = np.copy(ref_state[tt, 0:n_state] + x_hat.T).flatten()
        else:
            IC = np.copy(ref_state[tt, 0:n_state])
        # initial condition of STM is identity
        phi = np.identity(n_state)
        # t_span is from current time step to next
        t_span = np.copy(time[tt:(tt + 2)])
        print tt

        # for each beacon, retrieve the state at the filter considered time ( t + 1 )
        for ii in xrange( len( Y_obs['beacons'][tt+1] ) ):
            extras['obs_beacons'] =  Y_obs['beacons'][tt+1]
            key = extras['obs_beacons'][ii]
            SPICE_data_GH[key] = np.expand_dims(SPICE_data[key][:, tt + 1], axis=1)

        # input to the propagator takes the state and STM at t0, as well as the list of times
        prop_input = (IC, phi, t_span, extras)

        # execute reference propagation
        state = runRef(prop_input)

        # save the state into a very large array for later
        ref_state[tt + 1, :] = np.copy(state[-1, :])
        est_state[tt + 1, :] = np.copy(state[-1, :])

        # pull out the STM
        phi_step = np.reshape(state[-1, n_state:], (n_state, n_state))

        # linearly propagate x_hat
        if tt > tt_switch:
            x_bar = np.zeros(x_hat.shape)
        else:
            x_bar = np.dot(phi_step, x_hat)

        # linearly propagate covariance
        P_bar = np.dot(np.dot(phi_step, P), phi_step.T)

        # inputs for Y_refs (G) calculation include the ref state at the lastest time and
        # the contemporary SPICE data
        G_ref_inputs = (np.expand_dims(state[-1, 0:n_state], axis=0), SPICE_data_GH, extras)

        # calculate the estimated observables and organize into an array
        Y_refs = fncG(G_ref_inputs)

        # using the inputs of G, calculate the H matrix
        H_inputs = (np.expand_dims(state[-1, 0:n_state], axis=0), SPICE_data_GH, extras)

        H_tilde = fncH(H_inputs)

        # calculate the deviation of the observables ( Y - G )
        y = Y_obs['data'][tt+1, :] - Y_refs
        # print y
        y_array[tt+1, :] = y

        KR = aInv(np.dot(np.dot(H_tilde, P_bar), H_tilde.T))
        KL = np.dot(P_bar, H_tilde.T)
        K  = np.dot(KL, KR)

        # print K

        # perform least squares on the info and observation matrix to compute the residuals

        if tt > tt_switch:
            x_hat = np.dot(K, y.T)
        else:
            x_hat = x_bar + np.dot(K, (y.T - np.dot(H_tilde, x_bar)))

        # print x_hat
        # raw_input()

        x_hat_array[tt + 1, :] = x_hat.T

        postfit = (y.T - np.dot(H_tilde, x_hat))
        postfits[tt+1, :] = postfit.T

        P = np.dot((np.identity(n_state) - np.dot(K, H_tilde)), P_bar)
        # P_array[tt + 1, :] = np.resize(P, (1, n_state ** 2))
        P_array[tt+1, :, :] = P
        est_state[tt + 1, 0:n_state] += np.squeeze(x_hat)

    extra_data = {}
    extra_data['P_array'] = P_array
    extra_data['x_hat_array'] = x_hat_array
    extra_data['prefit residuals'] = y_array
    extra_data['postfit residuals'] = postfits
    extra_data['Y'] = Y_obs

    return ref_state, est_state, extra_data

def run_DBG(input):
    IC0        = input[0]
    phi0       = input[1]
    time       = input[2]
    SPICE_data = input[3]
    barP0      = input[4]
    observation_uncertainty = input[5]
    extras     = input[-1]

    # number of estimated states
    n_state = IC0.shape[0]

    # number of samples/observations
    n_samples = time.shape[0]

    # create very large state array for the reference as well as estimated states
    # this includes the STM values reshaped into 1D arrays
    ref_state = np.zeros((n_samples, n_state ** 2 + n_state))
    est_state = np.zeros((n_samples, n_state ** 2 + n_state))

    # copy the STM for the IC
    phi = np.resize(phi0, (1, n_state ** 2))

    # append STM to state
    IC = np.append(IC0, phi)

    # add the IC to both the est and reference state. NOTE: the first observation
    # will not be computed as it is considered t0
    ref_state[0, :] = IC
    est_state[0, :] = IC

    # loop through the observations. the -1 is because the first observation is at t0, but the
    # propagator starts at t0 to t1
    for tt in xrange(n_samples - 1):
        # SPICE data dict for input of data into estimated observation calculations
        SPICE_data_GH = {}
        # initial conditions of t0, where t0 = t_(tt-1)
        IC = np.copy(ref_state[tt, 0:n_state])
        # initial condition of STM is identity
        phi = np.identity(n_state)
        # t_span is from current time step to next
        t_span = np.copy(time[tt:(tt + 2)])
        print tt

        # for each beacon, retrieve the state at the filter considered time ( t + 1 )
        for ii in xrange( len( Y_obs['beacons'][tt+1] ) ):
            extras['obs_beacons'] =  Y_obs['beacons'][tt+1]
            key = extras['obs_beacons'][ii]
            SPICE_data_GH[key] = np.expand_dims(SPICE_data[key][:, tt + 1], axis=1)

        # input to the propagator takes the state and STM at t0, as well as the list of times
        prop_input = (IC, phi, t_span, extras)

        # execute reference propagation
        state = runRef(prop_input)

        # save the state into a very large array for later
        ref_state[tt + 1, :] = np.copy(state[-1, :])
        est_state[tt + 1, :] = np.copy(state[-1, :])

    extra_data = {}

    return ref_state, est_state, extra_data



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
