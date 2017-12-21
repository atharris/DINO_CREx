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
import matplotlib.pyplot as plt

filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))

sys.path.append(path)

import numpy as np

import pickle

import pdb

## \defgroup compare_SRP compare_SRP_test - a test of the filter with and without SRP
##   @{
## A script used to test the effectiveness of filters with and without SRP
#
# Overview {#overview}
# ====
#
# Purpose
# -----
# This module is a runable script for the purpose of testing the effectiveness of filters with and without SRP in the EOMs. Four variation are considered: runs of \ref init_vanilla with and without SRP in the EOMs and runs of \ref init_acc with and without SRP in the EOMS. These four variations are then compared
#
# Contents
# -----
# The `compare_SRP_test.py` script contains two secondary functions and a main function:
#
# - `main`
# - `rms`
# - `error_norm`
#
# Here, we note that `compare_SRP_test.py` mainly pulls already computed data from designated folders in the system path. These folders contain the outputs of batches run by \ref init_vanilla , \ref init_SRP_vanilla , \ref init_acc , and \ref init_SRP_acc .
#
# The Code
# =====
#
# `main`
# -----
# Because `compare_SRP_test.py` is a script meant to be run from the terminal, there are no inputs. 
#
# Looping through iterations, the script first creates folders for the comparison figures
# ~~~~~~~~~~~~~~~~{.py}
#    # Iteration Directory
#    dirIt   = 'Batch_Iteration' + str(itr)
#    saveDir = 'comparison_figures/'+dirIt
# ~~~~~~~~~~~~~~~~
#
# The file tags for nominal and SRP test pickle data are then defined. These tags are those used in the various `init` scripts to save the filter results
# ~~~~~~~~~~~~~~~~{.py}
#    # state the paths for nominal and SRP test data
#    nominal_path    = dirIt+'/nominal_data.pkl'
#    SRP_test_path   = dirIt+'/SRP_test_data.pkl'
# ~~~~~~~~~~~~~~~~
#
# Going through each of the four variations of filter results, the script then pulls pickle data, unpacks the necessary data and then calculates the root mean square (RMS) of error for position and velocity components, e.g., 
# ~~~~~~~~~~~~~~~~{.py}
#    # pull out the pickle data of the nominal vanilla run (with SRP)
#    nominal_vanilla_data_path = path + '/vanilla_pl/' + nominal_path
#    nominal_vanilla_file      = open( nominal_vanilla_data_path, 'rb' )
#    nominal_vanilla           = pickle.load( nominal_vanilla_file )
#    nominal_vanilla_file.close()
#
#    # unpack data such as truth ephemeris, estimated states, covariance and time span
#    trueEphemeris = nominal_vanilla['trueEphemeris']
#    est_states    = nominal_vanilla['states']
#    covar_vanilla = np.array([[np.sum(nominal_vanilla['covArray'][0:3,0:3]), 0.], [0.,np.sum(nominal_vanilla['covArray'][3:6,3:6]) ]])
#    timeSpan      = nominal_vanilla['timeSpan']
#    
#    # compute RMS absolute error of position and velocity
#    posErrNormNmnlVanilla, velErrNormNmnlVanilla = \
#        errorNorm( trueEphemeris['spacecraft'], est_states )
# ~~~~~~~~~~~~~~~~
#
#  The next significant block of code contains calculations of differences between methodologies. The first lines in this residual work initialize the array
# ~~~~~~~~~~~~~~~~{.py}
#    # initialize arrays for the storage of absolute and relative residuals. Note that the rows
#    # correspond to position and velocity
#    SRP_test_abs_residual = np.zeros( (2,len(timeSpan)) ) # R[pos,vel] C[length]
#    nominal_abs_residual  = np.zeros( (2,len(timeSpan)) ) # R[pos,vel] C[length]
#    SRP_test_rel_residual = np.zeros( (2,len(timeSpan)) ) # R[pos,vel] C[length]
#    nominal_rel_residual  = np.zeros( (2,len(timeSpan)) ) # R[pos,vel] C[length]
# ~~~~~~~~~~~~~~~~
# The "SRP_test" residuals will compare the differences between the two filters (vanilla and acc) without SRP in the reference EOMs, while the "nominal" residuals track the differences between the two filters with SRP in the EOMs. The computation of these differences are contained in the lines:
# ~~~~~~~~~~~~~~~~{.py}
#    # calculate the relative residuals between the two vanilla methods and the two unmodeled 
#    # acc methods
#    SRP_test_rel_residual[0,:] = np.divide(\
#                          np.array([posErrNormTestVanilla-posErrNormTestUnmodeledAcc]),\
#                          np.array(posErrNormTestVanilla))
#    SRP_test_rel_residual[1,:] = np.divide(\
#                          np.array([velErrNormTestVanilla-velErrNormTestUnmodeledAcc]),\
#                          np.array(velErrNormTestVanilla))
#
#    nominal_rel_residual[0,:] = np.divide(\
#                          np.array([posErrNormNmnlVanilla-posErrNormNmnlUnmodeledAcc]),\
#                          np.array(posErrNormNmnlVanilla))
#    nominal_rel_residual[1,:] = np.divide(\
#                          np.array([velErrNormNmnlVanilla-velErrNormNmnlUnmodeledAcc]),\
#                          np.array(velErrNormNmnlVanilla))
# ~~~~~~~~~~~~~~~~
#
# The remaining code plots out these and other more basic values, as well as writing RMS values to the terminal.
#
# `rms`
# -----
# A simple function used to calculate the root mean square value of an input numpy array
#
# `errorNorm`
# -----
# Here, we give the inputs:
# Name      | Description                                           | Size/Type
# -----     | -------------------                                   | -----
# truth | truth state | (N,d) numpy array  
# estimate | comparison state (most likely an estimate) | (N,d) numpy array
#
# The same is repeated for outputs:
#
# Name      | Description         | Size/Type
# -----     | ------------------- | -----
# posErrorNorm | norm of the position errors | (N,) numpy array  
# velErrorNorm | norm of the velocity errors | (N,) numpy array
#
# The `errorNorm` function computes the difference between two input arrays (truth and estimate) and calculates the norm of each entry in the array.
## @}
################################################################################
#                  S E C O N D A R Y     F U N C T I O N S:
################################################################################

# -------------------------------------------------------------------------------
def rms( array ):
    array_rms = np.sqrt( np.mean( array**2 ) ) 

    return array_rms

def errorNorm( truth, estimate ):
    posErrorNorm = np.zeros(len(estimate[:, 0]))
    velErrorNorm = np.zeros(len(estimate[:, 0]))    

    for ii in range(len(estimate[:,0])):
        posErrorNorm[ii] = np.linalg.norm(truth.T[ii, 0:3] - estimate[ii, 0:3])
        velErrorNorm[ii] = np.linalg.norm(truth.T[ii, 3:6] - estimate[ii, 3:6])   

    return posErrorNorm, velErrorNorm 

################################################################################
#                  P R I M A R Y     F U N C T I O N:
################################################################################

# -------------------------------------------------------------------------------
def main() :

  for ii in xrange(3):
    itr = ii+1
    print 'Iteration ' + str(itr)
    # Iteration Directory
    dirIt   = 'Batch_Iteration' + str(itr)
    saveDir = 'comparison_figures/'+dirIt

    # Make directory for the iterations
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)

    # state the paths for nominal and SRP test data
    nominal_path    = dirIt+'/nominal_data.pkl'
    SRP_test_path   = dirIt+'/SRP_test_data.pkl'

    # pull out the pickle data of the nominal vanilla run (with SRP)
    nominal_vanilla_data_path = path + '/vanilla_pl/' + nominal_path
    nominal_vanilla_file      = open( nominal_vanilla_data_path, 'rb' )
    nominal_vanilla           = pickle.load( nominal_vanilla_file )
    nominal_vanilla_file.close()

    # unpack data such as truth ephemeris, estimated states, covariance and time span
    trueEphemeris = nominal_vanilla['trueEphemeris']
    est_states    = nominal_vanilla['states']
    covar_vanilla = np.array([[np.sum(nominal_vanilla['covArray'][0:3,0:3]), 0.], [0.,np.sum(nominal_vanilla['covArray'][3:6,3:6]) ]])
    timeSpan      = nominal_vanilla['timeSpan']
    
    # compute RMS absolute error of position and velocity
    posErrNormNmnlVanilla, velErrNormNmnlVanilla = \
        errorNorm( trueEphemeris['spacecraft'], est_states )

    # pull out the pickle data of the SRP test vanilla run (without SRP)
    SRP_test_vanilla_data_path = path + '/vanilla_pl/' + SRP_test_path
    SRP_test_vanilla_file      = open( SRP_test_vanilla_data_path, 'rb' )
    SRP_test_vanilla           = pickle.load( SRP_test_vanilla_file )
    SRP_test_vanilla_file.close()

    # unpack data such as truth ephemeris, estimated states, covariance and time span
    trueEphemeris = SRP_test_vanilla['trueEphemeris']
    est_states = SRP_test_vanilla['states']
    covar_vanilla_SRP_test = np.array([[np.sum(SRP_test_vanilla['covArray'][0:3,0:3]), 0.], [0.,np.sum(SRP_test_vanilla['covArray'][3:6,3:6]) ]])
    
    # compute RMS absolute error of position and velocity
    posErrNormTestVanilla, velErrNormTestVanilla = \
        errorNorm( trueEphemeris['spacecraft'], est_states )

    # pull out the pickle data of the nominal unmodeled acceleration run (with SRP)
    nominal_unmodeled_acc_data_path = path + '/unmodeled_acc/' + nominal_path
    nominal_unmodeled_acc_file      = open( nominal_unmodeled_acc_data_path, 'rb' )
    nominal_unmodeled_acc           = pickle.load( nominal_unmodeled_acc_file )
    nominal_unmodeled_acc_file.close()

    # unpack data such as truth ephemeris, estimated states, covariance and time span
    trueEphemeris = nominal_unmodeled_acc['trueEphemeris']
    est_states = nominal_unmodeled_acc['states']
    covar_unmodeledAcc = np.array([[np.sum(nominal_unmodeled_acc['covArray'][0:3,0:3]), 0.], [0.,np.sum(nominal_unmodeled_acc['covArray'][3:6,3:6]) ]])
 
    # compute RMS absolute error of position and velocity   
    posErrNormNmnlUnmodeledAcc, velErrNormNmnlUnmodeledAcc = \
        errorNorm( trueEphemeris['spacecraft'], est_states )

    # pull out the pickle data of the SRP test unmodeled acceleration run (without SRP)
    SRP_test_unmodeled_acc_data_path = path + '/unmodeled_acc/' + SRP_test_path
    SRP_test_unmodeled_acc_file      = open( SRP_test_unmodeled_acc_data_path, 'rb' )
    SRP_test_unmodeled_acc           = pickle.load( SRP_test_unmodeled_acc_file )
    SRP_test_unmodeled_acc_file.close()

    # unpack data such as truth ephemeris, estimated states, covariance and time span
    trueEphemeris = SRP_test_unmodeled_acc['trueEphemeris']
    est_states = SRP_test_unmodeled_acc['states']
    covar_unmodeledAcc_SRP_test = np.array([[np.sum(SRP_test_unmodeled_acc['covArray'][0:3,0:3]), 0.], [0.,np.sum(SRP_test_unmodeled_acc['covArray'][3:6,3:6]) ]])

    # compute RMS absolute error of position and velocity   
    posErrNormTestUnmodeledAcc, velErrNormTestUnmodeledAcc = \
        errorNorm( trueEphemeris['spacecraft'], est_states )

    ###########################################################
    # 
    # DYNAMICS WORK
    #
    ###########################################################

    # assign the reference propagation of the vanilla nominal to a variable
    refState_nominal  = nominal_vanilla['referenceState']

    # assign the reference propagation of the SRP test nominal to a variable
    refState_SRP_test = SRP_test_vanilla['referenceState']

    # calculate an absolute difference
    refDifference     = refState_nominal.T - refState_SRP_test.T

    ###########################################################
    # 
    # RESIDUAL WORK
    #
    ###########################################################
   
    # initialize arrays for the storage of absolute and relative residuals. Note that the rows
    # correspond to position and velocity
    SRP_test_abs_residual = np.zeros( (2,len(timeSpan)) ) # R[pos,vel] C[length]
    nominal_abs_residual  = np.zeros( (2,len(timeSpan)) ) # R[pos,vel] C[length]
    SRP_test_rel_residual = np.zeros( (2,len(timeSpan)) ) # R[pos,vel] C[length]
    nominal_rel_residual  = np.zeros( (2,len(timeSpan)) ) # R[pos,vel] C[length]

    # calculate the relative residuals between the two vanilla methods and the two unmodeled 
    # acc methods
    SRP_test_rel_residual[0,:] = np.divide(\
                          np.array([posErrNormTestVanilla-posErrNormTestUnmodeledAcc]),\
                          np.array(posErrNormTestVanilla))
    SRP_test_rel_residual[1,:] = np.divide(\
                          np.array([velErrNormTestVanilla-velErrNormTestUnmodeledAcc]),\
                          np.array(velErrNormTestVanilla))

    nominal_rel_residual[0,:] = np.divide(\
                          np.array([posErrNormNmnlVanilla-posErrNormNmnlUnmodeledAcc]),\
                          np.array(posErrNormNmnlVanilla))
    nominal_rel_residual[1,:] = np.divide(\
                          np.array([velErrNormNmnlVanilla-velErrNormNmnlUnmodeledAcc]),\
                          np.array(velErrNormNmnlVanilla))

    # calculate the absolute residuals between the two vanilla methods and the two unmodeled
    # acc methods
    SRP_test_abs_residual[0,:] = np.array([posErrNormTestVanilla-posErrNormTestUnmodeledAcc])
    SRP_test_abs_residual[1,:] = np.array([velErrNormTestVanilla-velErrNormTestUnmodeledAcc])

    nominal_abs_residual[0,:] = np.array([posErrNormNmnlVanilla-posErrNormNmnlUnmodeledAcc])
    nominal_abs_residual[1,:] = np.array([velErrNormNmnlVanilla-velErrNormNmnlUnmodeledAcc])

    ###########################################################
    # 
    # PLOT WORK
    #
    ###########################################################


    # plot the absolute error of the acc method
    plt.subplot(121)
    plt.plot(timeSpan/timeSpan[-1],posErrNormTestUnmodeledAcc, 'b', label='W/O SRP EOM')
    # plt.plot(timeSpan / timeSpan[-1], posErrNormTestUnmodeledAcc + 0.1*np.sqrt(covar_unmodeledAcc_SRP_test[0,0]**2), 'b--')
    # plt.plot(timeSpan / timeSpan[-1], posErrNormTestUnmodeledAcc - 0.1*np.sqrt(covar_unmodeledAcc_SRP_test[0,0]**2), 'b--')
    plt.plot(timeSpan/timeSpan[-1],posErrNormNmnlUnmodeledAcc, 'r', label='With SRP EOM')
    # plt.plot(timeSpan / timeSpan[-1], posErrNormNmnlUnmodeledAcc + 0.1*np.sqrt(covar_unmodeledAcc[0,0]**2), 'r--')
    # plt.plot(timeSpan / timeSpan[-1], posErrNormNmnlUnmodeledAcc - 0.1*np.sqrt(covar_unmodeledAcc[0,0]**2), 'r--')
    # plt.ylabel('N/A')
    #plt.yticks(np.linspace(PosCovarNormMax.max(), PosCovarNormMin.min(), 12))
    plt.xticks([])
    ax = plt.gca()
    #ax.set_xticklabels(['t_0', 't_f'])
    #leg = ax.legend(loc=1,bbox_to_anchor=(.95,.95))
    plt.title('Position Error')
    # plt.ylim((-ymax, ymax))

    plt.subplot(122)
    plt.plot(timeSpan/timeSpan[-1],velErrNormTestUnmodeledAcc, 'b', label='W/O SRP EOM')
    # plt.plot(timeSpan / timeSpan[-1], velErrNormTestUnmodeledAcc + 0.1*np.sqrt(covar_unmodeledAcc_SRP_test[1,1]**2), 'b--')
    # plt.plot(timeSpan / timeSpan[-1], velErrNormTestUnmodeledAcc - 0.1*np.sqrt(covar_unmodeledAcc_SRP_test[1,1]**2), 'b--')
    plt.plot(timeSpan/timeSpan[-1],velErrNormNmnlUnmodeledAcc, 'r', label='With SRP EOM')
    # plt.plot(timeSpan / timeSpan[-1], velErrNormNmnlUnmodeledAcc + 0.1*np.sqrt(covar_unmodeledAcc[1,1]**2), 'r--')
    # plt.plot(timeSpan / timeSpan[-1], velErrNormNmnlUnmodeledAcc - 0.1*np.sqrt(covar_unmodeledAcc[1,1]**2), 'r--')
    #plt.ylabel('km/s')
    #plt.yticks(np.linspace(VelCovarNormMax.max(), VelCovarNormMin.min(), 12))
    plt.xticks([])
    ax = plt.gca()
    #ax.set_xticklabels(['t_0', 't_f'])
    leg = ax.legend(loc=1,bbox_to_anchor=(.99,.99))
    plt.title('Velocity Error')
    # plt.ylim((-ymax, ymax))

    plt.suptitle('Error wrt truth - un-modeled acc filter')
    plt.tight_layout()
    plt.subplots_adjust(top=.9)
    plt.savefig(saveDir + '/errorTruthUnmodeledAcc.png', dpi=300, format='png')
    plt.close()

    ####################################################################################

    # plot the absolute error of the vanilla method
    plt.subplot(121)
    plt.plot(timeSpan/timeSpan[-1],posErrNormTestVanilla, 'b', label='W/O SRP EOM')
    # plt.plot(timeSpan / timeSpan[-1], posErrNormTestVanilla + 0.1*np.sqrt(covar_vanilla_SRP_test[0,0]**2), 'b--')
    # plt.plot(timeSpan / timeSpan[-1], posErrNormTestVanilla - 0.1*np.sqrt(covar_vanilla_SRP_test[0,0]**2), 'b--')
    plt.plot(timeSpan/timeSpan[-1],posErrNormNmnlVanilla, 'r', label='With SRP EOM')
    # plt.plot(timeSpan / timeSpan[-1], posErrNormNmnlVanilla + 0.1*np.sqrt(covar_vanilla[0,0]**2), 'r--')
    # plt.plot(timeSpan / timeSpan[-1], posErrNormNmnlVanilla - 0.1*np.sqrt(covar_vanilla[0,0]**2), 'r--')
    # plt.ylabel('N/A')
    #plt.yticks(np.linspace(PosCovarNormMax.max(), PosCovarNormMin.min(), 12))
    plt.xticks([])
    ax = plt.gca()
    #ax.set_xticklabels(['t_0', 't_f'])
    #leg = ax.legend(loc=1,bbox_to_anchor=(.95,.95))
    plt.title('Position Error')
    # plt.ylim((-ymax, ymax))

    plt.subplot(122)
    plt.plot(timeSpan/timeSpan[-1],velErrNormTestVanilla, 'b', label='W/O SRP EOM')
    # plt.plot(timeSpan / timeSpan[-1], velErrNormTestVanilla + 0.1*np.sqrt(covar_vanilla_SRP_test[1,1]**2), 'b--')
    # plt.plot(timeSpan / timeSpan[-1], velErrNormTestVanilla - 0.1*np.sqrt(covar_vanilla_SRP_test[1,1]**2), 'b--')
    plt.plot(timeSpan/timeSpan[-1],velErrNormNmnlVanilla, 'r', label='With SRP EOM')
    # plt.plot(timeSpan / timeSpan[-1], velErrNormNmnlVanilla + 0.1*np.sqrt(covar_vanilla[1,1]**2), 'r--')
    # plt.plot(timeSpan / timeSpan[-1], velErrNormNmnlVanilla - 0.1*np.sqrt(covar_vanilla[1,1]**2), 'r--')
    #plt.ylabel('km/s')
    #plt.yticks(np.linspace(VelCovarNormMax.max(), VelCovarNormMin.min(), 12))
    plt.xticks([])
    ax = plt.gca()
    #ax.set_xticklabels(['t_0', 't_f'])
    leg = ax.legend(loc=1,bbox_to_anchor=(.99,.99))
    plt.title('Velocity Error')
    # plt.ylim((-ymax, ymax))

    plt.suptitle('Error wrt truth - vanilla filter')
    plt.tight_layout()
    plt.subplots_adjust(top=.9)
    plt.savefig(saveDir + '/errorTruthVanilla.png', dpi=300, format='png')
    plt.close()

    ####################################################################################

    # plot the relative error of the acc method when compared to nominal vanilla 
    plt.subplot(121)
    plt.plot(timeSpan/timeSpan[-1],(SRP_test_rel_residual[0,:]), 'b', label='W/O SRP EOM')
    plt.plot(timeSpan/timeSpan[-1],(nominal_rel_residual[0,:]), 'r', label='With SRP EOM')
    plt.ylabel('N/A')
    #plt.yticks(np.linspace(PosCovarNormMax.max(), PosCovarNormMin.min(), 12))
    plt.xticks([])
    ax = plt.gca()
    #ax.set_xticklabels(['t_0', 't_f'])
    #leg = ax.legend(loc=1,bbox_to_anchor=(.95,.95)) 
    plt.title('Relative Position Residuals')
    # plt.ylim((-ymax, ymax))

    plt.subplot(122)
    plt.plot(timeSpan/timeSpan[-1],(SRP_test_rel_residual[1,:]), 'b', label='W/O SRP EOM')
    plt.plot(timeSpan/timeSpan[-1],(nominal_rel_residual[1,:]), 'r', label='With SRP EOM')
    #plt.ylabel('km/s')
    #plt.yticks(np.linspace(VelCovarNormMax.max(), VelCovarNormMin.min(), 12))
    plt.xticks([])
    ax = plt.gca()
    #ax.set_xticklabels(['t_0', 't_f'])
    leg = ax.legend(loc=1,bbox_to_anchor=(.99,.99)) 
    plt.title('Relative Velocity Residuals')
    # plt.ylim((-ymax, ymax))

    plt.suptitle('Relative Residual Comparison (wrt vanilla)')
    plt.tight_layout()
    plt.subplots_adjust(top=.9)
    plt.savefig(saveDir + '/errorRelativeComparison.png', dpi=300, format='png')
    plt.close()

    ####################################################################################
    # plot the absolute error of the acc method when compared to nominal vanilla 
    plt.subplot(121)
    plt.plot(timeSpan/timeSpan[-1],(SRP_test_abs_residual[0,:]), 'b', label='W/O SRP EOM')
    plt.plot(timeSpan/timeSpan[-1],(nominal_abs_residual[0,:]), 'r', label='With SRP EOM')
    plt.ylabel('km')
    #plt.yticks(np.linspace(PosCovarNormMax.max(), PosCovarNormMin.min(), 12))
    plt.xticks([])
    ax = plt.gca()
    #ax.set_xticklabels(['t_0', 't_f'])
    #leg = ax.legend(loc=1,bbox_to_anchor=(.95,.95)) 
    plt.title('Absolute Position Residuals')
    # plt.ylim((-ymax, ymax))

    plt.subplot(122)
    plt.plot(timeSpan/timeSpan[-1],(SRP_test_abs_residual[1,:]), 'b', label='W/O SRP EOM')
    plt.plot(timeSpan/timeSpan[-1],(nominal_abs_residual[1,:]), 'r', label='With SRP EOM')
    plt.ylabel('km/s')
    #plt.yticks(np.linspace(VelCovarNormMax.max(), VelCovarNormMin.min(), 12))
    plt.xticks([])
    ax = plt.gca()
    #ax.set_xticklabels(['t_0', 't_f'])
    leg = ax.legend(loc=1,bbox_to_anchor=(.99,.99)) 
    plt.title('Absolute Velocity Residuals')
    # plt.ylim((-ymax, ymax))

    plt.suptitle('Absolute Residual Comparison')
    plt.tight_layout()
    plt.subplots_adjust(top=.9)
    plt.savefig(saveDir + '/errorAbsoluteComparison.png', dpi=300, format='png')
    plt.close()

    ####################################################################################
    # plot of the difference in reference trajectories
    plt.subplot(121)
    plt.plot(timeSpan/timeSpan[-1],refDifference[0,:], 'b', label='X')
    plt.plot(timeSpan/timeSpan[-1],refDifference[1,:], 'r', label='Y')
    plt.plot(timeSpan/timeSpan[-1],refDifference[2,:], 'g', label='Z')
    plt.ylabel('km')
    #plt.yticks(np.linspace(PosCovarNormMax.max(), PosCovarNormMin.min(), 12))
    plt.xticks([])
    ax = plt.gca()
    #ax.set_xticklabels(['t_0', 't_f'])
    leg = ax.legend(loc=1,bbox_to_anchor=(.99,.99)) 
    plt.title('Position Difference')
    # plt.ylim((-ymax, ymax))

    plt.subplot(122)
    plt.plot(timeSpan/timeSpan[-1],refDifference[3,:], 'b', label='$\dot{X}$')
    plt.plot(timeSpan/timeSpan[-1],refDifference[4,:], 'r', label='$\dot{Y}$')
    plt.plot(timeSpan/timeSpan[-1],refDifference[5,:], 'g', label='$\dot{Z}$')
    plt.ylabel('km/s')
    #plt.yticks(np.linspace(VelCovarNormMax.max(), VelCovarNormMin.min(), 12))
    plt.xticks([])
    ax = plt.gca()
    #ax.set_xticklabels(['t_0', 't_f'])
    leg = ax.legend(loc=1,bbox_to_anchor=(.99,.99)) 
    plt.title('Velocity Difference')
    # plt.ylim((-ymax, ymax))

    plt.suptitle('Reference Difference (vanilla)')
    plt.tight_layout()
    plt.subplots_adjust(top=.9)
    plt.savefig(saveDir + '/refDifferenceComparison.png', dpi=300, format='png')
    plt.close()

    ####################################################################################
    # plot all the reference trajectories and the truth
    plt.subplot(121)
    plt.plot(timeSpan/timeSpan[-1],nominal_vanilla['states'][:,0], 'b', label='$X_{vanilla}$')
    #plt.plot(timeSpan/timeSpan[-1],nominal_vanilla['states'][:,1], 'r', label='$Y_{vanilla}$')
    #plt.plot(timeSpan/timeSpan[-1],nominal_vanilla['states'][:,2], 'g', label='$Z_{vanilla}$')
    plt.plot(timeSpan/timeSpan[-1],nominal_unmodeled_acc['states'][:,0], 'b--', label='$X_{acc}$')
    plt.plot(timeSpan/timeSpan[-1],trueEphemeris['spacecraft'][0,:], 'r--', label='$X_{truth}$')
    #plt.plot(timeSpan/timeSpan[-1],nominal_unmodeled_acc['states'][:,1], 'r--', label='$Y_{acc}$')
    #plt.plot(timeSpan/timeSpan[-1],nominal_unmodeled_acc['states'][:,2], 'g--', label='$Z_{acc}$')
    plt.ylabel('km')
    #plt.yticks(np.linspace(PosCovarNormMax.max(), PosCovarNormMin.min(), 12))
    plt.xticks([])
    ax = plt.gca()
    #ax.set_xticklabels(['t_0', 't_f'])
    leg = ax.legend(loc=1,bbox_to_anchor=(.99,.99)) 
    plt.title('Position')
    # plt.ylim((-ymax, ymax))

    plt.subplot(122)
    plt.plot(timeSpan/timeSpan[-1],nominal_vanilla['states'][:,3], 'b', label='$\dot{X}_{vanilla}$')
    plt.plot(timeSpan/timeSpan[-1],nominal_unmodeled_acc['states'][:,3], 'b', label='$\dot{X}_{acc}$')
    plt.plot(timeSpan/timeSpan[-1],trueEphemeris['spacecraft'][3,:], 'r--', label='$X_{truth}$')
    #plt.plot(timeSpan/timeSpan[-1],refDifference[4,:], 'r', label='$\dot{Y}$')
    #plt.plot(timeSpan/timeSpan[-1],refDifference[5,:], 'g', label='$\dot{Z}$')
    plt.ylabel('km/s')
    #plt.yticks(np.linspace(VelCovarNormMax.max(), VelCovarNormMin.min(), 12))
    plt.xticks([])
    ax = plt.gca()
    #ax.set_xticklabels(['t_0', 't_f'])
    leg = ax.legend(loc=1,bbox_to_anchor=(.99,.99)) 
    plt.title('Velocity')
    # plt.ylim((-ymax, ymax))

    plt.suptitle('Reference Trajectories vs Truth')
    plt.tight_layout()
    plt.subplots_adjust(top=.9)
    plt.savefig(saveDir + '/refComparison.png', dpi=300, format='png')
    plt.close()

    ############################################################
    #
    # RMS WORK
    #
    ############################################################
    # calculate root mean square values for all errors and compare
    SRP_test_rms   = np.zeros( (2,2) ) # R[pos,vel] C[vanilla,unmodeled]
    nominal_rms    = np.zeros( (2,2) ) # R[pos,vel] C[vanilla,unmodeled]

    SRP_test_rms[0,:]   = np.array([rms(posErrNormTestVanilla),rms(posErrNormTestUnmodeledAcc)])
    SRP_test_rms[1,:]   = np.array([rms(velErrNormTestVanilla),rms(velErrNormTestUnmodeledAcc)])

    nominal_rms[0,:] = np.array([rms(posErrNormNmnlVanilla),rms(posErrNormNmnlUnmodeledAcc)])
    nominal_rms[1,:] = np.array([rms(velErrNormNmnlVanilla),rms(velErrNormNmnlUnmodeledAcc)])

    print SRP_test_rms
    print nominal_rms

    print 'Without SRP EOM, position vanilla/unmodeled'
    print SRP_test_rms[0,0]/SRP_test_rms[0,1]

    print 'With SRP EOM, position vanilla/unmodeled'
    print nominal_rms[0,0]/nominal_rms[0,1]

    print 'Without SRP EOM, velocity vanilla/unmodeled'
    print SRP_test_rms[0,0]/SRP_test_rms[0,1]

    print 'With SRP EOM, velocity vanilla/unmodeled'
    print nominal_rms[1,0]/nominal_rms[1,1]

  return


if __name__ == "__main__":
    main()
