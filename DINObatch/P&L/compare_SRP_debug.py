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
    nominal_path = 'Batch_Iteration3/nominal_data.pkl'
    debug_path   = 'Batch_Iteration3/debug_data.pkl'

    nominal_vanilla_data_path = path + '/vanilla_pl/' + nominal_path
    nominal_vanilla_file      = open( nominal_vanilla_data_path, 'rb' )
    nominal_vanilla           = pickle.load( nominal_vanilla_file )
    nominal_vanilla_file.close()

    true_ephem = nominal_vanilla['true_ephem']
    est_states = nominal_vanilla['states']
    t_span     = nominal_vanilla['t_span']
    
    posErrNormNmnlVanilla, velErrNormNmnlVanilla = \
        errorNorm( true_ephem['spacecraft'], est_states )

    debug_vanilla_data_path = path + '/vanilla_pl/' + debug_path
    debug_vanilla_file      = open( debug_vanilla_data_path, 'rb' )
    debug_vanilla           = pickle.load( debug_vanilla_file )
    debug_vanilla_file.close()

    true_ephem = debug_vanilla['true_ephem']
    est_states = debug_vanilla['states']
    
    posErrNormDbgVanilla, velErrNormDbgVanilla = \
        errorNorm( true_ephem['spacecraft'], est_states )

    nominal_unmodeled_acc_data_path = path + '/unmodeled_acc/' + nominal_path
    nominal_unmodeled_acc_file      = open( nominal_unmodeled_acc_data_path, 'rb' )
    nominal_unmodeled_acc           = pickle.load( nominal_unmodeled_acc_file )
    nominal_unmodeled_acc_file.close()

    true_ephem = nominal_unmodeled_acc['true_ephem']
    est_states = nominal_unmodeled_acc['states']
    
    posErrNormNmnlUnmodeledAcc, velErrNormNmnlUnmodeledAcc = \
        errorNorm( true_ephem['spacecraft'], est_states )

    debug_unmodeled_acc_data_path = path + '/unmodeled_acc/' + debug_path
    debug_unmodeled_acc_file      = open( debug_unmodeled_acc_data_path, 'rb' )
    debug_unmodeled_acc           = pickle.load( debug_unmodeled_acc_file )
    debug_unmodeled_acc_file.close()

    true_ephem = debug_unmodeled_acc['true_ephem']
    est_states = debug_unmodeled_acc['states']
    
    posErrNormDbgUnmodeledAcc, velErrNormDbgUnmodeledAcc = \
        errorNorm( true_ephem['spacecraft'], est_states )

    ###########################################################
    # 
    # RESIDUAL WORK
    #
    ###########################################################

    debug_abs_residual   = np.zeros( (2,len(t_span)) ) # R[pos,vel] C[length]
    nominal_abs_residual = np.zeros( (2,len(t_span)) ) # R[pos,vel] C[length]
    debug_rel_residual   = np.zeros( (2,len(t_span)) ) # R[pos,vel] C[length]
    nominal_rel_residual = np.zeros( (2,len(t_span)) ) # R[pos,vel] C[length]
    pdb.set_trace()
    debug_rel_residual[0,:] = np.divide(\
                          np.array([posErrNormDbgVanilla-posErrNormDbgUnmodeledAcc]),\
                          np.array(posErrNormDbgVanilla))
    debug_rel_residual[1,:] = np.divide(\
                          np.array([velErrNormDbgVanilla-velErrNormDbgUnmodeledAcc]),\
                          np.array(velErrNormDbgVanilla))

    nominal_rel_residual[0,:] = np.divide(\
                          np.array([posErrNormNmnlVanilla-posErrNormNmnlUnmodeledAcc]),\
                          np.array(posErrNormNmnlVanilla))
    nominal_rel_residual[1,:] = np.divide(\
                          np.array([velErrNormNmnlVanilla-velErrNormNmnlUnmodeledAcc]),\
                          np.array(velErrNormNmnlVanilla))

    debug_abs_residual[0,:] = np.array([posErrNormDbgVanilla-posErrNormDbgUnmodeledAcc])
    debug_abs_residual[1,:] = np.array([velErrNormDbgVanilla-velErrNormDbgUnmodeledAcc])

    nominal_abs_residual[0,:] = np.array([posErrNormNmnlVanilla-posErrNormNmnlUnmodeledAcc])
    nominal_abs_residual[1,:] = np.array([velErrNormNmnlVanilla-velErrNormNmnlUnmodeledAcc])

    plt.subplot(121)
    plt.plot(t_span/t_span[-1],np.abs(debug_rel_residual[0,:]), 'b', label='W/O SRP EOM')
    plt.plot(t_span/t_span[-1],np.abs(nominal_rel_residual[0,:]), 'r', label='With SRP EOM')
    plt.ylabel('N/A')
    #plt.yticks(np.linspace(PosCovarNormMax.max(), PosCovarNormMin.min(), 12))
    plt.xticks([])
    ax = plt.gca()
    #ax.set_xticklabels(['t_0', 't_f'])
    #leg = ax.legend(loc=1,bbox_to_anchor=(.95,.95)) 
    plt.title('Relative Position Residuals')
    # plt.ylim((-ymax, ymax))

    plt.subplot(122)
    plt.plot(t_span/t_span[-1],np.abs(debug_rel_residual[1,:]), 'b', label='W/O SRP EOM')
    plt.plot(t_span/t_span[-1],np.abs(nominal_rel_residual[1,:]), 'r', label='With SRP EOM')
    #plt.ylabel('km/s')
    #plt.yticks(np.linspace(VelCovarNormMax.max(), VelCovarNormMin.min(), 12))
    plt.xticks([])
    ax = plt.gca()
    #ax.set_xticklabels(['t_0', 't_f'])
    leg = ax.legend(loc=1,bbox_to_anchor=(.99,.99)) 
    plt.title('Relative Velocity Residuals')
    # plt.ylim((-ymax, ymax))

    plt.suptitle('Residual Comparison')
    plt.tight_layout()
    plt.subplots_adjust(top=.9)
    plt.savefig(path + '/errorRelativeComparison.png', dpi=300, format='png')
    plt.close()

    plt.subplot(121)
    plt.plot(t_span/t_span[-1],np.abs(debug_abs_residual[0,:]), 'b', label='W/O SRP EOM')
    plt.plot(t_span/t_span[-1],np.abs(nominal_abs_residual[0,:]), 'r', label='With SRP EOM')
    plt.ylabel('km')
    #plt.yticks(np.linspace(PosCovarNormMax.max(), PosCovarNormMin.min(), 12))
    plt.xticks([])
    ax = plt.gca()
    #ax.set_xticklabels(['t_0', 't_f'])
    #leg = ax.legend(loc=1,bbox_to_anchor=(.95,.95)) 
    plt.title('Relative Position Residuals')
    # plt.ylim((-ymax, ymax))

    plt.subplot(122)
    plt.plot(t_span/t_span[-1],np.abs(debug_abs_residual[1,:]), 'b', label='W/O SRP EOM')
    plt.plot(t_span/t_span[-1],np.abs(nominal_abs_residual[1,:]), 'r', label='With SRP EOM')
    plt.ylabel('km/s')
    #plt.yticks(np.linspace(VelCovarNormMax.max(), VelCovarNormMin.min(), 12))
    plt.xticks([])
    ax = plt.gca()
    #ax.set_xticklabels(['t_0', 't_f'])
    leg = ax.legend(loc=1,bbox_to_anchor=(.99,.99)) 
    plt.title('Relative Velocity Residuals')
    # plt.ylim((-ymax, ymax))

    plt.suptitle('Residual Comparison')
    plt.tight_layout()
    plt.subplots_adjust(top=.9)
    plt.savefig(path + '/errorAbsoluteComparison.png', dpi=300, format='png')

    ############################################################
    #
    # RMS WORK
    #
    ############################################################

    debug_rms   = np.zeros( (2,2) ) # R[pos,vel] C[vanilla,unmodeled]
    nominal_rms = np.zeros( (2,2) ) # R[pos,vel] C[vanilla,unmodeled]

    debug_rms[0,:]   = np.array([rms(posErrNormDbgVanilla),rms(posErrNormDbgUnmodeledAcc)])
    debug_rms[1,:]   = np.array([rms(velErrNormDbgVanilla),rms(velErrNormDbgUnmodeledAcc)])

    nominal_rms[0,:] = np.array([rms(posErrNormNmnlVanilla),rms(posErrNormNmnlUnmodeledAcc)])
    nominal_rms[1,:] = np.array([rms(velErrNormNmnlVanilla),rms(velErrNormNmnlUnmodeledAcc)])

    print debug_rms
    print nominal_rms

    print 'Without SRP EOM, position vanilla/unmodeled'
    print debug_rms[0,0]/debug_rms[0,1]

    print 'With SRP EOM, position vanilla/unmodeled'
    print nominal_rms[0,0]/nominal_rms[0,1]

    print 'Without SRP EOM, velocity vanilla/unmodeled'
    print debug_rms[0,0]/debug_rms[0,1]

    print 'With SRP EOM, velocity vanilla/unmodeled'
    print nominal_rms[1,0]/nominal_rms[1,1]

    return


if __name__ == "__main__":
    main()
