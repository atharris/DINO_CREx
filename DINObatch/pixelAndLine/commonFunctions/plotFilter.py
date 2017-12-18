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
import matplotlib.pyplot as plt
import pylab
import numpy as np

import pdb

################################################################################
#                  E X P O R T E D     F U N C T I O N S:
################################################################################

# -------------------------------------------------------------------------------
def plotFunction(dataDictionary):

    estimatedStates        = dataDictionary['states']
    stateDevHatArray     = dataDictionary['stateDevHatArray']
    covArray               = dataDictionary['covArray']
    prefits                = dataDictionary['prefit residuals']
    postfits               = dataDictionary['postfit residuals']
    beacon_list            = dataDictionary['beacon_list']
    observations           = dataDictionary['Y']
    trueObservations       = dataDictionary['truth']
    timeSpan               = dataDictionary['timeSpan']
    dirIt                  = dataDictionary['dirIt']
    err                    = dataDictionary['err']
    stateErrorHat          = dataDictionary['stateErrorHat']
    observationUncertainty = dataDictionary['obs_uncertainty']
    referenceState         = dataDictionary['referenceState']
    trueEphemeris          = dataDictionary['trueEphemeris']
    extras                 = dataDictionary['extras']
    postDelta              = dataDictionary['postfit changes']

    stand_devs = 3. * np.array([np.sqrt(np.fabs(np.diag(P))) for P in covArray])

    plt.figure(1)
    plt.subplot(321)
    plt.plot(timeSpan, stand_devs[:, 0], 'r--', timeSpan, -1 * stand_devs[:, 0], 'r--')
    plt.plot(timeSpan, stateErrorHat[:, 0], 'k.')
    plt.ylabel('x (km)')
    plt.xticks([])
    plt.title('Position')
    plt.ylim((-1.1*np.max(stand_devs[:, 0]), 1.1*np.max(stand_devs[:, 0])))

    plt.subplot(323)
    plt.plot(timeSpan, stand_devs[:, 1], 'r--', timeSpan, -1 * stand_devs[:, 1], 'r--')
    plt.plot(timeSpan, stateErrorHat[:, 1], 'k.')
    plt.ylabel('y (km)')
    plt.xticks([])
    plt.ylim((-1.1*np.max(stand_devs[:, 1]), 1.1*np.max(stand_devs[:, 1])))

    plt.subplot(325)
    plt.plot(timeSpan, stand_devs[:, 2], 'r--', timeSpan, -1 * stand_devs[:, 2], 'r--')
    plt.plot(timeSpan, stateErrorHat[:, 2], 'k.')
    plt.ylabel('z (km)')
    plt.xticks([min(timeSpan), max(timeSpan)], rotation=30, ha='right')
    ax = plt.gca()
    ax.set_xticklabels(['t_0', 't_f'])
    plt.ylim((-1.1*np.max(stand_devs[:, 2]), 1.1*np.max(stand_devs[:, 2])))

    plt.subplot(322)
    plt.plot(timeSpan, stand_devs[:, 3], 'r--', timeSpan, -1 * stand_devs[:, 3], 'r--')
    plt.plot(timeSpan, stateErrorHat[:, 3], 'k.')
    plt.ylabel('vx (km/s)')
    plt.xticks([])
    plt.title('Velocity')
    plt.ylim((-1.1*np.max(stand_devs[:, 3]), 1.1*np.max(stand_devs[:, 3])))

    plt.subplot(324)
    plt.plot(timeSpan, stand_devs[:, 4], 'r--', timeSpan, -1 * stand_devs[:, 4], 'r--')
    plt.plot(timeSpan, stateErrorHat[:, 4], 'k.')
    plt.ylabel('vy (km/s)')
    plt.xticks([])
    plt.ylim((-1.1*np.max(stand_devs[:, 4]), 1.1*np.max(stand_devs[:, 4])))

    plt.subplot(326)
    plt.plot(timeSpan, stand_devs[:, 5], 'r--', timeSpan, -1 * stand_devs[:, 5], 'r--')
    plt.plot(timeSpan, stateErrorHat[:, 5], 'k.')
    plt.ylabel('vz (km/s)')
    plt.xticks([min(timeSpan), max(timeSpan)], rotation=30, ha='right')
    ax = plt.gca()
    ax.set_xticklabels(['t_0', 't_f'])
    plt.ylim((-1.1*np.max(stand_devs[:, 5]), 1.1*np.max(stand_devs[:, 5])))

    plt.suptitle('Error Between Estimated State' + ' $\hat{X}$ ' + ' and Truth with Covariance Bounds')
    plt.tight_layout()
    plt.subplots_adjust(top=.9)
    plt.savefig(dirIt+'/State Errors and Covariance Bounds.png', dpi=300, format='png')


    if dataDictionary['acc_est'] == 'OFF':
      plt.figure(111)
      plt.subplot(321)
      plt.plot(timeSpan, estimatedStates[:, 0], 'b.')
      plt.ylabel('x (km)')
      plt.xticks([])
      plt.title('Position')

      plt.subplot(323)
      plt.plot(timeSpan, estimatedStates[:, 1], 'b.')
      plt.ylabel('y (km)')
      plt.xticks([])

      plt.subplot(325)
      plt.plot(timeSpan, estimatedStates[:, 2], 'b.')
      plt.ylabel('z (km)')
      plt.xticks([min(timeSpan), max(timeSpan)], rotation=30, ha='right')
      ax = plt.gca()
      ax.set_xticklabels(['t_0', 't_f'])

      plt.subplot(322)
      plt.plot(timeSpan, estimatedStates[:, 3], 'b.')
      plt.ylabel('vx (km/s)')
      plt.xticks([])
      plt.title('Velocity')


      plt.subplot(324)
      plt.plot(timeSpan, estimatedStates[:, 4], 'b.')
      plt.ylabel('vy (km/s)')
      plt.xticks([])

      plt.subplot(326)
      plt.plot(timeSpan, estimatedStates[:, 5], 'b.')
      plt.ylabel('vz (km/s)')
      plt.xticks([min(timeSpan), max(timeSpan)], rotation=30, ha='right')
      ax = plt.gca()
      ax.set_xticklabels(['t_0', 't_f'])

      plt.suptitle('States' + ' $X$' )
      plt.tight_layout()
      plt.subplots_adjust(top=.9)
      plt.savefig(dirIt+'/States.png', dpi=300, format='png')
    elif dataDictionary['acc_est'] == 'ON':
      plt.figure(111)
      plt.subplot(331)
      plt.plot(timeSpan, estimatedStates[:, 0], 'b.')
      plt.ylabel('x (km)')
      plt.xticks([])
      plt.title('Position')

      plt.subplot(334)
      plt.plot(timeSpan, estimatedStates[:, 1], 'b.')
      plt.ylabel('y (km)')
      plt.xticks([])

      plt.subplot(337)
      plt.plot(timeSpan, estimatedStates[:, 2], 'b.')
      plt.ylabel('z (km)')
      plt.xticks([min(timeSpan), max(timeSpan)], rotation=30, ha='right')
      ax = plt.gca()
      ax.set_xticklabels(['t_0', 't_f'])

      plt.subplot(332)
      plt.plot(timeSpan, estimatedStates[:, 3], 'b.')
      plt.ylabel('vx (km/s)')
      plt.xticks([])
      plt.title('Velocity')

      plt.subplot(335)
      plt.plot(timeSpan, estimatedStates[:, 4], 'b.')
      plt.ylabel('vy (km/s)')
      plt.xticks([])

      plt.subplot(338)
      plt.plot(timeSpan, estimatedStates[:, 5], 'b.')
      plt.ylabel('vz (km/s)')
      plt.xticks([min(timeSpan), max(timeSpan)], rotation=30, ha='right')
      ax = plt.gca()
      ax.set_xticklabels(['t_0', 't_f'])

      plt.subplot(333)
      plt.plot(timeSpan, estimatedStates[:, 6], 'b.')
      plt.ylabel('ax (km/s)')
      plt.ylim(-max(estimatedStates[:,6])*3,max(estimatedStates[:, 6])*3)
      plt.xticks([])
      plt.title('Accelerations')


      plt.subplot(336)
      plt.plot(timeSpan, estimatedStates[:, 7], 'b.')
      plt.ylabel('ay (km/s)')
      plt.ylim(-max(estimatedStates[:,7])*3,max(estimatedStates[:, 7])*3)
      plt.xticks([])

      plt.subplot(339)
      plt.plot(timeSpan, estimatedStates[:, 8], 'b.')
      plt.ylabel('az (km/s)')
      plt.xticks([min(timeSpan), max(timeSpan)], rotation=30, ha='right')
      plt.ylim(-max(estimatedStates[:,8])*3,max(estimatedStates[:, 8])*3)
      ax = plt.gca()
      ax.set_xticklabels(['t_0', 't_f'])

      plt.suptitle('States' + ' $X$' )
      plt.tight_layout()
      plt.subplots_adjust(top=.9)
      plt.savefig(dirIt+'/States.png', dpi=300, format='png')
    else:
      print 'Did not specify if plotting estimated accelerations was on or off'

    plt.figure(2)
    plt.subplot(321)
    plt.plot(timeSpan, stateDevHatArray[:, 0])
    plt.ylabel('x (km)')
    plt.xticks([])
    plt.title('Position')
    # plt.ylim((-ymax, ymax))

    plt.subplot(323)
    plt.plot(timeSpan, stateDevHatArray[:, 1])
    plt.ylabel('y (km)')
    plt.xticks([])
    # plt.ylim((-ymax, ymax))

    plt.subplot(325)
    plt.plot(timeSpan, stateDevHatArray[:, 2])
    plt.ylabel('z (km)')
    plt.xticks([min(timeSpan), max(timeSpan)], rotation=30, ha='right')
    ax = plt.gca()
    ax.set_xticklabels(['t_0', 't_f'])
    # plt.ylim((-ymax, ymax))

    plt.subplot(322)
    plt.plot(timeSpan, stateDevHatArray[:, 3])
    plt.ylabel('vx (km/s)')
    plt.xticks([])
    plt.title('Velocity')
    # plt.ylim((-ymax, ymax))

    plt.subplot(324)
    plt.plot(timeSpan, stateDevHatArray[:, 4])
    plt.ylabel('vy (km/s)')
    plt.xticks([])
    # plt.ylim((-ymax, ymax))

    plt.subplot(326)
    plt.plot(timeSpan, stateDevHatArray[:, 5])
    plt.ylabel('vz (km/s)')
    plt.xticks([min(timeSpan), max(timeSpan)], rotation=30, ha='right')
    ax = plt.gca()
    ax.set_xticklabels(['t_0', 't_f'])
    # plt.ylim((-ymax, ymax))
    plt.suptitle('State Deviation Estimates ' + '$\hat{x}$')
    plt.tight_layout()
    plt.subplots_adjust(top=.9)
    plt.savefig(dirIt + '/State Deviation Estimates.png', dpi=300, format='png')

    plt.figure(3)
    plt.subplot(321)
    plt.plot(timeSpan, err[:, 0])
    plt.xticks([])
    plt.ylabel('x (km)')
    plt.title('Position')

    plt.subplot(323)
    plt.plot(timeSpan, err[:, 1])
    plt.xticks([])
    plt.ylabel('y (km)')

    plt.subplot(325)
    plt.plot(timeSpan, err[:, 2])
    plt.ylabel('z (km)')
    plt.xticks([min(timeSpan), max(timeSpan)], rotation=30, ha='right')
    ax = plt.gca()
    ax.set_xticklabels(['t_0', 't_f'])

    plt.subplot(322)
    plt.plot(timeSpan, err[:, 3])
    plt.ylabel('vx (km/s)')
    plt.xticks([])
    plt.title('Velocity')

    plt.subplot(324)
    plt.plot(timeSpan, err[:, 4])
    plt.ylabel('vy (km/s)')
    plt.xticks([])
    plt.subplot(326)
    plt.plot(timeSpan, err[:, 5])
    plt.ylabel('vz (km/s)')
    plt.xticks([min(timeSpan), max(timeSpan)], rotation=30, ha='right')
    ax = plt.gca()
    ax.set_xticklabels(['t_0', 't_f'])
    plt.suptitle('Reference Trajectory Error vs Truth Trajectory: Dynamics Errors')
    plt.tight_layout()
    plt.subplots_adjust(top=.9)
    plt.savefig(dirIt + '/Reference Trajectory Error.png', dpi=300, format='png')

    plt.figure(4)
    plt.subplot(121)
    plt.plot(timeSpan, prefits[:, 0], '.')
    plt.plot(timeSpan, 3*np.ones(len(timeSpan))*observationUncertainty[0,0]**2, 'r--')
    plt.plot(timeSpan, -3*np.ones(len(timeSpan))*observationUncertainty[0,0]**2, 'r--')
    plt.xticks([])
    plt.title('Pixel (-)')
    plt.ylabel('Residual')
    plt.xticks(timeSpan, rotation=90, ha='right')
    ax = plt.gca()
    ax.set_xticklabels(beacon_list)

    plt.subplot(122)
    plt.plot(timeSpan, prefits[:, 1], '.')
    plt.plot(timeSpan, 3*np.ones(len(timeSpan))*observationUncertainty[1, 1], 'r--')
    plt.plot(timeSpan, -3*np.ones(len(timeSpan))*observationUncertainty[1, 1], 'r--')
    plt.xticks([])
    plt.title('Line (-)')
    plt.xticks(timeSpan, rotation=90, ha='right')
    ax = plt.gca()
    ax.set_xticklabels(beacon_list)

    plt.suptitle('Prefit Observation Residuals')
    plt.tight_layout()
    plt.subplots_adjust(top=.9)
    plt.savefig(dirIt + '/Prefit Observation Residuals.png', dpi=300, format='png')

    plt.figure(5)
    plt.subplot(121)
    plt.plot(timeSpan, postfits[:, 0], '.')
    plt.plot(timeSpan, 3*np.ones(len(timeSpan))*observationUncertainty[0,0], 'r--')
    plt.plot(timeSpan, -3*np.ones(len(timeSpan))*observationUncertainty[0,0], 'r--')
    plt.xticks([])
    plt.ylabel('Residual')
    plt.title('Pixel (-)')
    plt.xticks(timeSpan, rotation=90, ha='right')
    ax = plt.gca()
    ax.set_xticklabels(beacon_list)

    plt.subplot(122)
    plt.plot(timeSpan, postfits[:, 1], '.')
    plt.plot(timeSpan, 3*np.ones(len(timeSpan))*observationUncertainty[1, 1], 'r--')
    plt.plot(timeSpan, -3*np.ones(len(timeSpan))*observationUncertainty[1, 1], 'r--')
    plt.xticks([])
    plt.title('Line (-)')
    plt.xticks(timeSpan, rotation=90, ha='right')
    ax = plt.gca()
    ax.set_xticklabels(beacon_list)


    plt.suptitle('Postfit Observation Residuals')
    plt.tight_layout()
    plt.subplots_adjust(top=.9)
    plt.savefig(dirIt + '/Postfit Observation Residuals.png', dpi=300, format='png')

    plt.figure(55)
    plt.subplot(121)
    plt.plot(timeSpan, postDelta[:, 0], '.')
    plt.xticks([])
    plt.ylabel('Residual')
    plt.title('Pixel (-)')
    plt.xticks(timeSpan, rotation=90, ha='right')
    ax = plt.gca()
    ax.set_xticklabels(beacon_list)

    plt.subplot(122)
    plt.plot(timeSpan, postDelta[:, 1], '.')
    plt.xticks([])
    plt.title('Line (-)')
    plt.xticks(timeSpan, rotation=90, ha='right')
    ax = plt.gca()
    ax.set_xticklabels(beacon_list)


    plt.suptitle('Postfit Changes')
    plt.tight_layout()
    plt.subplots_adjust(top=.9)
    plt.savefig(dirIt + '/Postfit Delta.png', dpi=300, format='png')

    plt.figure(77)

    PosErrorNorm = np.zeros(len(estimatedStates[:,0]))
    VelErrorNorm = np.zeros(len(estimatedStates[:, 0]))
    PosCovarNormMax =  np.zeros(len(estimatedStates[:, 0]))
    PosCovarNormMin =  np.zeros(len(estimatedStates[:, 0]))
    VelCovarNormMax =  np.zeros(len(estimatedStates[:, 0]))
    VelCovarNormMin =  np.zeros(len(estimatedStates[:, 0]))
    for i in range(len(estimatedStates[:,0])):
        PosErrorNorm[i] = np.linalg.norm(trueEphemeris['spacecraft'].T[i, 0:3] - estimatedStates[i, 0:3])
        VelErrorNorm[i] = np.linalg.norm(trueEphemeris['spacecraft'].T[i, 3:6] - estimatedStates[i, 3:6])
        PosCovarNormMax[i] = PosErrorNorm[i] + np.sqrt(stand_devs[i,0]**2 + stand_devs[i,1]**2 + stand_devs[i,2]**2)/30.
        PosCovarNormMin[i] = PosErrorNorm[i] - np.sqrt(stand_devs[i,0]**2 + stand_devs[i,1]**2 + stand_devs[i,2]**2)/30.
        VelCovarNormMax[i] = VelErrorNorm[i] + np.sqrt(stand_devs[i,3]**2 + stand_devs[i,4]**2 + stand_devs[i,5]**2)/30.
        VelCovarNormMin[i] = VelErrorNorm[i] - np.sqrt(stand_devs[i,3]**2 + stand_devs[i,4]**2 + stand_devs[i,5]**2)/30.

    plt.subplot(121)
    plt.plot(timeSpan,PosErrorNorm, 'b')
    plt.plot(timeSpan,PosCovarNormMax, 'r--', label= '$\sigma$/10')
    plt.plot(timeSpan,PosCovarNormMin, 'r--')
    plt.ylabel('$\delta$ R (km)')
    plt.yticks(np.linspace(PosCovarNormMax.max(), PosCovarNormMin.min(), 12))
    plt.xticks([])
    ax = plt.gca()
    ax.set_xticklabels(['t_0', 't_f'])
    plt.title('Position Error')
    # plt.ylim((-ymax, ymax))

    plt.subplot(122)
    plt.plot(timeSpan, VelErrorNorm, 'b')
    plt.plot(timeSpan,VelCovarNormMax, 'r--',label= '$\sigma$/10')
    plt.plot(timeSpan,VelCovarNormMin, 'r--')
    plt.ylabel('$\delta$ V (km/s)')
    plt.yticks(np.linspace(VelCovarNormMax.max(), VelCovarNormMin.min(), 12))
    plt.xticks([])
    ax = plt.gca()
    ax.set_xticklabels(['t_0', 't_f'])
    plt.title('Velocity Error')
    # plt.ylim((-ymax, ymax))

    plt.suptitle('Estimated state Error Norm from Truth data')
    plt.tight_layout()
    plt.subplots_adjust(top=.9)
    plt.savefig(dirIt + '/Truth-EstimatedStateNorms.png', dpi=300, format='png')

    plt.figure(6)
    plt.scatter(referenceState[:, 1], referenceState[:, 2], color='r')
    plt.scatter(trueEphemeris['spacecraft'].T[:, 1], trueEphemeris['spacecraft'].T[:, 2], color='b')
    plt.suptitle('Reference Trajectory (r) vs Truth Trajectory (b)')
    plt.savefig(dirIt + '/Trajectories.png', dpi=300, format='png')

    jet = plt.get_cmap('jet')
    colors = iter(jet(np.linspace(0, 1, extras['n_unique_beacons'])))

    jet = plt.get_cmap('jet')
    colors = iter(jet(np.linspace(0, 1, extras['n_unique_beacons'])))

    plt.figure(8)
    plt.subplot(321)
    plt.plot(timeSpan, stand_devs[:, 0], 'r-', timeSpan, covArray[:, 0, 0], 'b-')
    plt.ylabel('x (km)')
    plt.xticks([])
    ax = plt.gca()
    ax.set_yscale('log')
    plt.title('Position')

    plt.subplot(323)
    plt.plot(timeSpan, stand_devs[:, 1], 'r-', timeSpan, covArray[:, 1, 1], 'b-')
    plt.ylabel('y (km)')
    plt.xticks([])
    ax = plt.gca()
    ax.set_yscale('log')

    plt.subplot(325)
    plt.plot(timeSpan, stand_devs[:, 2], 'r-', timeSpan, covArray[:, 2, 2], 'b-')
    plt.ylabel('z (km)')
    plt.xticks([min(timeSpan), max(timeSpan)], rotation=30, ha='right')
    ax = plt.gca()
    ax.set_yscale('log')
    ax.set_xticklabels(['t_0', 't_f'])

    plt.subplot(322)
    plt.plot(timeSpan, stand_devs[:, 3], 'r-', timeSpan, covArray[:, 3, 3], 'b-')
    plt.ylabel('vx (km/s)')
    plt.xticks([])
    ax = plt.gca()
    ax.set_yscale('log')
    plt.title('Velocity')

    plt.subplot(324)
    plt.plot(timeSpan, stand_devs[:, 4], 'r-', timeSpan, covArray[:, 4, 4], 'b-')
    plt.ylabel('vy (km/s)')
    ax = plt.gca()
    ax.set_yscale('log')
    plt.xticks([])

    plt.subplot(326)
    plt.plot(timeSpan, stand_devs[:, 5], 'r-', timeSpan, covArray[:, 5, 5], 'b-')
    plt.ylabel('vz (km/s)')
    plt.xticks([min(timeSpan), max(timeSpan)], rotation=30, ha='right')
    ax = plt.gca()
    ax.set_yscale('log')
    ax.set_xticklabels(['t_0', 't_f'])
    plt.suptitle('Standard Deviation and Covariance')
    plt.tight_layout()
    plt.subplots_adjust(top=.9)
    plt.savefig(dirIt + '/Standard Deviation and Covariance.png', dpi=300, format='png')

    plt.figure(9)
    plt.subplot(121)
    plt.plot(timeSpan, observations[:, 0] - trueObservations[:, 0])
    plt.xticks([])
    plt.title('Range (km)')
    plt.ylabel('Added Error')
    plt.xticks(timeSpan, rotation=90, ha='right')
    ax = plt.gca()
    ax.set_xticklabels(beacon_list)

    plt.subplot(122)
    plt.plot(timeSpan, observations[:, 1] - trueObservations[:, 1])
    plt.xticks([])
    plt.title('Range Rate (km/s)')
    plt.xticks(timeSpan, rotation=90, ha='right')
    ax = plt.gca()
    ax.set_xticklabels(beacon_list)

    plt.suptitle('Observation Errors: Measured vs Truth')
    plt.tight_layout()
    plt.subplots_adjust(top=.9)
    plt.savefig(dirIt + '/Observation Errors.png', dpi=300, format='png')

    plt.figure(10)
    plt.subplot(321)
    plt.plot(timeSpan, stateErrorHat[:, 0])
    plt.ylabel('x (km)')
    plt.xticks([])
    plt.title('Position')
    # plt.ylim((-ymax, ymax))

    plt.subplot(323)
    plt.plot(timeSpan, stateErrorHat[:, 1])
    plt.ylabel('y (km)')
    plt.xticks([])
    # plt.ylim((-ymax, ymax))

    plt.subplot(325)
    plt.plot(timeSpan, stateErrorHat[:, 2])
    plt.ylabel('z (km)')
    plt.xticks([min(timeSpan), max(timeSpan)], rotation=30, ha='right')
    ax = plt.gca()
    ax.set_xticklabels(['t_0', 't_f'])
    # plt.ylim((-ymax, ymax))

    plt.subplot(322)
    plt.plot(timeSpan, stateErrorHat[:, 3])
    plt.ylabel('vx (km/s)')
    plt.xticks([])
    plt.title('Velocity')
    # plt.ylim((-ymax, ymax))

    plt.subplot(324)
    plt.plot(timeSpan, stateErrorHat[:, 4])
    plt.ylabel('vy (km/s)')
    plt.xticks([])
    # plt.ylim((-ymax, ymax))

    plt.subplot(326)
    plt.plot(timeSpan, stateErrorHat[:, 5])
    plt.ylabel('vz (km/s)')
    plt.xticks([min(timeSpan), max(timeSpan)], rotation=30, ha='right')
    ax = plt.gca()
    ax.set_xticklabels(['t_0', 't_f'])
    # plt.ylim((-ymax, ymax))
    plt.suptitle('Estimated State Error WRT Truth')
    plt.tight_layout()
    plt.subplots_adjust(top=.9)
    plt.savefig(dirIt + '/Estimated State Error.png', dpi=300, format='png')

    plt.close('all')
    return

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

def main( input ):

  return


if __name__ == "__main__":
    main()
