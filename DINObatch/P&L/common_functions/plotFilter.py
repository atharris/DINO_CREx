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


################################################################################
#                  E X P O R T E D     F U N C T I O N S:
################################################################################

# -------------------------------------------------------------------------------
def plotFunction(data_dict):

    est_states = data_dict['states']
    x_hat_array = data_dict['x_hat_array']
    P_array = data_dict['P_array']
    prefits = data_dict['prefit residuals']
    postfits = data_dict['postfit residuals']
    beacon_list = data_dict['beacon_list']
    Y_observed = data_dict['Y']
    Y_truth = data_dict['truth']
    t_span  = data_dict['t_span']
    dirIt   = data_dict['dirIt']
    err     = data_dict['err']
    err_hat = data_dict['err_hat']
    observation_uncertainty = data_dict['obs_uncertainty']
    ref_state = data_dict['ref_state']
    true_ephem = data_dict['true_ephem']
    extras = data_dict['extras']
    postDelta = data_dict['postfit changes']

    stand_devs = 3. * np.array([np.sqrt(np.fabs(np.diag(P))) for P in P_array])

    if data_dict['acc_est'] == 'OFF':
      plt.figure(1)
      plt.subplot(321)
      plt.plot(t_span, stand_devs[:, 0], 'r--', t_span, -1 * stand_devs[:, 0], 'r--')
      plt.plot(t_span, x_hat_array[:, 0], 'k.')
      plt.ylabel('x (km)')
      plt.xticks([])
      plt.title('Position')
      plt.ylim((-1.1*np.max(stand_devs[:, 0]), 1.1*np.max(stand_devs[:, 0])))

      plt.subplot(323)
      plt.plot(t_span, stand_devs[:, 1], 'r--', t_span, -1 * stand_devs[:, 1], 'r--')
      plt.plot(t_span, x_hat_array[:, 1], 'k.')
      plt.ylabel('y (km)')
      plt.xticks([])
      plt.ylim((-1.1*np.max(stand_devs[:, 1]), 1.1*np.max(stand_devs[:, 1])))

      plt.subplot(325)
      plt.plot(t_span, stand_devs[:, 2], 'r--', t_span, -1 * stand_devs[:, 2], 'r--')
      plt.plot(t_span, x_hat_array[:, 2], 'k.')
      plt.ylabel('z (km)')
      plt.xticks([min(t_span), max(t_span)], rotation=30, ha='right')
      ax = plt.gca()
      ax.set_xticklabels(['t_0', 't_f'])
      plt.ylim((-1.1*np.max(stand_devs[:, 2]), 1.1*np.max(stand_devs[:, 2])))

      plt.subplot(322)
      plt.plot(t_span, stand_devs[:, 3], 'r--', t_span, -1 * stand_devs[:, 3], 'r--')
      plt.plot(t_span, x_hat_array[:, 3], 'k.')
      plt.ylabel('vx (km/s)')
      plt.xticks([])
      plt.title('Velocity')
      plt.ylim((-1.1*np.max(stand_devs[:, 3]), 1.1*np.max(stand_devs[:, 3])))


      plt.subplot(324)
      plt.plot(t_span, stand_devs[:, 4], 'r--', t_span, -1 * stand_devs[:, 4], 'r--')
      plt.plot(t_span, x_hat_array[:, 4], 'k.')
      plt.ylabel('vy (km/s)')
      plt.xticks([])
      plt.ylim((-1.1*np.max(stand_devs[:, 4]), 1.1*np.max(stand_devs[:, 4])))

      plt.subplot(326)
      plt.plot(t_span, stand_devs[:, 5], 'r--', t_span, -1 * stand_devs[:, 5], 'r--')
      plt.plot(t_span, x_hat_array[:, 5], 'k.')
      plt.ylabel('vz (km/s)')
      plt.xticks([min(t_span), max(t_span)], rotation=30, ha='right')
      ax = plt.gca()
      ax.set_xticklabels(['t_0', 't_f'])
      plt.ylim((-1.1*np.max(stand_devs[:, 5]), 1.1*np.max(stand_devs[:, 5])))

      plt.suptitle('State Errors' + ' $\hat{x}$' + ' and Covariance Bounds')
      plt.tight_layout()
      plt.subplots_adjust(top=.9)
      plt.savefig(dirIt+'/State Errors and Covariance Bounds.png', dpi=300, format='png')
    elif data_dict['acc_est'] == 'ON':
      plt.figure(1)
      plt.subplot(331)
      plt.plot(t_span, stand_devs[:, 0], 'r--', t_span, -1 * stand_devs[:, 0], 'r--')
      plt.plot(t_span, x_hat_array[:, 0], 'k.')
      plt.ylabel('x (km)')
      plt.xticks([])
      plt.title('Position')
      plt.ylim((-np.max(stand_devs[:, 0]), np.max(stand_devs[:, 0])))

      plt.subplot(334)
      plt.plot(t_span, stand_devs[:, 1], 'r--', t_span, -1 * stand_devs[:, 1], 'r--')
      plt.plot(t_span, x_hat_array[:, 1], 'k.')
      plt.ylabel('y (km)')
      plt.xticks([])
      plt.ylim((-np.max(stand_devs[:, 1]), np.max(stand_devs[:, 1])))

      plt.subplot(337)
      plt.plot(t_span, stand_devs[:, 2], 'r--', t_span, -1 * stand_devs[:, 2], 'r--')
      plt.plot(t_span, x_hat_array[:, 2], 'k.')
      plt.ylabel('z (km)')
      plt.xticks([min(t_span), max(t_span)], rotation=30, ha='right')
      ax = plt.gca()
      ax.set_xticklabels(['t_0', 't_f'])
      plt.ylim((-np.max(stand_devs[:, 2]), np.max(stand_devs[:, 2])))

      plt.subplot(332)
      plt.plot(t_span, stand_devs[:, 3], 'r--', t_span, -1 * stand_devs[:, 3], 'r--')
      plt.plot(t_span, x_hat_array[:, 3], 'k.')
      plt.ylabel('vx (km/s)')
      plt.xticks([])
      plt.title('Velocity')
      plt.ylim((-np.max(stand_devs[:, 3]), np.max(stand_devs[:, 3])))


      plt.subplot(335)
      plt.plot(t_span, stand_devs[:, 4], 'r--', t_span, -1 * stand_devs[:, 4], 'r--')
      plt.plot(t_span, x_hat_array[:, 4], 'k.')
      plt.ylabel('vy (km/s)')
      plt.xticks([])
      plt.ylim((-np.max(stand_devs[:, 4]), np.max(stand_devs[:, 4])))

      plt.subplot(338)
      plt.plot(t_span, stand_devs[:, 5], 'r--', t_span, -1 * stand_devs[:, 5], 'r--')
      plt.plot(t_span, x_hat_array[:, 5], 'k.')
      plt.ylabel('vz (km/s)')
      plt.xticks([min(t_span), max(t_span)], rotation=30, ha='right')
      ax = plt.gca()
      ax.set_xticklabels(['t_0', 't_f'])
      plt.ylim((-np.max(stand_devs[:, 5]), np.max(stand_devs[:, 5])))

      plt.subplot(333)
      plt.plot(t_span, stand_devs[:, 6], 'r--', t_span, -1 * stand_devs[:, 6], 'r--')
      plt.plot(t_span, x_hat_array[:, 6], 'k.')
      plt.ylabel('ax (km/s)')
      plt.xticks([])
      plt.title('Accelerations')
      plt.ylim((-1.1*np.max(stand_devs[:, 6]), 1.1*np.max(stand_devs[:, 6])))


      plt.subplot(336)
      plt.plot(t_span, stand_devs[:, 7], 'r--', t_span, -1 * stand_devs[:, 7], 'r--')
      plt.plot(t_span, x_hat_array[:, 7], 'k.')
      plt.ylabel('ay (km/s)')
      plt.xticks([])
      plt.ylim((-1.1*np.max(stand_devs[:, 7]), 1.1*np.max(stand_devs[:, 7])))

      plt.subplot(339)
      plt.plot(t_span, stand_devs[:, 8], 'r--', t_span, -1 * stand_devs[:, 8], 'r--')
      plt.plot(t_span, x_hat_array[:, 8], 'k.')
      plt.ylabel('az (km/s)')
      plt.xticks([min(t_span), max(t_span)], rotation=30, ha='right')
      ax = plt.gca()
      ax.set_xticklabels(['t_0', 't_f'])
      plt.ylim((-1.1*np.max(stand_devs[:, 8]), 1.1*np.max(stand_devs[:, 8])))

      plt.suptitle('State Errors' + ' $\hat{x}$' + ' and Covariance Bounds')
      plt.tight_layout()
      plt.subplots_adjust(top=.9)
      plt.savefig(dirIt+'/State Errors and Covariance Bounds.png', dpi=300, format='png')
    else:
      print 'Did not specify if plotting estimated accelerations was on or off'


    if data_dict['acc_est'] == 'OFF':
      plt.figure(111)
      plt.subplot(321)
      plt.plot(t_span, est_states[:, 0], 'b.')
      plt.ylabel('x (km)')
      plt.xticks([])
      plt.title('Position')

      plt.subplot(323)
      plt.plot(t_span, est_states[:, 1], 'b.')
      plt.ylabel('y (km)')
      plt.xticks([])

      plt.subplot(325)
      plt.plot(t_span, est_states[:, 2], 'b.')
      plt.ylabel('z (km)')
      plt.xticks([min(t_span), max(t_span)], rotation=30, ha='right')
      ax = plt.gca()
      ax.set_xticklabels(['t_0', 't_f'])

      plt.subplot(322)
      plt.plot(t_span, est_states[:, 3], 'b.')
      plt.ylabel('vx (km/s)')
      plt.xticks([])
      plt.title('Velocity')


      plt.subplot(324)
      plt.plot(t_span, est_states[:, 4], 'b.')
      plt.ylabel('vy (km/s)')
      plt.xticks([])

      plt.subplot(326)
      plt.plot(t_span, est_states[:, 5], 'b.')
      plt.ylabel('vz (km/s)')
      plt.xticks([min(t_span), max(t_span)], rotation=30, ha='right')
      ax = plt.gca()
      ax.set_xticklabels(['t_0', 't_f'])

      plt.suptitle('States' + ' $X$' )
      plt.tight_layout()
      plt.subplots_adjust(top=.9)
      plt.savefig(dirIt+'/States.png', dpi=300, format='png')
    elif data_dict['acc_est'] == 'ON':
      plt.figure(111)
      plt.subplot(331)
      plt.plot(t_span, est_states[:, 0], 'b.')
      plt.ylabel('x (km)')
      plt.xticks([])
      plt.title('Position')

      plt.subplot(334)
      plt.plot(t_span, est_states[:, 1], 'b.')
      plt.ylabel('y (km)')
      plt.xticks([])

      plt.subplot(337)
      plt.plot(t_span, est_states[:, 2], 'b.')
      plt.ylabel('z (km)')
      plt.xticks([min(t_span), max(t_span)], rotation=30, ha='right')
      ax = plt.gca()
      ax.set_xticklabels(['t_0', 't_f'])

      plt.subplot(332)
      plt.plot(t_span, est_states[:, 3], 'b.')
      plt.ylabel('vx (km/s)')
      plt.xticks([])
      plt.title('Velocity')

      plt.subplot(335)
      plt.plot(t_span, est_states[:, 4], 'b.')
      plt.ylabel('vy (km/s)')
      plt.xticks([])

      plt.subplot(338)
      plt.plot(t_span, est_states[:, 5], 'b.')
      plt.ylabel('vz (km/s)')
      plt.xticks([min(t_span), max(t_span)], rotation=30, ha='right')
      ax = plt.gca()
      ax.set_xticklabels(['t_0', 't_f'])

      plt.subplot(333)
      plt.plot(t_span, est_states[:, 6], 'b.')
      plt.ylabel('ax (km/s)')
      plt.ylim(-max(est_states[:,6])*3,max(est_states[:, 6])*3)
      plt.xticks([])
      plt.title('Accelerations')


      plt.subplot(336)
      plt.plot(t_span, est_states[:, 7], 'b.')
      plt.ylabel('ay (km/s)')
      plt.ylim(-max(est_states[:,7])*3,max(est_states[:, 7])*3)
      plt.xticks([])

      plt.subplot(339)
      plt.plot(t_span, est_states[:, 8], 'b.')
      plt.ylabel('az (km/s)')
      plt.xticks([min(t_span), max(t_span)], rotation=30, ha='right')
      plt.ylim(-max(est_states[:,8])*3,max(est_states[:, 8])*3)
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
    plt.plot(t_span, x_hat_array[:, 0])
    plt.ylabel('x (km)')
    plt.xticks([])
    plt.title('Position')
    # plt.ylim((-ymax, ymax))

    plt.subplot(323)
    plt.plot(t_span, x_hat_array[:, 1])
    plt.ylabel('y (km)')
    plt.xticks([])
    # plt.ylim((-ymax, ymax))

    plt.subplot(325)
    plt.plot(t_span, x_hat_array[:, 2])
    plt.ylabel('z (km)')
    plt.xticks([min(t_span), max(t_span)], rotation=30, ha='right')
    ax = plt.gca()
    ax.set_xticklabels(['t_0', 't_f'])
    # plt.ylim((-ymax, ymax))

    plt.subplot(322)
    plt.plot(t_span, x_hat_array[:, 3])
    plt.ylabel('vx (km/s)')
    plt.xticks([])
    plt.title('Velocity')
    # plt.ylim((-ymax, ymax))

    plt.subplot(324)
    plt.plot(t_span, x_hat_array[:, 4])
    plt.ylabel('vy (km/s)')
    plt.xticks([])
    # plt.ylim((-ymax, ymax))

    plt.subplot(326)
    plt.plot(t_span, x_hat_array[:, 5])
    plt.ylabel('vz (km/s)')
    plt.xticks([min(t_span), max(t_span)], rotation=30, ha='right')
    ax = plt.gca()
    ax.set_xticklabels(['t_0', 't_f'])
    # plt.ylim((-ymax, ymax))
    plt.suptitle('State Error Estimates ' + '$\hat{x}$')
    plt.tight_layout()
    plt.subplots_adjust(top=.9)
    plt.savefig(dirIt + '/State Error Estimates.png', dpi=300, format='png')

    plt.figure(3)
    plt.subplot(321)
    plt.plot(t_span, err[:, 0])
    plt.xticks([])
    plt.ylabel('x (km)')
    plt.title('Position')

    plt.subplot(323)
    plt.plot(t_span, err[:, 1])
    plt.xticks([])
    plt.ylabel('y (km)')

    plt.subplot(325)
    plt.plot(t_span, err[:, 2])
    plt.ylabel('z (km)')
    plt.xticks([min(t_span), max(t_span)], rotation=30, ha='right')
    ax = plt.gca()
    ax.set_xticklabels(['t_0', 't_f'])

    plt.subplot(322)
    plt.plot(t_span, err[:, 3])
    plt.ylabel('vx (km/s)')
    plt.xticks([])
    plt.title('Velocity')

    plt.subplot(324)
    plt.plot(t_span, err[:, 4])
    plt.ylabel('vy (km/s)')
    plt.xticks([])
    plt.subplot(326)
    plt.plot(t_span, err[:, 5])
    plt.ylabel('vz (km/s)')
    plt.xticks([min(t_span), max(t_span)], rotation=30, ha='right')
    ax = plt.gca()
    ax.set_xticklabels(['t_0', 't_f'])
    plt.suptitle('Reference Trajectory Error vs Spice Trajectory: Dynamics Errors')
    plt.tight_layout()
    plt.subplots_adjust(top=.9)
    plt.savefig(dirIt + '/Reference Trajectory Error.png', dpi=300, format='png')

    plt.figure(4)
    plt.subplot(121)
    plt.plot(t_span, prefits[:, 0], '.')
    plt.plot(t_span, 3*np.ones(len(t_span))*observation_uncertainty[0,0]**2, 'r--')
    plt.plot(t_span, -3*np.ones(len(t_span))*observation_uncertainty[0,0]**2, 'r--')
    plt.xticks([])
    plt.title('Pixel (-)')
    plt.ylabel('Residual')
    plt.xticks(t_span, rotation=90, ha='right')
    ax = plt.gca()
    ax.set_xticklabels(beacon_list)

    plt.subplot(122)
    plt.plot(t_span, prefits[:, 1], '.')
    plt.plot(t_span, 3*np.ones(len(t_span))*observation_uncertainty[1, 1], 'r--')
    plt.plot(t_span, -3*np.ones(len(t_span))*observation_uncertainty[1, 1], 'r--')
    plt.xticks([])
    plt.title('Line (-)')
    plt.xticks(t_span, rotation=90, ha='right')
    ax = plt.gca()
    ax.set_xticklabels(beacon_list)

    plt.suptitle('Prefit Observation Residuals')
    plt.tight_layout()
    plt.subplots_adjust(top=.9)
    plt.savefig(dirIt + '/Prefit Observation Residuals.png', dpi=300, format='png')

    plt.figure(5)
    plt.subplot(121)
    plt.plot(t_span, postfits[:, 0], '.')
    plt.plot(t_span, 3*np.ones(len(t_span))*observation_uncertainty[0,0], 'r--')
    plt.plot(t_span, -3*np.ones(len(t_span))*observation_uncertainty[0,0], 'r--')
    plt.xticks([])
    plt.ylabel('Residual')
    plt.title('Pixel (-)')
    plt.xticks(t_span, rotation=90, ha='right')
    ax = plt.gca()
    ax.set_xticklabels(beacon_list)

    plt.subplot(122)
    plt.plot(t_span, postfits[:, 1], '.')
    plt.plot(t_span, 3*np.ones(len(t_span))*observation_uncertainty[1, 1], 'r--')
    plt.plot(t_span, -3*np.ones(len(t_span))*observation_uncertainty[1, 1], 'r--')
    plt.xticks([])
    plt.title('Line (-)')
    plt.xticks(t_span, rotation=90, ha='right')
    ax = plt.gca()
    ax.set_xticklabels(beacon_list)


    plt.suptitle('Postfit Observation Residuals')
    plt.tight_layout()
    plt.subplots_adjust(top=.9)
    plt.savefig(dirIt + '/Postfit Observation Residuals.png', dpi=300, format='png')

    plt.figure(55)
    plt.subplot(121)
    plt.plot(t_span, postDelta[:, 0], '.')
    plt.xticks([])
    plt.ylabel('Residual')
    plt.title('Pixel (-)')
    plt.xticks(t_span, rotation=90, ha='right')
    ax = plt.gca()
    ax.set_xticklabels(beacon_list)

    plt.subplot(122)
    plt.plot(t_span, postDelta[:, 1], '.')
    plt.xticks([])
    plt.title('Line (-)')
    plt.xticks(t_span, rotation=90, ha='right')
    ax = plt.gca()
    ax.set_xticklabels(beacon_list)


    plt.suptitle('Postfit Changes')
    plt.tight_layout()
    plt.subplots_adjust(top=.9)
    plt.savefig(dirIt + '/Postfit Delta.png', dpi=300, format='png')

    plt.figure(7)
    plt.subplot(321)
    plt.plot(t_span, true_ephem['spacecraft'].T[:, 0] - est_states[:, 0])
    plt.ylabel('$\delta$ x (km)')
    plt.xticks([])
    plt.title('Position')
    # plt.ylim((-ymax, ymax))

    plt.subplot(323)
    plt.plot(t_span, true_ephem['spacecraft'].T[:, 1] - est_states[:, 1])
    plt.ylabel('$\delta$ y (km)')
    plt.xticks([])
    # plt.ylim((-ymax, ymax))

    plt.subplot(325)
    plt.plot(t_span, true_ephem['spacecraft'].T[:, 2] - est_states[:, 2])
    plt.ylabel('$\delta$ z (km)')
    plt.xticks([min(t_span), max(t_span)], rotation=30, ha='right')
    ax = plt.gca()
    ax.set_xticklabels(['t_0', 't_f'])
    # plt.ylim((-ymax, ymax))

    plt.subplot(322)
    plt.plot(t_span, true_ephem['spacecraft'].T[:, 3] - est_states[:, 3])
    plt.ylabel('$\delta$ vx (km/s)')
    plt.xticks([])
    plt.title('Velocity')
    # plt.ylim((-ymax, ymax))

    plt.subplot(324)
    plt.plot(t_span, true_ephem['spacecraft'].T[:, 4] - est_states[:, 4])
    plt.ylabel('$\delta$ vy (km/s)')
    plt.xticks([])
    # plt.ylim((-ymax, ymax))

    plt.subplot(326)
    plt.plot(t_span, true_ephem['spacecraft'].T[:, 5] - est_states[:, 5])
    plt.ylabel('$\delta$ vz (km/s)')
    plt.xticks([min(t_span), max(t_span)], rotation=30, ha='right')
    ax = plt.gca()
    ax.set_xticklabels(['t_0', 't_f'])
    # plt.ylim((-ymax, ymax))
    plt.suptitle('Estimated state deviation from true Spice data')
    plt.tight_layout()
    plt.subplots_adjust(top=.9)
    plt.savefig(dirIt + '/Truth-EstimatedState.png', dpi=300, format='png')



    plt.figure(77)

    PosErrorNorm = np.zeros(len(est_states[:,0]))
    VelErrorNorm = np.zeros(len(est_states[:, 0]))
    PosCovarNormMax =  np.zeros(len(est_states[:, 0]))
    PosCovarNormMin =  np.zeros(len(est_states[:, 0]))
    VelCovarNormMax =  np.zeros(len(est_states[:, 0]))
    VelCovarNormMin =  np.zeros(len(est_states[:, 0]))
    for i in range(len(est_states[:,0])):
        PosErrorNorm[i] = np.linalg.norm(true_ephem['spacecraft'].T[i, 0:3] - est_states[i, 0:3])
        VelErrorNorm[i] = np.linalg.norm(true_ephem['spacecraft'].T[i, 3:6] - est_states[i, 3:6])
        PosCovarNormMax[i] = PosErrorNorm[i] + np.sqrt(stand_devs[i,0]**2 + stand_devs[i,1]**2 + stand_devs[i,2]**2)/30.
        PosCovarNormMin[i] = PosErrorNorm[i] - np.sqrt(stand_devs[i,0]**2 + stand_devs[i,1]**2 + stand_devs[i,2]**2)/30.
        VelCovarNormMax[i] = VelErrorNorm[i] + np.sqrt(stand_devs[i,3]**2 + stand_devs[i,4]**2 + stand_devs[i,5]**2)/30.
        VelCovarNormMin[i] = VelErrorNorm[i] - np.sqrt(stand_devs[i,3]**2 + stand_devs[i,4]**2 + stand_devs[i,5]**2)/30.

    n_PosMin = PosErrorNorm.argmax()
    n_VelMin = VelErrorNorm.argmax()

    plt.subplot(121)
    plt.plot(t_span,PosErrorNorm, 'b')
    plt.plot(t_span,PosCovarNormMax, 'r--', label= '$\sigma$/10')
    plt.plot(t_span,PosCovarNormMin, 'r--')
    plt.ylabel('$\delta$ R (km)')
    plt.yticks(np.linspace(PosCovarNormMax.max(), PosCovarNormMin.min(), 12))
    plt.xticks([])
    ax = plt.gca()
    ax.set_xticklabels(['t_0', 't_f'])
    plt.title('Position Error')
    # plt.ylim((-ymax, ymax))

    plt.subplot(122)
    plt.plot(t_span, VelErrorNorm, 'b')
    plt.plot(t_span,VelCovarNormMax, 'r--',label= '$\sigma$/10')
    plt.plot(t_span,VelCovarNormMin, 'r--')
    plt.ylabel('$\delta$ V (km/s)')
    plt.yticks(np.linspace(VelCovarNormMax.max(), VelCovarNormMin.min(), 12))
    plt.xticks([])
    ax = plt.gca()
    ax.set_xticklabels(['t_0', 't_f'])
    plt.title('Velocity Error')
    # plt.ylim((-ymax, ymax))

    plt.suptitle('Estimated state Norm from true Spice data')
    plt.tight_layout()
    plt.subplots_adjust(top=.9)
    plt.savefig(dirIt + '/Truth-EstimatedStateNorms.png', dpi=300, format='png')

    plt.figure(6)
    plt.scatter(ref_state[:, 1], ref_state[:, 2], color='r')
    plt.scatter(true_ephem['spacecraft'].T[:, 1], true_ephem['spacecraft'].T[:, 2], color='b')
    plt.suptitle('Reference Trajectory (r) vs Spice Trajectory (b)')
    plt.savefig(dirIt + '/Trajectories.png', dpi=300, format='png')

    jet = plt.get_cmap('jet')
    colors = iter(jet(np.linspace(0, 1, extras['n_unique_beacons'])))

    jet = plt.get_cmap('jet')
    colors = iter(jet(np.linspace(0, 1, extras['n_unique_beacons'])))

    plt.figure(8)
    plt.subplot(321)
    plt.plot(t_span, stand_devs[:, 0], 'r-', t_span, P_array[:, 0, 0], 'b-')
    plt.ylabel('x (km)')
    plt.xticks([])
    ax = plt.gca()
    ax.set_yscale('log')
    plt.title('Position')

    plt.subplot(323)
    plt.plot(t_span, stand_devs[:, 1], 'r-', t_span, P_array[:, 1, 1], 'b-')
    plt.ylabel('y (km)')
    plt.xticks([])
    ax = plt.gca()
    ax.set_yscale('log')

    plt.subplot(325)
    plt.plot(t_span, stand_devs[:, 2], 'r-', t_span, P_array[:, 2, 2], 'b-')
    plt.ylabel('z (km)')
    plt.xticks([min(t_span), max(t_span)], rotation=30, ha='right')
    ax = plt.gca()
    ax.set_yscale('log')
    ax.set_xticklabels(['t_0', 't_f'])

    plt.subplot(322)
    plt.plot(t_span, stand_devs[:, 3], 'r-', t_span, P_array[:, 3, 3], 'b-')
    plt.ylabel('vx (km/s)')
    plt.xticks([])
    ax = plt.gca()
    ax.set_yscale('log')
    plt.title('Velocity')

    plt.subplot(324)
    plt.plot(t_span, stand_devs[:, 4], 'r-', t_span, P_array[:, 4, 4], 'b-')
    plt.ylabel('vy (km/s)')
    ax = plt.gca()
    ax.set_yscale('log')
    plt.xticks([])

    plt.subplot(326)
    plt.plot(t_span, stand_devs[:, 5], 'r-', t_span, P_array[:, 5, 5], 'b-')
    plt.ylabel('vz (km/s)')
    plt.xticks([min(t_span), max(t_span)], rotation=30, ha='right')
    ax = plt.gca()
    ax.set_yscale('log')
    ax.set_xticklabels(['t_0', 't_f'])
    plt.suptitle('Standard Deviation and Covariance')
    plt.tight_layout()
    plt.subplots_adjust(top=.9)
    plt.savefig(dirIt + '/Standard Deviation and Covariance.png', dpi=300, format='png')

    plt.figure(9)
    plt.subplot(121)
    plt.plot(t_span, Y_observed[:, 0] - Y_truth[:, 0])
    plt.xticks([])
    plt.title('Range (km)')
    plt.ylabel('Added Error')
    plt.xticks(t_span, rotation=90, ha='right')
    ax = plt.gca()
    ax.set_xticklabels(beacon_list)

    plt.subplot(122)
    plt.plot(t_span, Y_observed[:, 1] - Y_truth[:, 1])
    plt.xticks([])
    plt.title('Range Rate (km/s)')
    plt.xticks(t_span, rotation=90, ha='right')
    ax = plt.gca()
    ax.set_xticklabels(beacon_list)

    plt.suptitle('Observation Errors: $Y_{obs} - Y_{true}$')
    plt.tight_layout()
    plt.subplots_adjust(top=.9)
    plt.savefig(dirIt + '/Observation Errors.png', dpi=300, format='png')

    plt.figure(10)
    plt.subplot(321)
    plt.plot(t_span, err_hat[:, 0])
    plt.ylabel('x (km)')
    plt.xticks([])
    plt.title('Position')
    # plt.ylim((-ymax, ymax))

    plt.subplot(323)
    plt.plot(t_span, err_hat[:, 1])
    plt.ylabel('y (km)')
    plt.xticks([])
    # plt.ylim((-ymax, ymax))

    plt.subplot(325)
    plt.plot(t_span, err_hat[:, 2])
    plt.ylabel('z (km)')
    plt.xticks([min(t_span), max(t_span)], rotation=30, ha='right')
    ax = plt.gca()
    ax.set_xticklabels(['t_0', 't_f'])
    # plt.ylim((-ymax, ymax))

    plt.subplot(322)
    plt.plot(t_span, err_hat[:, 3])
    plt.ylabel('vx (km/s)')
    plt.xticks([])
    plt.title('Velocity')
    # plt.ylim((-ymax, ymax))

    plt.subplot(324)
    plt.plot(t_span, err_hat[:, 4])
    plt.ylabel('vy (km/s)')
    plt.xticks([])
    # plt.ylim((-ymax, ymax))

    plt.subplot(326)
    plt.plot(t_span, err_hat[:, 5])
    plt.ylabel('vz (km/s)')
    plt.xticks([min(t_span), max(t_span)], rotation=30, ha='right')
    ax = plt.gca()
    ax.set_xticklabels(['t_0', 't_f'])
    # plt.ylim((-ymax, ymax))
    plt.suptitle('Estimated State minus Spice Ephemeris')
    plt.tight_layout()
    plt.subplots_adjust(top=.9)
    plt.savefig(dirIt + '/State Deviation Est - Spice.png', dpi=300, format='png')

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
