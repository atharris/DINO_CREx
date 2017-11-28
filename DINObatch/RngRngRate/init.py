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
import pdb
filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))
path2 = os.path.dirname(os.path.abspath(filename))
bskName = 'Basilisk'
dinoName = 'DINO_CREx'
splitPath = path.split(dinoName)
dinoSpicePath = splitPath[0] + dinoName + '/DINObatch/SPICE/'
bskSpicePath = splitPath[0] + bskName + '/External/EphemerisData/'
bskPath = splitPath[0] + bskName + '/'
sys.path.append(bskPath + 'modules')
sys.path.append(bskPath + 'PythonModules')
sys.path.append(dinoSpicePath)


import pyswice
import numpy as np
import matplotlib.pyplot as plt
from batchFilter import run_batch
import data_generation as dg
################################################################################
#                  S E C O N D A R Y     F U N C T I O N S:
################################################################################

# -------------------------------------------------------------------------------
def norm(input):
    norm = np.sqrt(sum(np.square(input)))
    return norm

def writingText(itr, ref_state, est_state, true_ephem, extra_data, position_error , velocity_error):
    # calculate the difference between the perturbed reference and true trajectories: reference state errors
    err = ref_state[:, 0:6] - true_ephem['spacecraft'].T

    # compare the estimated and true trajectories: estimated state errors
    err_hat = est_state[:, 0:6] - true_ephem['spacecraft'].T

    ResultString = ''

    ResultString += '---------------------------------------------------' + '\n'
    ResultString += 'Iteration number '+ str(itr) + '\n'
    ResultString += '---------------------------------------------------'+ '\n'
    ResultString += '\n'
    ResultString += 'Estimated x_hat_0 = ' + str(extra_data['x_hat_0'])+ '\n'
    ResultString += 'Actual Error = ' + str(position_error) + str(velocity_error) + '\n'
    ResultString += '\n'

    ResultString += 'Last X Pos err = ' + str(err[-1, 0]) + '\n'
    ResultString += 'Last Y Pos err = ' + str(err[-1, 1]) + '\n'
    ResultString += 'Last Z Pos err = ' + str(err[-1, 2]) + '\n'
    ResultString += 'Last X Vel err = ' + str(err[-1, 3]) + '\n'
    ResultString += 'Last Y Vel err = ' + str(err[-1, 4]) + '\n'
    ResultString += 'Last Z Vel err = ' + str(err[-1, 5]) + '\n'
    ResultString += '\n'
    ResultString += 'Last X Pos err = ' + str(err_hat[-1, 0]) + '\n'
    ResultString += 'Last Y Pos err = ' + str(err_hat[-1, 1]) + '\n'
    ResultString += 'Last Z Pos err = ' + str(err_hat[-1, 2]) + '\n'
    ResultString += 'Last X Vel err = ' + str(err_hat[-1, 3]) + '\n'
    ResultString += 'Last Y Vel err = ' + str(err_hat[-1, 4]) + '\n'
    ResultString += 'Last Z Vel err = ' + str(err_hat[-1, 5]) + '\n'
    ResultString += '\n'

    print ResultString

    text_file = open('Batch_Iteration' + str(itr) + "/Batch" + str(itr) + ".txt", "w")
    text_file.write(ResultString)
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

def main():

    ###########################################
    #  S  P  I  C  E  C  O  D  E
    #  S  P  I  C  E  C  O  D  E
    ##########################################
    pyswice.furnsh_c(bskSpicePath + 'de430.bsp')
    pyswice.furnsh_c(dinoSpicePath + 'naif0011.tls')
    pyswice.furnsh_c(dinoSpicePath + 'DINO_kernel.bsp')
    # SP.furnsh('SPICE/jup310.bsp')
    # SP.furnsh('SPICE/mar097.bsp')

    DINO_kernel = dinoSpicePath + 'DINO_kernel.bsp'
    body_int = -100#SP.spkobj(DINO_kernel)
    body_id_str = str(body_int)

    # search_window = pyswice.new_doubleArray(2)
    # pyswice.spkcov_c(DINO_kernel, body_int, search_window)
    # list_of_events = pyswice.wnfetd_c(search_window, 0)
    # tBSP_Start = list_of_events[0]
    # tBSP_End = list_of_events[1]

    ###########################################
    # Initial condition for spacecraft
    # data = io.loadmat('saves/obsData.mat')
    # true_ephem = {}
    # reference of sun to sc
    # true_ephem['spacecraft'] = np.copy(data['stateS'])
    # # reference of sun to Earth
    # true_ephem['S2E'] = np.copy(data['stateE'])
    # # reference of sun to Mars
    # true_ephem['S2M'] = np.copy(data['stateM'])

    # time span
    # t_span = data['etT'].flatten()
    #Filtering End Epochs
    start_et = pyswice.new_doubleArray(1)
    end_et=pyswice.new_doubleArray(1)
    pyswice.utc2et_c('23 JUL 2020 17:00:00', start_et)
    pyswice.utc2et_c('30 JUL 2020 17:00:00', end_et)

    start_et = pyswice.doubleArray_getitem(start_et, 0)
    end_et = pyswice.doubleArray_getitem(end_et, 0)


    # extras dictionary for importing to functions
    extras = {}
    # body vector for SUN, EARTH, MARS
    # CODE RELIES ON SUN BEING INDEXED AS 0
    extras['bodies'] = ['SUN', '3', '399']

    # specify primary and secondary
    extras['primary'] = 0
    extras['secondary'] = [1, 2]

    # respective GP vector
    extras['mu'] = [1.32712428 * 10 ** 11, 3.986004415 * 10 ** 5, 4.305 * 10 ** 4]

    # abcorr for spkzer
    extras['abcorr'] = 'NONE'

    # reference frame
    extras['ref_frame'] = 'J2000'

    # SRP parameter
    extras['SRP'] = -20 * 4.57 * 10 ** (-6)

    # coefficient of reflectivity
    extras['cR'] = 1

    # number of observations per beacon until moving to the next
    extras['n_obs'] = 1

    # SNC coefficient
    extras['SNC'] = (2 * 10 ** (-4)) ** 3

    # Number of batch iterations
    extras['iterations'] = 5

    # Are we using the real dynamics for the ref or the trueData
    extras['realData']= 'OFF'

    # Initializing the error
    extras['x_hat_0'] = 0
    # Get Observation Times and Ephemerides
    true_ephem, t_span = dg.generate_data(sc_ephem_file=DINO_kernel,
                                          planet_beacons = ['earth','mars barycenter'],
                                          beacon_ids=[],
                                          n_observations=12,
                                          start_et=start_et,
                                          end_et=end_et,
                                          extras = extras,
                                          realData = extras['realData'])

    tt_switch = 5




    # number and keys of beacons. note that the true ephem is going to have one spot for the
    # sun, which in NOT a beacon
    beacon_names = true_ephem.keys()
    beacon_names.remove('spacecraft')
    extras['beacons'] = beacon_names
    extras['n_beacons'] = len(beacon_names)

    pdb.set_trace()

    ##################################################################################
    #
    # BLOCK A page 196
    #
    ##################################################################################

    # copy the initial conditions as the first sun to SC ref_states from the SPICE file
    IC = np.copy(true_ephem['spacecraft'][:, 0])

    print 'IC', IC
    spice_derived_state = pyswice.new_doubleArray(6)
    lt = pyswice.new_doubleArray(1)
    pyswice.spkezr_c(body_id_str, t_span[0], 'J2000', 'None', 'Sun', spice_derived_state, lt)

    # a priori uncertainty for the ref_states
    P_bar = np.zeros((IC.shape[0], IC.shape[0]))
    P_bar[0, 0] = 10000**2
    P_bar[1, 1] = 10000**2
    P_bar[2, 2] = 10000**2
    P_bar[3, 3] = .1**2
    P_bar[4, 4] = .1**2
    P_bar[5, 5] = .1**2

    # position_error = np.zeros(3)
    # velocity_error = np.zeros(3)

    # add uncertainty to the IC
    position_error = 5000 * np.divide(IC[0:3], norm(IC[0:3]))
    velocity_error = 0.05 * np.divide(IC[3:6], norm(IC[3:6]))

    IC += np.append(position_error, velocity_error)

    # uncertainty to be added in the form of noise to the measurables. Takes the form of variance
    observation_uncertainty = np.identity(2)
    observation_uncertainty[0, 0] = 1**2.
    observation_uncertainty[1, 1] = .005**2.

    # the initial STM is an identity matrix
    phi0 = np.identity(IC.shape[0])

    # initiate a priori deviation
    x_bar = np.zeros(IC.shape)

    # initiate a filter output dictionary
    filter_outputs = {}

    # run the filter and output the reference ref_states (including STMs), est states and extra data
    for itr in xrange(extras['iterations']):

        if itr > 0:
            IC = est_state[0, :]
            x_bar -= extra_data['x_hat_array'][0, :]

        # the arguments for the filter are the IC, the first STM, the time span, the observables
        # data dictionary, a priori uncertainty, and the measurables' uncertainty,
        # as well as any extras
        filter_inputs = (IC, phi0, t_span, true_ephem, P_bar, observation_uncertainty, x_bar, extras)
        # run filter function
        ref_state, est_state, extra_data = run_batch(filter_inputs)

        # save all outputs into the dictionary with a name associated with the iteration
        filter_outputs[str(itr)] = {}
        filter_outputs[str(itr)]['ref_state'] = ref_state
        filter_outputs[str(itr)]['est_state'] = est_state
        filter_outputs[str(itr)]['extra_data'] = extra_data

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
        writingText(itr+1, ref_state, est_state, true_ephem, extra_data, position_error , velocity_error)

        # calculate the difference between the perturbed reference and true trajectories: reference state errors
        err = ref_state[:, 0:6] - true_ephem['spacecraft'].T

        # compare the estimated and true trajectories: estimated state errors
        err_hat = est_state[:, 0:6] - true_ephem['spacecraft'].T

        x_hat_array = extra_data['x_hat_array']
        P_array = extra_data['P_array']
        prefits = extra_data['prefit residuals']
        postfits = extra_data['postfit residuals']
        beacon_list = extra_data['Y']['beacons']
        Y_observed = extra_data['Y']['data']
        Y_truth = extra_data['Y']['truth']

        stand_devs = 3. * np.array([np.sqrt(np.fabs(np.diag(P))) for P in P_array])

        plt.figure(1)
        plt.subplot(321)
        plt.plot(t_span, stand_devs[:, 0], 'r--', t_span, -1 * stand_devs[:, 0], 'r--')
        plt.plot(t_span, x_hat_array[:, 0], 'k.')
        plt.ylabel('x (km)')
        plt.xticks([])
        plt.title('Position')
        plt.ylim((-np.max(stand_devs[:, 0]), np.max(stand_devs[:, 0])))

        plt.subplot(323)
        plt.plot(t_span, stand_devs[:, 1], 'r--', t_span, -1 * stand_devs[:, 1], 'r--')
        plt.plot(t_span, x_hat_array[:, 1], 'k.')
        plt.ylabel('y (km)')
        plt.xticks([])
        plt.ylim((-np.max(stand_devs[:, 1]), np.max(stand_devs[:, 1])))

        plt.subplot(325)
        plt.plot(t_span, stand_devs[:, 2], 'r--', t_span, -1 * stand_devs[:, 2], 'r--')
        plt.plot(t_span, x_hat_array[:, 2], 'k.')
        plt.ylabel('z (km)')
        plt.xticks([min(t_span), max(t_span)], rotation=30, ha='right')
        ax = plt.gca()
        ax.set_xticklabels(['t_0', 't_f'])
        plt.ylim((-np.max(stand_devs[:, 2]), np.max(stand_devs[:, 2])))

        plt.subplot(322)
        plt.plot(t_span, stand_devs[:, 3], 'r--', t_span, -1 * stand_devs[:, 3], 'r--')
        plt.plot(t_span, x_hat_array[:, 3], 'k.')
        plt.ylabel('vx (km/s)')
        plt.xticks([])
        plt.title('Velocity')
        plt.ylim((-np.max(stand_devs[:, 3]), np.max(stand_devs[:, 3])))


        plt.subplot(324)
        plt.plot(t_span, stand_devs[:, 4], 'r--', t_span, -1 * stand_devs[:, 4], 'r--')
        plt.plot(t_span, x_hat_array[:, 4], 'k.')
        plt.ylabel('vy (km/s)')
        plt.xticks([])
        plt.ylim((-np.max(stand_devs[:, 4]), np.max(stand_devs[:, 4])))

        plt.subplot(326)
        plt.plot(t_span, stand_devs[:, 5], 'r--', t_span, -1 * stand_devs[:, 5], 'r--')
        plt.plot(t_span, x_hat_array[:, 5], 'k.')
        plt.ylabel('vz (km/s)')
        plt.xticks([min(t_span), max(t_span)], rotation=30, ha='right')
        ax = plt.gca()
        ax.set_xticklabels(['t_0', 't_f'])
        plt.ylim((-np.max(stand_devs[:, 5]), np.max(stand_devs[:, 5])))

        plt.suptitle('State Errors' + ' $\hat{x}$' + ' and Covariance Bounds')
        plt.tight_layout()
        plt.subplots_adjust(top=.9)
        plt.savefig(dirIt+'/State Errors and Covariance Bounds.png', dpi=300, format='png')

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
        plt.title('Range (km)')
        plt.ylabel('Residual')
        plt.xticks(t_span, rotation=90, ha='right')
        ax = plt.gca()
        ax.set_xticklabels(beacon_list)

        plt.subplot(122)
        plt.plot(t_span, prefits[:, 1], '.')
        plt.plot(t_span, 3*np.ones(len(t_span))*np.sqrt(observation_uncertainty[1, 1]**2.), 'r--')
        plt.plot(t_span, -3*np.ones(len(t_span))*np.sqrt(observation_uncertainty[1, 1]**2.), 'r--')
        plt.xticks([])
        plt.title('Range Rate (km/s)')
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
        plt.plot(t_span, 3*np.ones(len(t_span))*np.sqrt(observation_uncertainty[0,0]**2.), 'r--')
        plt.plot(t_span, -3*np.ones(len(t_span))*np.sqrt(observation_uncertainty[0,0]**2.), 'r--')
        plt.xticks([])
        plt.ylabel('Residual')
        plt.title('Range (km)')
        plt.xticks(t_span, rotation=90, ha='right')
        ax = plt.gca()
        ax.set_xticklabels(beacon_list)

        plt.subplot(122)
        plt.plot(t_span, postfits[:, 1], '.')
        plt.plot(t_span, 3*np.ones(len(t_span))*np.sqrt(observation_uncertainty[1, 1]**2.), 'r--')
        plt.plot(t_span, -3*np.ones(len(t_span))*np.sqrt(observation_uncertainty[1, 1]**2.), 'r--')
        plt.xticks([])
        plt.title('Range Rate (km/s)')
        plt.xticks(t_span, rotation=90, ha='right')
        ax = plt.gca()
        ax.set_xticklabels(beacon_list)


        plt.suptitle('Postfit Observation Residuals')
        plt.tight_layout()
        plt.subplots_adjust(top=.9)
        plt.savefig(dirIt + '/Postfit Observation Residuals.png', dpi=300, format='png')

        plt.figure(6)
        plt.scatter(ref_state[:, 1], ref_state[:, 2], color='r')
        plt.scatter(true_ephem['spacecraft'].T[:, 1], true_ephem['spacecraft'].T[:, 2], color='b')
        plt.suptitle('Reference Trajectory (r) vs Spice Trajectory (b)')
        plt.savefig(dirIt + '/Trajectories.png', dpi=300, format='png')

        jet = plt.get_cmap('jet')
        colors = iter(jet(np.linspace(0, 1, extras['n_beacons'])))

        jet = plt.get_cmap('jet')
        colors = iter(jet(np.linspace(0, 1, extras['n_beacons'])))

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


if __name__ == "__main__":
    main()
