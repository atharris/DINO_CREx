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

import pdb
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

    ResultString += 'Ref X Pos err = ' + str(err[-1, 0]) + '\n'
    ResultString += 'Ref Y Pos err = ' + str(err[-1, 1]) + '\n'
    ResultString += 'Ref Z Pos err = ' + str(err[-1, 2]) + '\n'
    ResultString += 'Ref X Vel err = ' + str(err[-1, 3]) + '\n'
    ResultString += 'Ref Y Vel err = ' + str(err[-1, 4]) + '\n'
    ResultString += 'Ref Z Vel err = ' + str(err[-1, 5]) + '\n'
    ResultString += '\n'
    ResultString += 'Est X Pos err = ' + str(err_hat[-1, 0]) + '\n'
    ResultString += 'Est Y Pos err = ' + str(err_hat[-1, 1]) + '\n'
    ResultString += 'Est Z Pos err = ' + str(err_hat[-1, 2]) + '\n'
    ResultString += 'Est X Vel err = ' + str(err_hat[-1, 3]) + '\n'
    ResultString += 'Est Y Vel err = ' + str(err_hat[-1, 4]) + '\n'
    ResultString += 'Est Z Vel err = ' + str(err_hat[-1, 5]) + '\n'
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

    # extras dictionary for importing to functions
    extras = {}

    ###########################################
    #
    #  S  P  I  C  E  C  O  D  E
    #  
    ##########################################

    # basic .bsp filename (generic, such as de430, etc)
    extras['basic_bsp']   = 'de430.bsp'
    # .bsp filename for mission
    extras['mission_bsp'] = 'DINO_kernel.bsp'
    # .tls filename 
    extras['tls']         = 'naif0011.tls'

    # prep pyswice for the extraction of initial data
    # is the only reason that we do this is for lines 165 and 166?
    pyswice.furnsh_c(bskSpicePath  + 'de430.bsp')
    pyswice.furnsh_c(dinoSpicePath + 'naif0011.tls')
    pyswice.furnsh_c(dinoSpicePath + 'DINO_kernel.bsp')

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
    end_et   = pyswice.doubleArray_getitem(end_et, 0)



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
    # A/M ratio multiplied by solar pressure constant at 1 AU with adjustments
    extras['SRP'] = 0.3**2/14. * 149597870.**2 * 1358. / 299792458. / 1000. # turboprop document Eq (64)

    # coefficient of reflectivity
    extras['cR'] = 1

    # number of observations per beacon until moving to the next
    extras['repeat_obs'] = 1

    # SNC coefficient
    extras['SNC'] = (2 * 10 ** (-4)) ** 3

    # Number of batch iterations
    extras['iterations'] = 3

    # Initializing the error
    extras['x_hat_0'] = 0


    ##################################################################################
    #
    # Camera/P&L Parameters
    #
    ##################################################################################

    # Focal Length (mm)
    extras['FoL'] = 100.
    extras['DCM_BI'] = np.eye(3)
    extras['DCM_TVB'] = np.eye(3)

    # Camera resolution (pixels)
    extras['resolution'] = [1024., 1024.]

    # width and height of pixels in camera
    extras['pixel_width'] = 5.
    extras['pixel_height'] = 5.

    # direction coefficient of pixel and line axes
    extras['pixel_direction'] = 1.
    extras['line_direction'] = 1.

    # Are we using the real dynamics for the ref or the trueData
    extras['realData']= 'ON'

    # Add anomaly detection parameters
    extras['anomaly']= False
    extras['anomaly_num'] = 0
    extras['anomaly_threshold'] = 4

    ##################################################################################

    # Get Observation Times and Ephemerides. This outputs a full data set that is not
    # parsed in any way. Ephemerides for all objects at all times are given.
    true_ephem, t_span = dg.generate_data(sc_ephem_file=DINO_kernel,
                                          planet_beacons = ['earth','mars barycenter'],
                                          beacon_ids=[],
                                          n_observations=24,
                                          start_et=start_et,
                                          end_et=end_et,
                                          extras = extras,
                                          realData = extras['realData'])

    tt_switch = 5



    # number and keys of beacons. note that the true ephem is going to have one spot for the
    # sun, which in NOT a beacon. These are used in beaconBinSPICE. 
    beacon_names = true_ephem.keys()
    beacon_names.remove('spacecraft')
    extras['unique_beacon_IDs'] = beacon_names
    extras['n_unique_beacons'] = len(beacon_names)

    ##################################################################################
    #
    # BLOCK A page 196
    #
    ##################################################################################

    # copy the initial conditions as the first sun to SC ref_states from the SPICE file
    IC = np.copy(true_ephem['spacecraft'][:, 0])

    print 'IC', IC

    # spice_derived_state is only referenced here. Should these be axed?
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

    # add uncertainty to the IC
    position_error = 1000 * np.divide(IC[0:3], norm(IC[0:3]))
    velocity_error = 0.01 * np.divide(IC[3:6], norm(IC[3:6]))

    IC[0:6] += np.append(position_error, velocity_error)

    # uncertainty to be added in the form of noise to the measurables. 
    # Takes the form of variance. Currently, the same value is used in both
    # the creation of the measurements as well as the weighting of the filter (W)
    observation_uncertainty = np.identity(2)
    observation_uncertainty[0, 0] = 0.2 ** 2
    observation_uncertainty[1, 1] = 0.2 ** 2

    # the initial STM is an identity matrix
    phi0 = np.identity(IC.shape[0])

    # initiate a priori deviation
    x_bar = np.zeros(IC.shape)

    # initiate a filter output dictionary
    filter_outputs = {}

    ##################################################################################
    #
    # Get the noisy observations
    #
    ##################################################################################
        
    # observation inputs
    obs_inputs = (true_ephem, observation_uncertainty, extras)

    # Get the observation data (obs_data). This dictionary contains the SPICE data
    # from which values are calculated (key = 'SPICE'), the true observations before
    # uncertainty is added (key = 'truth') and the measured observations (key = 'data').
    # These are the 'data' values that are now simulating an actual observation, 
    # and they are to be processed by the filter. 
    # The dictionary also contains the list of beacons by name and order of processing. 
    # This list of strings (key = 'beacons') is needed for 
    # the filter's own beacon position generator
    obs_data = getObs(obs_inputs)

    # create dictionary for observation data to be inputs in filter
    obs_filter = {}
    obs_filter['measurements'] = obs_data['data']
    obs_filter['beaconIDs']    = obs_data['beacons']

    ##################################################################################
    #
    # Run the Filter
    #
    ##################################################################################

   # run the filter and output the reference ref_states (including STMs), est states and extra data
    for itr in xrange(extras['iterations']):

        if itr > 0:
            IC = est_state[0, :]
            x_bar -= extra_data['x_hat_array'][0, :]

        if itr==0:
            extras['oldPost'] = np.zeros([len(t_span), 2])

        # the arguments for the filter are the IC, the first STM, the time span, the observables
        # data dictionary, a priori uncertainty, and the measurables' uncertainty,
        # as well as any extras
        filter_inputs = (IC, phi0, t_span, obs_filter,\
                         P_bar, observation_uncertainty, x_bar, extras)
        # run filter function
        ref_state, est_state, extra_data = run_batch(filter_inputs)
        extras['oldPost'] = extra_data['postfit residuals']

        # Check for anomaly:

        [anomaly_bool , anomaly_num] = extra_data['anomaly_detected']
        if anomaly_bool == True:
            print '**********************************************************'
            print 'Anomaly Detected - Estimates are not to be trusted'
            print '**********************************************************'
            print anomaly_num, 'Residuals out of bounds'
            return

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

        plot_data = extra_data

        plot_data['postfit delta']= extra_data['postfit changes']
        plot_data['states']      = est_state
        plot_data['truth']       = obs_data['truth']
        plot_data['beacon_list'] = obs_data['beacons']
        plot_data['t_span']      = t_span
        plot_data['dirIt']       = dirIt
        plot_data['err']         = err
        plot_data['err_hat']     = err_hat
        plot_data['obs_uncertainty'] = observation_uncertainty
        plot_data['ref_state']   = ref_state
        plot_data['true_ephem']  = true_ephem
        plot_data['extras']      = extras
        plot_data['acc_est']     = 'OFF'
        PF( plot_data )

    return


if __name__ == "__main__":
    main()
