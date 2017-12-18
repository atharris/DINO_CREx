import sys, os, inspect
import matplotlib.pyplot as plt
from numpy import linalg as la
import numpy as np
import math

# filename = inspect.getframeinfo(inspect.currentframe()).filename
# path = os.path.dirname(os.path.abspath(filename))
bskName = 'Basilisk'
bskPath = '../..' + '/' + bskName + '/'
sys.path.append(bskPath + 'modules')
sys.path.append(bskPath + 'PythonModules')
sys.path.append('../dinoModels/SimCode/opnavCamera/')
sys.path.append('../dinoModels/SimCode/opnavCamera/dependencies')
sys.path.append('../dinoModels/fswAlgorithms/imageProcessing/dependencies/')
sys.path.append('../DINObatch/pixelAndLine/commonFunctions/')
sys.path.append('../DINObatch/pixelAndLine/unmodeledAcc/')
sys.path.append('../DINObatch/SPICE/')

import BSK_plotting as BSKPlt
from DINO_scenarioUtils import *

try:
    import macros as mc
    import unitTestSupport as sp
    import orbitalMotion as om
    import RigidBodyKinematics as rbk
except ImportError:
    import Basilisk.utilities.macros as mc
    import Basilisk.utilities.unitTestSupport as sp
    import Basilisk.utilities.orbitalMotion as om
    import Basilisk.utilities.RigidBodyKinematics as rbk
    from Basilisk.fswAlgorithms import *

import camera
import imageProcessingExecutive as ip
from initBatch import initBatchFnc

# ------------------------------------- DATA LOGGING ------------------------------------------------------ #

def log_DynCelestialOutputs(TheDynSim, samplingTime):
    TheDynSim.TotalSim.logThisMessage(TheDynSim.DynClass.scObject.scStateOutMsgName, samplingTime)
    TheDynSim.TotalSim.logThisMessage(TheDynSim.DynClass.earthGravBody.bodyInMsgName, samplingTime)
    TheDynSim.TotalSim.logThisMessage(TheDynSim.DynClass.marsGravBody.bodyInMsgName, samplingTime)
    TheDynSim.TotalSim.logThisMessage(TheDynSim.DynClass.sunGravBody.bodyInMsgName, samplingTime)
    TheDynSim.TotalSim.logThisMessage(TheDynSim.DynClass.moonGravBody.bodyInMsgName, samplingTime)
    TheDynSim.TotalSim.logThisMessage(TheDynSim.DynClass.gyroModel.OutputDataMsg, samplingTime)
    TheDynSim.TotalSim.logThisMessage(TheDynSim.DynClass.starTracker.outputStateMessage, samplingTime)
    #    for ind in range(0,len(TheDynSim.DynClass.beaconList)):
    #        TheDynSim.TotalSim.logThisMessage(TheDynSim.DynClass.beaconList[ind].scStateOutMsgName, samplingTime)

    return


def log_DynOutputs(TheBSKSim, samplingTime):
    TheBSKSim.TotalSim.logThisMessage(TheBSKSim.DynClass.scObject.scStateOutMsgName, samplingTime)
    TheBSKSim.TotalSim.logThisMessage(TheBSKSim.DynClass.simpleNavObject.outputAttName, samplingTime)
    return


def log_aekfOutputs(TheBskSim, samplingTime):
    print "Att output msg name:", TheBskSim.FSWClass.attFilter.outputMsgName
    print "Att filter msg name:", TheBskSim.FSWClass.attFilter.filterMsgName
    TheBskSim.TotalSim.logThisMessage(TheBskSim.FSWClass.attFilter.outputMsgName, samplingTime)
    TheBskSim.TotalSim.logThisMessage(TheBskSim.FSWClass.attFilter.filterMsgName, samplingTime)

    return


def log_FSWOutputs(TheBSKSim, samplingTime):
    TheBSKSim.TotalSim.logThisMessage(TheBSKSim.FSWClass.attErrorConfig.outputDataName, samplingTime)
    TheBSKSim.TotalSim.logThisMessage(TheBSKSim.FSWClass.mrpControlConfig.outputDataName, samplingTime)
    TheBSKSim.TotalSim.logThisMessage(TheBSKSim.FSWClass.attGuideConfig.outputDataName, samplingTime)
    TheBSKSim.TotalSim.logThisMessage(TheBSKSim.FSWClass.attGuideConfig.inputCelMessName, samplingTime)
    TheBSKSim.TotalSim.logThisMessage(TheBSKSim.FSWClass.attGuideConfig.inputSecMessName, samplingTime)
    return


# ------------------------------------- DATA PULLING ------------------------------------------------------ #

def pull_DynCelestialOutputs(TheDynSim, plots=True):
    r_BN = TheDynSim.pullMessageLogData(TheDynSim.DynClass.scObject.scStateOutMsgName + '.r_BN_N', range(3))
    r_earth = TheDynSim.pullMessageLogData(TheDynSim.DynClass.earthGravBody.bodyInMsgName + '.PositionVector', range(3))
    r_sun = TheDynSim.pullMessageLogData(TheDynSim.DynClass.sunGravBody.bodyInMsgName + '.PositionVector', range(3))
    r_mars = TheDynSim.pullMessageLogData(TheDynSim.DynClass.marsGravBody.bodyInMsgName + '.PositionVector', range(3))
    r_moon = TheDynSim.pullMessageLogData(TheDynSim.DynClass.moonGravBody.bodyInMsgName + '.PositionVector', range(3))
    r_beacons = []
    #    for ind in range(0,len(TheDynSim.DynClass.beaconList)):
    #        r_beacons.append(TheDynSim.pullMessageLogData(TheDynSim.DynClass.beaconList[ind].scStateOutMsgName + '.r_BN_N', range(3)))

    # Print Dyn Celestial Outputs
    print '\n\n'
    print 'CELESTIAL:'
    print 'r_BN = ', la.norm(r_BN[-1:, 1:])
    print 'r_earth = ', la.norm(r_earth[-1:, 1:])
    print 'r_sun = ', la.norm(r_sun[-1:, 1:])
    print 'r_mars = ', la.norm(r_mars[-1:, 1:])
    print 'r_moon = ', la.norm(r_moon[-1:, 1:])
    for ind in range(0, len(r_beacons)):
        print 'r_beacon' + str(ind) + ':', la.norm(r_beacons[ind][-1:, 1:])

    dict_data_color = {
        'moon': [r_moon, 'cyan'],
        'earth': [r_earth, 'dodgerblue'],
        'mars': [r_mars, 'r'],
        'sun': [r_sun, 'orange'],
        # 'r_beacon_0': [r_beacons[0], 'g'],
        # 'r_beacon_1': [r_beacons[1], 'g'],
        # 'r_beacon_2': [r_beacons[2], 'g'],
        # 'r_beacon_3': [r_beacons[3], 'g'],
        # 'r_beacon_4': [r_beacons[4], 'g'],
        # 'r_beacon_5': [r_beacons[5], 'g'],
        # 'r_beacon_6': [r_beacons[6], 'g'],
        # 'r_beacon_7': [r_beacons[7], 'g'],
        # 'r_beacon_8': [r_beacons[8], 'g']
    }
    BSKPlt.plot_multi_orbit_0(dict_data_color)
    BSKPlt.plot_spacecraft_orbit_0(dict_data_color, r_BN)
    if plots == True:
        BSKPlt.plot_multi_orbit_0(dict_data_color)
        BSKPlt.plot_spacecraft_orbit_0(dict_data_color, r_BN)

    sc_dict_data_color = {
        'moon': [r_moon, 'cyan'],
        'earth': [r_earth, 'dodgerblue'],
        'mars': [r_mars, 'r'],
        #        'r_beacon_0': [r_beacons[0], 'g'],
        # 'r_beacon_1': [r_beacons[1], 'g'],
        # 'r_beacon_2': [r_beacons[2], 'g'],
        # 'r_beacon_3': [r_beacons[3], 'g'],
        # 'r_beacon_4': [r_beacons[4], 'g'],
        # 'r_beacon_5': [r_beacons[5], 'g'],
        # 'r_beacon_6': [r_beacons[6], 'g'],
        # 'r_beacon_7': [r_beacons[7], 'g'],
        # 'r_beacon_8': [r_beacons[8], 'g']
    }
    if plots == True:
        BSKPlt.plot_spacecraft_orbit(sc_dict_data_color, r_BN)
    print "Please make this work."
    return r_sun, r_earth, r_moon, r_mars, r_beacons


def pull_DynOutputs(TheBSKSim, plots=True):
    # Pull Dyn Outputs
    r_BN = TheBSKSim.pullMessageLogData(TheBSKSim.DynClass.scObject.scStateOutMsgName + '.r_BN_N', range(3))
    v_BN = TheBSKSim.pullMessageLogData(TheBSKSim.DynClass.scObject.scStateOutMsgName + '.v_BN_N', range(3))
    sigma_BN = TheBSKSim.pullMessageLogData(TheBSKSim.DynClass.simpleNavObject.outputAttName + ".sigma_BN", range(3))
    omega_BN_B = TheBSKSim.pullMessageLogData(TheBSKSim.DynClass.simpleNavObject.outputAttName + ".omega_BN_B",
                                              range(3))

    # Print Dyn Outputs
    if plots == True:
        print '\n\n'
        print 'DYNAMICS:'
        print 'sigma_BN = ', sigma_BN[-3:, 1:], '\n'
        print 'omega_BN_B = ', omega_BN_B[-3:, 1:], '\n'
        testRBN, testVBN = define_dino_earthSOI()
        print "Final position:", r_BN[-1, 1:] / 1000.0, '\n'
        print "Desired final position:", testRBN / 1000.0, '\n'
        print "Position percent error:", np.subtract(testRBN, r_BN[-1, 1:]) / la.norm(testRBN) * 100.0
        print "Final velocity", v_BN[-1, 1:] / 1000.0, '\n'
        print "Des final velocity:", testVBN / 1000.0, '\n'
        print "Velocity percent error:", np.subtract(testVBN, v_BN[-1, 1:]) / la.norm(testVBN) * 100.0

        print "Final time in sec since sim start:", (r_BN[-1, 0]) / 1.0e9
        # Plot Relevant Dyn Outputs
        # BSKPlt.plot_orbit(r_BN)
        BSKPlt.plot_rotationalNav(sigma_BN, omega_BN_B)

    return r_BN, v_BN, sigma_BN, omega_BN_B


def pull_senseOutputs(TheBSKSim, plots=True):
    # Pull Dyn Outputs
    beta_tilde_BN = TheBSKSim.pullMessageLogData(TheBSKSim.DynClass.starTracker.outputStateMessage + '.qInrtl2Case',
                                                 range(4))
    omega_tilde_BN = TheBSKSim.pullMessageLogData(TheBSKSim.DynClass.gyroModel.OutputDataMsg + '.AngVelPlatform',
                                                  range(3))

    numInds = beta_tilde_BN.shape[0]
    sigma_tilde_BN = np.zeros([numInds, 4])

    for ind in range(1, beta_tilde_BN.shape[0]):
        sigma_tilde_BN[ind, 0] = beta_tilde_BN[ind, 0]
        sigma_tilde_BN[ind, 1:] = rbk.EP2MRP(beta_tilde_BN[ind, 1:])

    if plots == True:
        # Print Dyn Outputs
        print '\n\n'
        print 'DYNAMICS:'
        print 'sigma_tilde_BN = ', sigma_tilde_BN[-3:, 1:], '\n'
        print 'omega_tilde_BN = ', omega_tilde_BN[-3:, 1:], '\n'
        testRBN, testVBN = define_dino_earthSOI()

        # Plot Relevant Dyn Outputs
        # BSKPlt.plot_orbit(r_BN)
        BSKPlt.plot_rotationalNav(sigma_tilde_BN, omega_tilde_BN)
    return sigma_tilde_BN, omega_tilde_BN


def pull_aekfOutputs(TheBSKSim, plots=True):
    # Pull Dyn OutputsTheBSKSim.DynClass.simpleNavObject.outputAttName + ".sigma_BN", range(3)
    sigma_hat_BN = TheBSKSim.pullMessageLogData(TheBSKSim.FSWClass.attFilter.outputMsgName + '.sigma_BN', range(3))
    omega_hat_BN = TheBSKSim.pullMessageLogData(TheBSKSim.FSWClass.attFilter.outputMsgName + '.omega_BN_B', range(3))

    # Pull true outputs in order to debug and plot fitler plots
    # Note that these are taken at a 5x higher sample rate. Can't figure out why...
    sigma_BN = TheBSKSim.pullMessageLogData(TheBSKSim.DynClass.simpleNavObject.outputAttName + ".sigma_BN", range(3))
    omega_BN_B = TheBSKSim.pullMessageLogData(TheBSKSim.DynClass.simpleNavObject.outputAttName + ".omega_BN_B",
                                              range(3))

    # Pull filter msg data
    covarLog1 = TheBSKSim.pullMessageLogData(TheBSKSim.FSWClass.attFilter.filterMsgName + '.sigma_BN', range(3))
    covarLog2 = TheBSKSim.pullMessageLogData(TheBSKSim.FSWClass.attFilter.filterMsgName + '.vehSunPntBdy', range(3))
    postFitLog = TheBSKSim.pullMessageLogData(TheBSKSim.FSWClass.attFilter.filterMsgName + '.omega_BN_B', range(3))

    # Print Dyn Outputs
    print '\n\n'
    print 'DYNAMICS:'
    print 'sigma_hat_BN = ', sigma_hat_BN[-3:, 1:], '\n'
    print 'omega_hat_BN = ', omega_hat_BN[-3:, 1:], '\n'
    testRBN, testVBN = define_dino_earthSOI()

    sigma_err = np.copy(sigma_hat_BN)
    omega_err = np.copy(omega_hat_BN)
    sigma_err[:, 1:4] = np.array(sigma_hat_BN)[:, 1:4] - np.array(sigma_BN)[:-1, 1:4]
    omega_err[:, 1:4] = np.array(omega_hat_BN)[:, 1:4] - np.array(omega_BN_B)[:-1, 1:4]
    covarLog = np.zeros([np.shape(sigma_err)[0], 7])
    covarLog[:, 0:4] = covarLog1
    covarLog[:, 4:7] = covarLog2[:, 1:4]
    # Plot Relevant Dyn Outputs
    # BSKPlt.plot_orbit(r_BN)
    BSKPlt.plot_rotationalNav(sigma_hat_BN, omega_hat_BN)
    BSKPlt.plot_filterOut(sigma_err, omega_err, covarLog)
    BSKPlt.plot_filterPostFits(postFitLog, 0.001 * np.identity(3))

    return sigma_hat_BN, omega_hat_BN


def pull_FSWOutputs(TheBSKSim, plots=True):
    sigma_RN = TheBSKSim.pullMessageLogData(TheBSKSim.FSWClass.attGuideConfig.outputDataName + ".sigma_RN", range(3))
    omega_RN_N = TheBSKSim.pullMessageLogData(TheBSKSim.FSWClass.attGuideConfig.outputDataName + ".omega_RN_N",
                                              range(3))

    pos1 = TheBSKSim.pullMessageLogData(TheBSKSim.FSWClass.attGuideConfig.inputNavDataName + ".r_BN_N", range(3))
    pos2 = TheBSKSim.pullMessageLogData(TheBSKSim.FSWClass.attGuideConfig.inputCelMessName + ".PositionVector",
                                        range(3))
    pos3 = TheBSKSim.pullMessageLogData(TheBSKSim.FSWClass.attGuideConfig.inputSecMessName + ".PositionVector",
                                        range(3))

    sigma_BR = TheBSKSim.pullMessageLogData(TheBSKSim.FSWClass.attErrorConfig.outputDataName + ".sigma_BR", range(3))
    omega_BR_B = TheBSKSim.pullMessageLogData(TheBSKSim.FSWClass.attErrorConfig.outputDataName + ".omega_BR_B",
                                              range(3))
    Lr = TheBSKSim.pullMessageLogData(TheBSKSim.FSWClass.mrpControlConfig.outputDataName + ".torqueRequestBody",
                                      range(3))

    if plots == True:
        print '\n\n'
        print 'FSW:'
        print 'sigma_RN = ', sigma_RN[-3:, 1:], '\n'
        print 'sigma_BR = ', sigma_BR[-3:, 1:], '\n'
        print 'Lr = ', Lr[:9, 1:], '\n'

        BSKPlt.plot_trackingError(sigma_BR, omega_BR_B)
        BSKPlt.plot_attitudeGuidance(sigma_RN, omega_RN_N)
        BSKPlt.plot_controlTorque(Lr)

    return Lr


# ------------------------------------- DATA HANDLING ------------------------------------------------------ #
def scenario_logResults(TheBSKSim, samplingTime):
    log_DynOutputs(TheBSKSim, samplingTime)
    log_FSWOutputs(TheBSKSim, samplingTime)


def scenario_plotResults(TheBSKSim):
    pull_DynOutputs(TheBSKSim)
    pull_FSWOutputs(TheBSKSim)