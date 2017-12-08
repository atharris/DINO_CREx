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

import BSK_plotting as BSKPlt

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
    TheBSKSim.TotalSim.logThisMessage(TheBSKSim.FSWClass.trackingErrorData.outputDataName, samplingTime)
    TheBSKSim.TotalSim.logThisMessage(TheBSKSim.FSWClass.mrpFeedbackData.outputDataName, samplingTime)
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

    return r_sc, r_sun, r_earth, r_moon, r_mars, r_beacons

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
    sigma_hat_BN = TheBSKSim.pullMessageLogData(TheBSKSim.FSWClass.attFilter.outputMsgName+ '.sigma_BN', range(3))
    omega_hat_BN = TheBSKSim.pullMessageLogData(TheBSKSim.FSWClass.attFilter.outputMsgName+ '.omega_BN_B', range(3))

    # Pull true outputs in order to debug and plot fitler plots
    # Note that these are taken at a 5x higher sample rate. Can't figure out why...
    sigma_BN = TheBSKSim.pullMessageLogData(TheBSKSim.DynClass.simpleNavObject.outputAttName + ".sigma_BN", range(3))
    omega_BN_B = TheBSKSim.pullMessageLogData(TheBSKSim.DynClass.simpleNavObject.outputAttName + ".omega_BN_B", range(3))

    # Pull filter msg data
    covarLog1 = TheBSKSim.pullMessageLogData(TheBSKSim.FSWClass.attFilter.filterMsgName+ '.sigma_BN', range(3))
    covarLog2 = TheBSKSim.pullMessageLogData(TheBSKSim.FSWClass.attFilter.filterMsgName+ '.vehSunPntBdy', range(3))
    postFitLog = TheBSKSim.pullMessageLogData(TheBSKSim.FSWClass.attFilter.filterMsgName+ '.omega_BN_B', range(3))

    # Print Dyn Outputs
    print '\n\n'
    print 'DYNAMICS:'
    print 'sigma_hat_BN = ', sigma_hat_BN[-3:, 1:], '\n'
    print 'omega_hat_BN = ', omega_hat_BN[-3:, 1:], '\n'
    testRBN, testVBN = define_dino_earthSOI()

    sigma_err = np.copy(sigma_hat_BN)
    omega_err =  np.copy(omega_hat_BN)
    sigma_err[:,1:4] = np.array(sigma_hat_BN)[:,1:4]-np.array(sigma_BN)[:-1,1:4]
    omega_err[:,1:4] = np.array(omega_hat_BN)[:,1:4]-np.array(omega_BN_B)[:-1,1:4]
    covarLog = np.zeros([np.shape(sigma_err)[0],7])
    covarLog[:,0:4] = covarLog1
    covarLog[:,4:7] = covarLog2[:,1:4]
    # Plot Relevant Dyn Outputs
    # BSKPlt.plot_orbit(r_BN)
    BSKPlt.plot_rotationalNav(sigma_hat_BN, omega_hat_BN)
    BSKPlt.plot_filterOut(sigma_err, omega_err, covarLog)
    BSKPlt.plot_filterPostFits(postFitLog, 0.001*np.identity(3))

    return sigma_hat_BN, omega_hat_BN

def pull_FSWOutputs(TheBSKSim, plots=True):
    sigma_RN = TheBSKSim.pullMessageLogData(TheBSKSim.FSWClass.trackingErrorData.inputRefName + ".sigma_RN", range(3))
    omega_RN_N = TheBSKSim.pullMessageLogData(TheBSKSim.FSWClass.trackingErrorData.inputRefName + ".omega_RN_N",
                                              range(3))
    sigma_BR = TheBSKSim.pullMessageLogData(TheBSKSim.FSWClass.trackingErrorData.outputDataName + ".sigma_BR", range(3))
    omega_BR_B = TheBSKSim.pullMessageLogData(TheBSKSim.FSWClass.trackingErrorData.outputDataName + ".omega_BR_B",
                                              range(3))
    Lr = TheBSKSim.pullMessageLogData(TheBSKSim.FSWClass.mrpFeedbackData.outputDataName + ".torqueRequestBody",
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


# ------------------------------------- SCENARIOS ------------------------------------------------------ #

def define_default_orbit():
    """
    Defines a default orbit for generic scenario
    :return: oe [structure] orbit elements set.
    :oe.a: Semi-major axis.
    :oe.e: Eccentricity
    :oe.i: Inclination
    :oe.Omega: Right Ascension of Ascending Node (raan)
    :oe.omega: Argument of Periapse (w)
    :oe.f: True Anomaly
    """
    oe = om.ClassicElements()
    oe.a = 10000.0 * 1000.0  # meters
    oe.e = 0.2  # 0.01
    oe.i = 0.0 * mc.D2R
    oe.Omega = 0.0 * mc.D2R
    oe.omega = 0.0 * mc.D2R
    oe.f = 280.0 * mc.D2R
    return oe


def define_dino_postTMI():
    sun_RN = np.array([49229056.29654553, -143897342.4521846, 9308.31227068603]) * 1000.
    sun_VN = np.array([33.23061608351387, 13.97929484231197, -0.5046624215893907]) * 1000.
    return sun_RN, sun_VN


def define_dino_earthSOI():
    sun_RN = np.array([53850779.24415602, -141892420.9599383, -61159.29323671013]) * 1000.
    sun_BN = np.array([32.95818246003602, 14.73457852769852, -0.5045251942794007]) * 1000.
    return sun_RN, sun_BN


def basicOrbit_dynScenario(TheDynSim):
    """
    Executes a default scenario for stand-alone dynamic simulations
    :params: TheDynSim: instantiation of class DINO_DynSim
    :return: None
    """

    # Initialize Simulation
    TheDynSim.InitializeSimulation()

    # Set up the orbit using classical orbit elements
    oe = define_default_orbit()
    mu = TheDynSim.DynClass.earthGravBody.mu
    rN, vN = om.elem2rv(mu, oe)
    om.rv2elem(mu, rN, vN)

    orbPeriod = period = 2 * np.pi * np.sqrt((oe.a ** 3.) / mu)

    # Log data for post-processing and plotting
    simulationTime = mc.sec2nano(orbPeriod)
    numDataPoints = int(orbPeriod)
    samplingTime = simulationTime / (numDataPoints - 1)
    log_DynOutputs(TheDynSim, samplingTime)

    # Initialize Spacecraft States within the state manager (after initialization)
    posRef = TheDynSim.DynClass.scObject.dynManager.getStateObject("hubPosition")
    velRef = TheDynSim.DynClass.scObject.dynManager.getStateObject("hubVelocity")
    sigmaRef = TheDynSim.DynClass.scObject.dynManager.getStateObject("hubSigma")
    omegaRef = TheDynSim.DynClass.scObject.dynManager.getStateObject("hubOmega")
    posRef.setState(sp.np2EigenVectorXd(rN))  # r_BN_N [m]
    velRef.setState(sp.np2EigenVectorXd(vN))  # r_BN_N [m]
    sigmaRef.setState([[0.1], [0.2], [-0.3]])  # sigma_BN_B
    omegaRef.setState([[0.001], [-0.01], [0.03]])  # omega_BN_B [rad/s]

    # Configure a simulation stop time time and execute the simulation run
    TheDynSim.ConfigureStopTime(simulationTime)
    TheDynSim.ExecuteSimulation()

    # Pull data for post-processing and plotting
    pull_DynOutputs(TheDynSim)
    plt.show()


def multiOrbitBeacons_dynScenario(TheDynSim):
    """
    Executes a default scenario for stand-alone dynamic simulations
    :params: TheDynSim: instantiation of class DINO_DynSim
    :return: None
    """
    # Log data for post-processing and plotting
    #   Set length of simulation in nanoseconds from the simulation start.
    simulationTime = mc.sec2nano(50.)
    #   Set the number of data points to be logged, and therefore the sampling frequency
    numDataPoints = 50
    samplingTime = simulationTime / (numDataPoints - 1)
    log_DynOutputs(TheDynSim, samplingTime)
    log_DynCelestialOutputs(TheDynSim, samplingTime)

    # Initialize Simulation
    TheDynSim.InitializeSimulation()

    # Set up the orbit using classical orbit elements
    # oe = define_default_orbit()
    mu = TheDynSim.DynClass.mu
    rN, vN = define_dino_postTMI()
    om.rv2elem(mu, rN, vN)

    # Initialize Spacecraft States within the state manager (after initialization)
    posRef = TheDynSim.DynClass.scObject.dynManager.getStateObject("hubPosition")
    velRef = TheDynSim.DynClass.scObject.dynManager.getStateObject("hubVelocity")
    sigmaRef = TheDynSim.DynClass.scObject.dynManager.getStateObject("hubSigma")
    omegaRef = TheDynSim.DynClass.scObject.dynManager.getStateObject("hubOmega")

    #   Set the spacecraft initial position, velocity, attitude parameters
    posRef.setState(sp.np2EigenVectorXd(rN))  # r_BN_N [m]
    velRef.setState(sp.np2EigenVectorXd(vN))  # r_BN_N [m]
    sigmaRef.setState([[0.1], [0.2], [-0.3]])  # sigma_BN_B
    omegaRef.setState([[0.001], [-0.01], [0.03]])  # omega_BN_B [rad/s]

    #   Set up beacon orbits using COEs loaded from a file
    ephemFile = open("observations_log.txt", 'r')

    lines = ephemFile.readlines()
    beaconOEs = []

    for ind in range(1, 1):
        vals = lines[ind].split(',')
        desVals = vals[8:]
        newOE = om.ClassicElements()
        newOE.a = float(desVals[0]) * 1000.0
        newOE.e = float(desVals[1])
        newOE.i = float(desVals[2])
        newOE.omega = float(desVals[3])
        newOE.Omega = float(desVals[4])
        newOE.f = float(desVals[5])
        beaconOEs.append(newOE)
        rN, vN = om.elem2rv(mu, newOE)
        print "Beacon", ind, " rN:", rN, " vN:", vN

        posRef = TheDynSim.DynClass.beaconList[ind - 1].dynManager.getStateObject("hubPosition")
        velRef = TheDynSim.DynClass.beaconList[ind - 1].dynManager.getStateObject("hubVelocity")
        posRef.setState(sp.np2EigenVectorXd(rN))  # r_BN_N [m]
        velRef.setState(sp.np2EigenVectorXd(vN))  # r_BN_N [m]

    ephemFile.close()

    # Configure a simulation stop time time and execute the simulation run
    TheDynSim.ConfigureStopTime(simulationTime)
    TheDynSim.ExecuteSimulation()

    # Pull data for post-processing and plotting
    r_BN, v_BN, sigma_BN, omega_BN_B = pull_DynOutputs(TheDynSim)
    sigma_tilde_BN, omega_tilde_BN = pull_senseOutputs(TheDynSim)
    r_sun, r_earth, r_moon, r_mars, r_beacons = pull_DynCelestialOutputs(TheDynSim)

    ##  Post-Process sim data using camera, image processing, batch filter DINO modules

    earth = camera.beacon()
    earth.r_eq = 6378.137
    earth.id = 'Earth'
    earth.albedo = 0.434

    mars = camera.beacon()
    mars.r_eq = 3396.2
    mars.id = 'Mars'
    mars.albedo = 0.17

    moon = camera.beacon()
    moon.r_eq = 1738.1
    moon.id = 'Earth'
    moon.albedo = 0.12

    beacons = [earth, moon, mars]

    # need loop to define asteroids, too

    

    # can kill these once I change the way camera is initialized
    takeImage = 0
    scState = -1
    scDCM = -1

    cam = camera.camera(
        0.01,  # detector_height
        0.01,  # detector_width
        0.05,  # focal_length
        512,  # resolution_height
        512,  # resolution_width
        np.identity(3),  # body2cameraDCM
        1000,  # maximum magnitude
        -1000,  # minimum magnitude (for debugging)
        qe,
        tc,
        1,
        0.01 ** 2,  # effective area in m^2
        100,  # dark current in electrons per second
        100,  # std for read noise in electrons
        1000,  # bin size
        2 ** 32,  # max bin depth
        1,
        0.01,  # simulation timestep
        scState,  # position state of s/c
        scDCM,  # intertal 2 body DCM for s/c
        beacons,  # bodies to track in images
        takeImage,  # takeImage message
        db='../dinoModels/SimCode/opnavCamera/db/tycho.db'  # stellar database
    )

    # this is spoofing the output of the nav exec
    # telling the camera when to take an image.
    takeImage = np.zeros(len(r_BN))
    takeImage[0] = 1
    takeImage[20] = 1
    takeImage[30] = 1
    takeImage[40] = 1
    takeImage[47] = 1

    lastTakeImage = 0
    for i in range(0, len(r_BN)):
        cam.scState = r_BN[i][1:4]
        earth.state = r_earth[i][1:4]
        moon.state = r_moon[i][1:4]
        mars.state = r_mars[i][1:4]
        # also need a loop here for
        # updating beacon position once they're added
        cam.scDCM = rbk.MRP2C(sigma_BN[i][1:4])

        # test that forces camera to point at cental star
        # in orion's belt
        # cam.scDCM = rbk.euler3212C(
        #     np.array([
        #         np.deg2rad(84.05338572),
        #         np.deg2rad(-1.20191725),
        #         np.deg2rad(0)]))

        cam.takeImage = takeImage[i]
        cam.imgTime = r_BN[i][0]
        cam.updateState()
    imgTimes = []
    detectorArrays = []
    imgPos = []
    imgMRP = []
    imgBeaconPos = []

    detectorArrays = []
    imgTimes = []
    imgPos = []
    imgMRP = []
    imgBeaconPos = []

    for i in range(0,len(cam.images)):
        detectorArrays.append(cam.images[i].detectorArray)
        imgTimes.append(cam.images[i].imgTime)
        imgPos.append(cam.images[i].imgPos)
        imgMRP.append(rbk.C2MRP(cam.images[i].imgDCM))
        imgBeaconPos.append(cam.images[i].imgBeaconPos)

        plt.figure()
        plt.imshow(cam.images[i].detectorArray)
    
    print('########################### END Image Generation ###########################')
    import pdb 
    pdb.set_trace()   


    # Run the Image Processing Module

    # required parameters from defineParameters function
    camParamIP = ipParam[0]
    beaconIDs = ipParam[1]
    beaconRadius = ipParam[2]

    imgTimesFound = []
    beaconIDsFound = []
    beaconPLFound = []
    imgMRPFound = []                # will have 'None' entries when not able to detect enough objects
    imgMRPFoundPassThrough = []     # identical attitude as input into the image processing module

    for indList in range(len(imgTimes)):
        currentBeaconIDs, currentPL, currentMRP = ip.imageProcessing(detectorArrays[indList],
                                                                     camParamIP,
                                                                     imgPos[indList],
                                                                     imgMRP[indList],
                                                                     imgBeaconPos[indList],
                                                                     beaconIDs,
                                                                     beaconRadius,
                                                                     makePlots=False,
                                                                     debugMode=True)

        for indBeacon in range(len(currentBeaconIDs)):
            imgTimesFound.append(imgTimes[indList])
            beaconIDsFound.append(currentBeaconIDs[indBeacon])
            beaconPLFound.append(currentPL[indBeacon])

            # pass through attitude estimate for navigation module
            imgMRPFoundPassThrough.append(imgMRP[indList])

            # attitude output of image processing logged for informational purposes only
            # (nav module to use sim attitude filter output)

            if currentMRP is not None:
                imgMRPFound.append(currentMRP)

        print '\nImage Processing Output: '
        print 'Image#: ', indList

        print '\nFound Beacon IDs, P/L, MRP'
        print currentBeaconIDs, currentPL, currentMRP
        print '\nInitial Estimate MRP: ', imgMRP[indList]
        print 'Initial Estimate DCM: '
        print cam.images[i].imgDCM


    # Generate inputs for navigation modulec
    numNavInputs = len(imgTimesFound)
    imgTimesNav = np.reshape(imgTimesFound, (numNavInputs, 1))
    beaconIDsNav = np.reshape(beaconIDsFound, (numNavInputs, 1))
    beaconPLNav = np.reshape(beaconPLFound, (numNavInputs, 2))
    print 'imgMRPFOundPassThrough: ', imgMRPFoundPassThrough
    imgMRPNav = np.reshape(imgMRPFoundPassThrough, (numNavInputs, 3))

    print beaconIDsNav
    print beaconPLNav
    print imgMRPNav

    import pdb
    pdb.set_trace()

    # Run the Navigation Module
    print('########################### END Image Processing ###########################')

    pdb.set_trace()
    plt.show()


    # Run the Image Processing Module

    imgTimesFound = []
    beaconIDsFound = []
    beaconPLFound = []
    imgMRPFound = []                # will have 'None' entries when not able to detect enough objects
    imgMRPFoundPassThrough = []     # identical attitude as input into the image processing module

    for indList in range(len(imgTimes)):
        currentBeaconIDs, currentPL, currentMRP = ip.imageProcessing(detectorArrays[indList],
                                                                    cameraParametersIP,
                                                                    imgPos[indList],
                                                                    imgMRP[indList],
                                                                    imgBeaconPos[indList],
                                                                    beaconIDs,
                                                                    beaconRadius)
        if currentBeaconIDs is not None:
            for indBeacon in range(len(currentBeaconIDs)):
                imgTimesFound.append(imgTimes[indList])
                beaconIDsFound.append(currentBeaconIDs[indBeacon])
                beaconPLFound.append(currentPL[indBeacon])

                # pass through attitude estimate for navigation module
                imgMRPFoundPassThrough.append(imgMRP[indList])

                # attitude output of image processing logged for informational purposes only
                # (nav module to use sim attitude filter output)
                imgMRPFound.append(currentMRP)


    # Generate inputs for navigation modulec
    numNavInputs = len(imgTimesFound)
    imgTimesNav = np.reshape(imgTimesFound, (numNavInputs, 1))
    beaconIDsNav = np.reshape(beaconIDsFound, (numNavInputs, 1))
    beaconPLNav = np.reshape(beaconPLFound, (numNavInputs, 2))
    imgMRPNav = np.reshape(imgMRPFoundPassThrough, (numNavInputs, 3))


    # Run the Navigation Module



def attFilter_dynScenario(TheDynSim):
    """
     Executes a default scenario for stand-alone dynamic simulations
    :params: TheDynSim: instantiation of class DINO_DynSim
    :return: None
    """
    # Log data for post-processing and plotting
    #   Set length of simulation in nanoseconds from the simulation start.
    simulationTime = mc.sec2nano(100)
    #   Set the number of data points to be logged, and therefore the sampling frequency
    numDataPoints = 10000
    samplingTime = simulationTime / (numDataPoints - 1)
    log_DynOutputs(TheDynSim, samplingTime)
    log_DynCelestialOutputs(TheDynSim, samplingTime)
    log_aekfOutputs(TheDynSim, samplingTime)

    # Initialize Simulation
    TheDynSim.InitializeSimulationAndDiscover()

    # Set up the orbit using classical orbit elements
    #oe = define_default_orbit()
    mu = TheDynSim.DynClass.mu
    rN, vN = define_dino_postTMI()
    om.rv2elem(mu, rN, vN)

    # Initialize Spacecraft States within the state manager (after initialization)
    posRef = TheDynSim.DynClass.scObject.dynManager.getStateObject("hubPosition")
    velRef = TheDynSim.DynClass.scObject.dynManager.getStateObject("hubVelocity")
    sigmaRef = TheDynSim.DynClass.scObject.dynManager.getStateObject("hubSigma")
    omegaRef = TheDynSim.DynClass.scObject.dynManager.getStateObject("hubOmega")

    #   Set the spacecraft initial position, velocity, attitude parameters
    posRef.setState(sp.np2EigenVectorXd(rN))  # r_BN_N [m]
    velRef.setState(sp.np2EigenVectorXd(vN))  # r_BN_N [m]
    sigmaRef.setState([[0.1], [0.2], [-0.3]])  # sigma_BN_B
    omegaRef.setState([[0.001], [-0.01], [0.03]])  # omega_BN_B [rad/s]

    # Configure a simulation stop time time and execute the simulation run
    TheDynSim.ConfigureStopTime(simulationTime)
    TheDynSim.ExecuteSimulation()

    # Pull data for post-processing and plotting
    pull_DynOutputs(TheDynSim)
    pull_senseOutputs(TheDynSim)
    pull_aekfOutputs(TheDynSim)
    #pull_DynCelestialOutputs(TheDynSim)
    plt.show()

def opnavCamera_dynScenario(TheDynSim):
    """
    Executes a default scenario for stand-alone camera simulation
    :params: TheDynSim: instantiation of class DINO_DynSim
    :return: None
    """
    # Log data for post-processing and plotting
    #   Set length of simulation in nanoseconds from the simulation start.
    simulationTime = mc.sec2nano(1000)
    #   Set the number of data points to be logged, and therefore the sampling frequency
    numDataPoints = 10000
    samplingTime = simulationTime / (numDataPoints - 1)
    log_DynOutputs(TheDynSim, samplingTime)
    log_DynCelestialOutputs(TheDynSim, samplingTime)

    # Initialize Simulation
    TheDynSim.InitializeSimulation()

    # Set up the orbit using classical orbit elements
    # oe = define_default_orbit()
    mu = TheDynSim.DynClass.mu
    rN, vN = define_dino_postTMI()
    om.rv2elem(mu, rN, vN)

    # Initialize Spacecraft States within the state manager (after initialization)
    posRef = TheDynSim.DynClass.scObject.dynManager.getStateObject("hubPosition")
    velRef = TheDynSim.DynClass.scObject.dynManager.getStateObject("hubVelocity")
    sigmaRef = TheDynSim.DynClass.scObject.dynManager.getStateObject("hubSigma")
    omegaRef = TheDynSim.DynClass.scObject.dynManager.getStateObject("hubOmega")

    #   Set the spacecraft initial position, velocity, attitude parameters
    posRef.setState(sp.np2EigenVectorXd(rN))  # r_BN_N [m]
    velRef.setState(sp.np2EigenVectorXd(vN))  # r_BN_N [m]
    sigmaRef.setState([[0.1], [0.2], [-0.3]])  # sigma_BN_B
    omegaRef.setState([[0.001], [-0.01], [0.03]])  # omega_BN_B [rad/s]

    # Configure a simulation stop time time and execute the simulation run
    TheDynSim.ConfigureStopTime(simulationTime)
    TheDynSim.ExecuteSimulation()

    # Pull data for post-processing and plotting
    pull_DynOutputs(TheDynSim)
    pull_senseOutputs(TheDynSim)
    # pull_DynCelestialOutputs(TheDynSim)
    plt.show()


def defineParameters(
    camResolution, 
    camFocalLength, 
    camSensorSize, 
    beacons,
    tc,
    qe,
    lambdaBinSize,
    effectiveArea,
    darkCurrent,
    readSTD,
    binSize,
    maxBinDepth,
    psfSTD,
    simTimeStep
    ):
    """
    Generates formatted inputs for camera, image processing, and navigation modules
    :params: camResolution      : (horizontal x vertical) camera resolution
    :params: camFocalLength     : camera focal length [m]
    :params: camSensorSize      : (horizontal x vertical) camera sensor size [m]
    :params: beacons            : N length list of beacon objects
    :params: tc                 : transmission curve dictionary (see SERs 4.3/4.3b)
    :params: qe                 : transmission curve dictionary (see SERs 4.3/4.3b)
    :params: lambdaBinSize      : bin size for lambda functions [nm]
    :params: effectiveArea      : effective area of camera [m^2]
    :params: darkCurrent        : dark current [electrons/s/pixel]
    :params: readSTD            : standard deviation of read noise [electrons/pixel]
    :params: binSize            : bin size [DN]
    :params: maxBinDepth        : saturation depth [DN]
    :params: psfSTD             : point spred funtion standard deviation [pixels]
    :params: simTimeStep        : simulation time step [s]

    :return: camInputs          : list of inputs for camera module
    :return: ipInputs           : list of inputs for image processing
    :return: navInputs          : lsit of inputs for navigation module
    """

    beaconIDs = []
    beaconRadius = []
    for each in beacons:
        beaconIDs.append(each.id)
        beaconRadius.append(each.r_eq)


    #init values for camera that will be set later.
    scState = -1
    scDCM = -1
    takeImage = 0


    cam = camera.camera(
        camSensorSize[0], #detector width in m
        camSensorSize[1], #detector height in m
        camFocalLength,     #focal lenght in m
        camResolution[0], #detector resolution (width direction)
        camResolution[1], #detector resolution (height direction)
        np.identity(3),  #body2cameraDCM
        1000,            #maximum magnitude (for debugging)
        -1000,           #minimum magnitude (for debugging)
        qe,              #quantum efficiency dictionary
        tc,              #transmission curve dictionary
        lambdaBinSize,   #lambda bin size
        effectiveArea,   #effective area in m^2
        darkCurrent,     #dark current in electrons per second
        readSTD,         #std for read noise in electrons
        binSize,         #bin size
        maxBinDepth,     #max bin depth
        psfSTD,          #std for psf
        simTimeStep,     #simulation timestep
        scState,         # position state of s/c
        scDCM,           # intertal 2 body DCM for s/c
        beacons,         # bodies to track in images
        takeImage,       # takeImage message
        db='../dinoModels/SimCode/opnavCamera/db/tycho.db'  # stellar database
    )


    # Camera Module Parameter Creation

    camInputs = cam

    # Image Processing Module Parameter Creation
    ipCamParam = {}
    ipCamParam['resolution'] = (cam.resolutionWidth,cam.resolutionHeight)
    ipCamParam['focal length'] = cam.focalLength
    ipCamParam['sensor size'] = (cam.detectorWidth,cam.detectorHeight)
    ipCamParam['pixel size'] = (
        cam.detectorWidth/cam.resolutionWidth, 
        cam.detectorHeight/cam.resolutionHeight)
    ipCamParam['field of view'] = (
        cam.angularWidth,
        cam.angularHeight)
    ipInputs = [ipCamParam, beaconIDs, beaconRadius]



    # Nav Module Parameter Creation

    navParams = {}

    # SPICE Parameters

    # basic .bsp filename (generic, such as de430, etc)
    navParams['basic_bsp']   = 'de430.bsp'
    # .bsp filename for mission
    navParams['mission_bsp'] = 'DINO_kernel.bsp'
    # .tls filename 
    navParams['tls']         = 'naif0011.tls'
    # abcorr for spkzer
    navParams['abcorr'] = 'NONE'
    # reference frame
    navParams['ref_frame'] = 'J2000'

    # Force Parameters
    
    #   Gravity
    # body vector for primary and secondary gravitational bodies
    navParams['bodies']    = ['SUN', '3', '399']
    # specify primary and secondary indices
    navParams['primary']   = 0
    navParams['secondary'] = [1, 2]
    # respective GP vector
    navParams['mu']        = [1.32712428 * 10 ** 11, 3.986004415 * 10 ** 5, 4.305 * 10 ** 4]
    #   SRP
    # A/M ratio multiplied by solar pressure constant at 1 AU with adjustments
    # Turboprop document Eq (64)
    navParams['SRP']       = 0.3**2/14. * 149597870.**2 * 1358. / 299792458. / 1000. 
    # coefficient of reflectivity
    navParams['cR']        = 1.

    # Camera/P&L Parameters

    # Focal Length (mm)
    navParams['FoL']             = 100.
    # default inertial to camera transformation matrices
    navParams['DCM_BI']          = np.eye(3)
    navParams['DCM_TVB']         = np.eye(3)
    # Camera resolution (pixels)
    navParams['resolution']      = [1024., 1024.]
    # width and height of pixels in camera
    navParams['pixel_width']     = 5.
    navParams['pixel_height']    = 5.
    # direction coefficient of pixel and line axes
    navParams['pixel_direction'] = 1.
    navParams['line_direction']  = 1.

    navInputs = navParams

    return camInputs, ipInputs, navInputs
