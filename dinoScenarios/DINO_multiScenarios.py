import sys, os, inspect
import matplotlib.pyplot as plt
from numpy import linalg as la
import numpy as np

# filename = inspect.getframeinfo(inspect.currentframe()).filename
# path = os.path.dirname(os.path.abspath(filename))
bskName = 'Basilisk'
bskPath = '../..' + '/' + bskName + '/'
sys.path.append(bskPath + 'modules')
sys.path.append(bskPath + 'PythonModules')
sys.path.append('../dinoModels/SimCode/opnavCamera/')
sys.path.append('../dinoModels/SimCode/opnavCamera/dependencies')


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
#    print "Att filter msg name:", TheBskSim.FSWClass.attFilter.outputMsgName
#    TheBskSim.TotalSim.logThisMessage(TheBskSim.FSWClass.attFilter.outputMsgName, samplingTime)
    return

def log_FSWOutputs(TheBSKSim, samplingTime):
    TheBSKSim.TotalSim.logThisMessage(TheBSKSim.FSWClass.trackingErrorData.outputDataName, samplingTime)
    TheBSKSim.TotalSim.logThisMessage(TheBSKSim.FSWClass.mrpFeedbackData.outputDataName, samplingTime)
    return

# ------------------------------------- DATA PULLING ------------------------------------------------------ #

def pull_DynCelestialOutputs(TheDynSim):
    r_sc = TheDynSim.pullMessageLogData(TheDynSim.DynClass.scObject.scStateOutMsgName + '.r_BN_N', range(3))
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
    print 'r_sc = ', la.norm(r_sc[-1:, 1:])
    print 'r_earth = ', la.norm(r_earth[-1:, 1:])
    print 'r_sun = ', la.norm(r_sun[-1:, 1:])
    print 'r_mars = ', la.norm(r_mars[-1:, 1:])
    print 'r_moon = ', la.norm(r_moon[-1:, 1:])
    for ind in range(0,len(r_beacons)):
        print 'r_beacon'+str(ind)+':', la.norm(r_beacons[ind][-1:,1:])

    dict_data_color = {
        'moon': [r_moon, 'cyan'],
        'earth': [r_earth, 'dodgerblue'],
        'mars': [r_mars, 'r'],
        'sun': [r_sun, 'orange'],
        #'r_beacon_0': [r_beacons[0], 'g'],
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
    BSKPlt.plot_spacecraft_orbit_0(dict_data_color, r_sc)


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
    if plots==True:
        BSKPlt.plot_spacecraft_orbit(sc_dict_data_color, r_sc)
    return r_sun, r_earth, r_moon, r_mars, r_beacons



def pull_DynOutputs(TheBSKSim):
    # Pull Dyn Outputs
    r_BN = TheBSKSim.pullMessageLogData(TheBSKSim.DynClass.scObject.scStateOutMsgName + '.r_BN_N', range(3))
    v_BN = TheBSKSim.pullMessageLogData(TheBSKSim.DynClass.scObject.scStateOutMsgName + '.v_BN_N', range(3))
    sigma_BN = TheBSKSim.pullMessageLogData(TheBSKSim.DynClass.simpleNavObject.outputAttName + ".sigma_BN", range(3))
    omega_BN_B = TheBSKSim.pullMessageLogData(TheBSKSim.DynClass.simpleNavObject.outputAttName + ".omega_BN_B", range(3))

    # Print Dyn Outputs
    print '\n\n'
    print 'DYNAMICS:'
    print 'sigma_BN = ', sigma_BN[-3:, 1:], '\n'
    print 'omega_BN_B = ', omega_BN_B[-3:, 1:], '\n'
    testRBN, testVBN = define_dino_earthSOI()
    print "Final position:", r_BN[-1,1:]/1000.0, '\n'
    print "Desired final position:", testRBN/1000.0, '\n'
    print "Position percent error:", np.subtract(testRBN, r_BN[-1, 1:])/la.norm(testRBN) * 100.0
    print "Final velocity", v_BN[-1,1:]/1000.0,  '\n'
    print "Des final velocity:", testVBN/1000.0, '\n'
    print "Velocity percent error:", np.subtract(testVBN, v_BN[-1,1:])/la.norm(testVBN) * 100.0

    print "Final time in sec since sim start:", (r_BN[-1,0])/1.0e9
    # Plot Relevant Dyn Outputs
    #BSKPlt.plot_orbit(r_BN)
    BSKPlt.plot_rotationalNav(sigma_BN, omega_BN_B)

    return r_BN, v_BN, sigma_BN, omega_BN_B



def pull_senseOutputs(TheBSKSim):
    # Pull Dyn Outputs
    beta_tilde_BN = TheBSKSim.pullMessageLogData(TheBSKSim.DynClass.starTracker.outputStateMessage + '.qInrtl2Case', range(4))
    omega_tilde_BN = TheBSKSim.pullMessageLogData(TheBSKSim.DynClass.gyroModel.OutputDataMsg + '.AngVelPlatform', range(3))

    numInds = beta_tilde_BN.shape[0]
    sigma_tilde_BN = np.zeros([numInds,4])

    for ind in range(1,beta_tilde_BN.shape[0]):
        sigma_tilde_BN[ind,0] = beta_tilde_BN[ind,0]
        sigma_tilde_BN[ind,1:] = rbk.EP2MRP(beta_tilde_BN[ind,1:])


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

def pull_aekfOutputs(TheBSKSim):
    # Pull Dyn Outputs
    sigma_hat_BN = TheBSKSim.pullMessageLogData(TheBSKSim.FSWClass.attFilter.outputMsgName+ '.sigma_BN', range(3))
    omega_hat_BN = TheBSKSim.pullMessageLogData(TheBSKSim.FSWClass.attFilter.outputMsgName+ '.omega_BN_B', range(3))


    # Print Dyn Outputs
    print '\n\n'
    print 'DYNAMICS:'
    print 'sigma_hat_BN = ', sigma_hat_BN[-3:, 1:], '\n'
    print 'omega_hat_BN = ', omega_hat_BN[-3:, 1:], '\n'
    testRBN, testVBN = define_dino_earthSOI()

    # Plot Relevant Dyn Outputs
    # BSKPlt.plot_orbit(r_BN)
    BSKPlt.plot_rotationalNav(sigma_hat_BN, omega_hat_BN)

    return sigma_hat_BN, omega_hat_BN

def pull_FSWOutputs(TheBSKSim):
    sigma_RN = TheBSKSim.pullMessageLogData(TheBSKSim.FSWClass.trackingErrorData.inputRefName + ".sigma_RN", range(3))
    omega_RN_N = TheBSKSim.pullMessageLogData(TheBSKSim.FSWClass.trackingErrorData.inputRefName + ".omega_RN_N", range(3))
    sigma_BR = TheBSKSim.pullMessageLogData(TheBSKSim.FSWClass.trackingErrorData.outputDataName + ".sigma_BR", range(3))
    omega_BR_B = TheBSKSim.pullMessageLogData(TheBSKSim.FSWClass.trackingErrorData.outputDataName + ".omega_BR_B", range(3))
    Lr = TheBSKSim.pullMessageLogData(TheBSKSim.FSWClass.mrpFeedbackData.outputDataName + ".torqueRequestBody", range(3))

    print '\n\n'
    print 'FSW:'
    print 'sigma_RN = ', sigma_RN[-3:, 1:], '\n'
    print 'sigma_BR = ', sigma_BR[-3:, 1:], '\n'
    print 'Lr = ', Lr[:9, 1:], '\n'

    BSKPlt.plot_trackingError(sigma_BR, omega_BR_B)
    BSKPlt.plot_attitudeGuidance(sigma_RN, omega_RN_N)
    BSKPlt.plot_controlTorque(Lr)

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
    sun_RN = np.array([49229056.29654553   ,      -143897342.4521846    ,    9308.31227068603 ]) *1000.
    sun_VN = np.array([33.23061608351387   ,      13.97929484231197   ,      -0.5046624215893907]) * 1000.
    return sun_RN, sun_VN

def define_dino_earthSOI():
    sun_RN = np.array([53850779.24415602    ,     -141892420.9599383  ,      -61159.29323671013]) * 1000.
    sun_BN = np.array([32.95818246003602    ,     14.73457852769852   ,      -0.5045251942794007]) * 1000.
    return sun_RN, sun_BN

def basicOrbit_dynScenario(TheDynSim):
    """
    Executes a default scenario for stand-alone dynamic simulations
    :params: TheDynSim: instantiation of class DINO_DynSim
    :return: None
    """
    # Log data for post-processing and plotting
    simulationTime = mc.sec2nano(139643.532)
    numDataPoints = 100
    samplingTime = simulationTime / (numDataPoints - 1)
    log_DynOutputs(TheDynSim, samplingTime)

    # Initialize Simulation
    TheDynSim.InitializeSimulation()

    # Set up the orbit using classical orbit elements
    oe = define_default_orbit()
    mu = TheDynSim.DynClass.earthGravBody.mu
    rN, vN = om.elem2rv(mu, oe)
    om.rv2elem(mu, rN, vN)

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
    simulationTime = mc.sec2nano(10)
    #   Set the number of data points to be logged, and therefore the sampling frequency
    numDataPoints = 100000
    samplingTime = simulationTime / (numDataPoints - 1)
    log_DynOutputs(TheDynSim, samplingTime)
    log_DynCelestialOutputs(TheDynSim, samplingTime)

    # Initialize Simulation
    TheDynSim.InitializeSimulation()

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

    #   Set up beacon orbits using COEs loaded from a file
    ephemFile = open("observations_log.txt", 'r')

    lines = ephemFile.readlines()
    beaconOEs = []

    for ind in range(1,1):
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

        posRef = TheDynSim.DynClass.beaconList[ind-1].dynManager.getStateObject("hubPosition")
        velRef = TheDynSim.DynClass.beaconList[ind-1].dynManager.getStateObject("hubVelocity")
        posRef.setState(sp.np2EigenVectorXd(rN))  # r_BN_N [m]
        velRef.setState(sp.np2EigenVectorXd(vN))  # r_BN_N [m]

    ephemFile.close()

    # Configure a simulation stop time time and execute the simulation run
    TheDynSim.ConfigureStopTime(simulationTime)
    TheDynSim.ExecuteSimulation()

    # Pull data for post-processing and plotting
    r_BN, v_BN, sigma_BN, omega_BN_B = pull_DynOutputs(TheDynSim)
    sigma_tilde_BN, omega_tilde_BN = pull_senseOutputs(TheDynSim)
    r_sc, r_sun, r_earth, r_moon, r_mars, r_beacons = pull_DynCelestialOutputs(TheDynSim)

    ##  Post-Process sim data using camera, image processing, batch filter DINO modules



    ###############################################################################
    #
    #       Pull in canned QE and Transmission curves from DINO C-REx files
    #
    ###############################################################################

    # load tranmission curve for Canon 20D
    tc = np.load('../dinoModels/SimCode/opnavCamera/tc/20D.npz')

    # load QE curve for Hubble Space Telecope Advanced Camera for Surveys SITe CCD
    qe = np.load('../dinoModels/SimCode/opnavCamera/qe/ACS.npz')

    ###############################################################################
    #
    #       Initialize camera
    #
    ###############################################################################

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

    #need loop to define asteroids, too


    #can kill these once I change the way camera is initialized
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

    #this is spoofing the output of the nav exec
    #telling the camera when to take an image.
    takeImage = np.zeros(len(r_sc))
    takeImage[100] = 1
    takeImage[200] = 1
    takeImage[300] = 1
    takeImage[400] = 1

    lastTakeImage = 0
    for i in range(0,len(r_sc)):
        cam.scState = r_sc[i][1:4]
        earth.state = r_earth[i][1:4]
        moon.state = r_moon[i][1:4]
        mars.state = r_mars[i][1:4]
        #also need a loop here for 
        #updating beacon position once they're added
        cam.scDCM = rbk.MRP2C(sigma_BN[i][1:4])
        cam.takeImage = takeImage[i]
        cam.imgTime = r_sc[i][0]
        cam.updateState()

    import pdb 
    pdb.set_trace()
    for i in range(0,len(cam.images)):
        plt.figure()
        plt.imshow(cam.images[i].detectorArray)
    
    plt.show()

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
    #pull_DynCelestialOutputs(TheDynSim)
    plt.show()

def defineParams(inputs):
    ##  Camera Setup

    ##  Image Processing Setup

    ##  Batch Extras Setup


    return cameraObj, cameraParams, batchExtras