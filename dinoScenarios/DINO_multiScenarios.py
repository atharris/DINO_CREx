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

import macros as mc
import unitTestSupport as sp
import orbitalMotion as om
import BSK_plotting as BSKPlt

# ------------------------------------- DATA LOGGING ------------------------------------------------------ #

def log_DynCelestialOutputs(TheDynSim, samplingTime):
    TheDynSim.TotalSim.logThisMessage(TheDynSim.DynClass.scObject.scStateOutMsgName, samplingTime)
    TheDynSim.TotalSim.logThisMessage(TheDynSim.DynClass.earthGravBody.bodyInMsgName, samplingTime)
    TheDynSim.TotalSim.logThisMessage(TheDynSim.DynClass.marsGravBody.bodyInMsgName, samplingTime)
    TheDynSim.TotalSim.logThisMessage(TheDynSim.DynClass.sunGravBody.bodyInMsgName, samplingTime)
    TheDynSim.TotalSim.logThisMessage(TheDynSim.DynClass.moonGravBody.bodyInMsgName, samplingTime)
    for ind in range(0,len(TheDynSim.DynClass.beaconList)):
        TheDynSim.TotalSim.logThisMessage(TheDynSim.DynClass.beaconList[ind].scStateOutMsgName, samplingTime)

    return

def log_DynOutputs(TheBSKSim, samplingTime):
    TheBSKSim.TotalSim.logThisMessage(TheBSKSim.DynClass.scObject.scStateOutMsgName, samplingTime)
    TheBSKSim.TotalSim.logThisMessage(TheBSKSim.DynClass.simpleNavObject.outputAttName, samplingTime)
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
    for ind in range(0,len(TheDynSim.DynClass.beaconList)):
        r_beacons.append(TheDynSim.pullMessageLogData(TheDynSim.DynClass.beaconList[ind].scStateOutMsgName + '.r_BN_N', range(3)))

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
        'r_beacon_0': [r_beacons[0], 'g'],
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
        'r_beacon_0': [r_beacons[0], 'g'],
        # 'r_beacon_1': [r_beacons[1], 'g'],
        # 'r_beacon_2': [r_beacons[2], 'g'],
        # 'r_beacon_3': [r_beacons[3], 'g'],
        # 'r_beacon_4': [r_beacons[4], 'g'],
        # 'r_beacon_5': [r_beacons[5], 'g'],
        # 'r_beacon_6': [r_beacons[6], 'g'],
        # 'r_beacon_7': [r_beacons[7], 'g'],
        # 'r_beacon_8': [r_beacons[8], 'g']
    }
    BSKPlt.plot_spacecraft_orbit(sc_dict_data_color, r_sc)



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
    simulationTime = mc.sec2nano(139643.532)
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
    pull_DynOutputs(TheDynSim)
    pull_DynCelestialOutputs(TheDynSim)
    plt.show()