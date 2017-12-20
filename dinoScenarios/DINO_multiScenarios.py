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
from DINO_logAndPlot import *
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


# ------------------------------------- SCENARIOS ------------------------------------------------------ #

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
    simulationTime = mc.sec2nano(10000.)
    numDataPoints = int(orbPeriod)
    samplingTime = simulationTime / (numDataPoints - 1)
    log_DynOutputs(TheDynSim, samplingTime)
    log_FSWOutputs(TheDynSim, samplingTime)

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
    pull_FSWOutputs(TheDynSim)
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
    # earth.albedo = 1e30

    mars = camera.beacon()
    mars.r_eq = 3396.2
    mars.id = 'Mars'
    mars.albedo = 0.17

    moon = camera.beacon()
    moon.r_eq = 1738.1
    moon.id = 'Moon'
    # moon.albedo = 0.12
    moon.albedo = .7

    beacons = [earth, mars, moon]
    # need loop to define asteroids, too

    cam, ipParam, navParam = defineParameters(
        (512, 512),  # camera resolution, width then height
        0.05,  # focal length in m
        (0.01, 0.01),  # detector dimensions in m, with then height
        beacons,  # list of beacons
        # transmission curve dict
        np.load('../dinoModels/SimCode/opnavCamera/tc/20D.npz'),
        # quantum efficiency curve dict
        np.load('../dinoModels/SimCode/opnavCamera/qe/ACS.npz'),
        1,  # bin size for wavelength functions (in nm)
        0.01 ** 2,  # effective area (m^2)
        100,  # dark current electrons/s/pixel
        100,  # read noise STD (in electrons per pixel)
        100,  # bin size
        2 ** 32,  # saturation depth
        1,  # Standard deviation for PSF (in Pixels)
        0.01  # simulation timestep
    )

    # this is spoofing the output of the nav exec
    # telling the camera when to take an image.
    takeImage = np.zeros(len(r_BN))
    takeImage[100] = 1
    takeImage[200] = 1
    takeImage[300] = 1

            #cam, beaconList, r_BN, r_earth, r_moon, r_mars

    detectorArrays, imgTimes, imgPos, imgMRP, imgBeaconPos = genCamImages(cam, beacons, r_BN, r_earth, r_moon, r_mars, takeImage)
    # Run the Image Processing Module


    beaconPLNav, beaconIDsNav, imgMRPNav, imgTimesNav, numNavInputs = genImgProc(detectorArrays, imgPos, imgMRP, imgBeaconPos, imgTimes, ipParam)

    print "*******Image Processing Outputs*******"
    print "Beacon PL:"
    print beaconPLNav
    print "beaconIDs:"
    print beaconIDsNav

    filterOutputs = genNavOutputs(beaconPLNav, beaconIDsNav, imgTimesNav, imgMRPNav, r_BN, v_BN, navParam)
    print "********** Filter Run Complete ************"

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
    pull_aekfOutputs(TheDynSim)
    # pull_DynCelestialOutputs(TheDynSim)
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