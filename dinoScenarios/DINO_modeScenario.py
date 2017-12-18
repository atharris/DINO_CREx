
##  Standard Python library imports
import sys, os, inspect
import matplotlib.pyplot as plt
from numpy import linalg as la
import numpy as np
import math
import datetime

##  Set path for BSK/DINO imports; assumes BSK path is in same root folder as DINO
bskName = 'Basilisk'
bskPath = '../..' + '/' + bskName + '/'
sys.path.append(bskPath + 'modules')
sys.path.append(bskPath + 'PythonModules')
sys.path.append('../DINObatch/P&L/common_functions/')
sys.path.append('../dinoModels/SimCode/opnavCamera/')
sys.path.append('../dinoModels/SimCode/opnavCamera/dependencies')
sys.path.append('../dinoModels/fswAlgorithms/imageProcessing/dependencies/')
sys.path.append('../DINObatch/pixelAndLine/commonFunctions/')
sys.path.append('../DINObatch/pixelAndLine/unmodeledAcc/')
sys.path.append('../DINObatch/SPICE/')

##  BSK imports
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

##  Import DINO modules
import camera
import imageProcessingExecutive as ip
from initBatch import initBatchFnc

##  Import relevant DINO scenario support modules
import DINO_main
import BSK_plotting as BSKPlt
from DINO_multiScenarios import *
from DINO_logAndPlot import *
from DINO_scenarioUtils import *


def runSimSegment(TheDynSim, simTime, initialState, initialAtt, timeStr):
    '''

    :param TheDynSim:
    :param simTime:
    :param initialState:
    :param initialAtt:
    :param timeStr:
    :return:
    '''
    simTime = mc.sec2nano(simTime)
    samplingTime = mc.sec2nano(min([TheDynSim.fswUpdateRate, TheDynSim.dynUpdateRate]))
    numDataPoints = simTime / samplingTime
    log_DynOutputs(TheDynSim, samplingTime)

    #   Ensure that the observation mode starts with the right sc positions, time
    posRef = TheDynSim.DynClass.scObject.dynManager.getStateObject("hubPosition")
    velRef = TheDynSim.DynClass.scObject.dynManager.getStateObject("hubVelocity")
    sigmaRef = TheDynSim.DynClass.scObject.dynManager.getStateObject("hubSigma")
    omegaRef = TheDynSim.DynClass.scObject.dynManager.getStateObject("hubOmega")
    posRef.setState(sp.np2EigenVectorXd(initialState[0:3]))  # r_BN_N [m]
    velRef.setState(sp.np2EigenVectorXd(initialState[3:6]))  # v_BN_N [m/s]
    sigmaRef.setState(sp.np2EigenVectorXd(initialAtt[0:3]))  # sigma_BN_B
    omegaRef.setState(sp.np2EigenVectorXd(initialAtt[3:6]))  # omega_BN_B [rad/s]

    TheDynSim.DynClass.spiceObject.UTCCalInit = timeStr

    TheDynSim.ConfigureStopTime(simTime)
    TheDynSim.ExecuteSimulation()



    return TheDynSim

def propAndObs_Scenario(useNavOutputs, genPlots):
    """
    Executes a default scenario for stand-alone dynamic simulations
    :params: None
    :return: None
    """

    propSim = DINO_main.DINO_DynSim(100., 100.)
    obsSim = DINO_main.DINO_DynSim(0.01, 0.01)

    propSim.InitializeSimulation()
    obsSim.InitializeSimulation()

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

    ##   Define Mode Sequence Parameters.

    propDurations = [1000, 10*139643.532]
    obsDurations = [100.0]

    modeSeq = [0, 1, 0]

    propInd = 0
    obsInd = 0

    ##  Define mode setup parameters.
    startTime = "11 Jul 2020 00:00:37.034" ##   Time given as UTC coordinated time for SPICE
    utcObj= datetime.datetime.strptime(startTime, "%d %b %Y %H:%M:%S.%f") ## datetime object used to store, update time.

    ##  Define initial trajectory, maneuver parameters.
    # Set up the orbit using classical orbit elements
 #   oe = define_default_orbit()
    mu = propSim.DynClass.earthGravBody.mu
    r_BN = np.zeros([1,4])
    v_BN = np.zeros([1,4])
    r_BN[0,1:4], v_BN[0,1:4] = define_dino_earthSOI()#om.elem2rv(mu, oe)

    r_BN_batch = np.zeros([1,4])
    v_BN_batch = np.zeros([1,4])
    r_BN_batch[0,1:4], v_BN_batch[0,1:4] = define_dino_earthSOI()

    sigma_BN = np.array([[0, 0.1, 0.2, 0.3]])
    omega_BN_B = np.array([[0, 0.001, 0.002, 0.003]])

    r_sun = np.zeros([1,4])
    r_earth = np.zeros([1,4])
    r_moon = np.zeros([1,4])
    r_mars = np.zeros([1,4])

    for mode in modeSeq:
        timeStr = utcObj.strftime("%d %b %Y %H:%M:%S.%f")

        if mode == 0:
            ##  Run the propagation sim
            print "*****************************"
            print "* Propagation Mode Starting *"
            print "*****************************"
            execTime = propDurations[propInd]
            runSimSegment(propSim, propDurations[propInd], np.hstack([r_BN[-1, 1:4], v_BN[-1, 1:4]]).tolist(), np.hstack([sigma_BN[-1, 1:4],omega_BN_B[-1, 1:4]]).tolist() , timeStr)
            r_BN_temp, v_BN_temp, sigma_BN_temp, omega_BN_B_temp = pull_DynOutputs(propSim,plots=False)
            #sigma_tilde_BN_temp, omega_tilde_BN_temp = pull_senseOutputs(propSim,plots=False)
            #sigma_hat_BN_temp, omega_hat_BN_temp = pull_aekfOutputs(propSim,plots=False)

            r_sun_temp, r_earth_temp, r_moon_temp, r_mars_temp, r_beacons_temp = pull_DynCelestialOutputs(propSim, plots=False)

        else:
            ##  Run the observation sim
            print "*****************************"
            print "* Observation Mode Starting *"
            print "*****************************"
            execTime = obsDurations[obsInd]
            runSimSegment(obsSim, obsDurations[obsInd], np.hstack([r_BN[-1, 1:4], v_BN[-1, 1:4]]).tolist(), np.hstack([sigma_BN[-1, 1:4],omega_BN_B[-1, 1:4]]).tolist(), timeStr)
            r_BN_temp, v_BN_temp, sigma_BN_temp, omega_BN_B_temp = pull_DynOutputs(obsSim,plots=False)
            sigma_tilde_BN_temp, omega_tilde_BN_temp = pull_senseOutputs(obsSim,plots=False)
            sigma_hat_BN_temp, omega_hat_BN_temp = pull_aekfOutputs(obsSim,plots=False)
            #Lr = pull_FSWOutputs(obsSim)
            r_sun_temp, r_earth_temp, r_moon_temp, r_mars_temp, r_beacons_temp = pull_DynCelestialOutputs(obsSim,plots=False)
            takeImage = np.ones([1,len(r_BN_temp)])

            batchIC = np.hstack([r_BN_batch[-1,1:],v_BN_batch[-1,1:]])
            phi0 = np.eye(6)
            timeSpan = r_BN_temp[0,:] * mc.NANO2SEC
            propagatorInput = (batchIC, phi0, timeSpan, navParam)

            outState = bf.runRef(propagatorInput)
            r_BN_batch_temp = outState[0:3,:]
            v_BN_batch_temp = outState[3:,:]

            if useNavOutputs:
                camInputPos = r_BN_batch_temp
            else:
                camInputPos = r_BN_temp

            detectorArrays, imgTimes, imgPos, imgMRP, imgBeaconPos = genCamImages(cam, beacons, camInputPos, r_earth, r_moon,
                                                                                  r_mars, takeImage)
            # Run the Image Processing Module


            beaconPLNav, beaconIDsNav, imgMRPNav, imgTimesNav, numNavInputs = genImgProc(detectorArrays, imgPos, imgMRP,
                                                                                         imgBeaconPos, imgTimes,
                                                                                         ipParam)

            print "*******Image Processing Outputs*******"
            print "Beacon PL:"
            print beaconPLNav
            print "beaconIDs:"
            print beaconIDsNav
            propInd = propInd + 1
            obsInd = obsInd+1



        utcObj = utcObj + datetime.timedelta(seconds=execTime)


        #   Update celestial positions
        r_sun = np.concatenate((r_sun, r_sun_temp))
        r_earth = np.concatenate((r_earth, r_earth_temp))
        r_moon = np.concatenate((r_moon, r_moon_temp))
        r_mars = np.concatenate((r_mars, r_mars_temp))
        #r_beacons_temp = np.concatenate((r_beacons, r_beacons_temp))

        #   Update SC dynamic params
        r_BN = np.concatenate((r_BN, r_BN_temp))
        v_BN = np.concatenate((v_BN, v_BN_temp))
        sigma_BN = np.concatenate((sigma_BN, sigma_BN_temp))
        omega_BN_B = np.concatenate((omega_BN_B, omega_BN_B_temp))

        #   Update fsw dynaic params
        r_BN_batch = np.concatenate((r_BN_batch, r_BN_batch_temp))
        v_BN_batch = np.concatenate((v_BN_batch, v_BN_batch_temp))

    ##  Plotting of results
    celes_data_dict = {
        'moon': [r_moon, 'cyan'],
        'earth': [r_earth, 'dodgerblue'],
        'mars': [r_mars, 'r'],
        'sun': [r_sun, 'orange']
    }

    BSKPlt.plot_spacecraft_orbit(celes_data_dict, r_BN)
    plt.show()


if __name__ == "__main__":
    propAndObs_Scenario(False, False)
