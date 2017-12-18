
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

def propAndObs_Scenario():
    '''

    :return:
    '''

    propSim = DINO_main.DINO_DynSim(10000., 100.)
    obsSim = DINO_main.DINO_DynSim(0.01, 0.01)

    propSim.InitializeSimulation()
    obsSim.InitializeSimulation()

    ##   Define Mode Sequence Parameters.

    propDurations = [139643.532, 10*139643.532]
    obsDurations = [1000.0]

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

    sigma_BN = np.array([[0, 0.1, 0.2, 0.3]])
    omega_BN_B = np.array([[0, 0.001, 0.002, 0.003]])

    r_sun = np.zeros([1,4])
    r_earth = np.zeros([1,4])
    r_moon = np.zeros([1,4])
    r_mars = np.zeros([1,4])
    r_beacons = np.zeros([obsSim.DynClass.numBeacons,4])


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
            sigma_tilde_BN_temp, omega_tilde_BN_temp = pull_senseOutputs(propSim,plots=False)
            sigma_hat_BN_temp, omega_hat_BN_temp = pull_aekfOutputs(propSim,plots=False)
            #Lr = pull_FSWOutputs(propSim)

            r_sun_temp, r_earth_temp, r_moon_temp, r_mars_temp, r_beacons_temp = pull_DynCelestialOutputs(propSim, plots=False)

            propInd = propInd+1

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
    propAndObs_Scenario()
