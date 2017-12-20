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

def genCamImages(cam, beaconList, r_BN, sigma_BN, omega_BN, r_earth, r_moon, r_mars, takeImage):
    """
    Generates camera images given spacecraft and celestial body positions. Intended to return all necessary/relevant
    information for genImgProc().
    :params: cam                : Camera object
    :params: beaconList         : List of nav beacon objects ordered from Earth
    :params: r_BN      : spacecraft inertial position vector over observation window
    :params: sigma_BN : spacecraft body-inertial attitude over observation window
    :params: omega_BN : spacecraft body-inertial angular rate over observation window
    :params: r_earth           : earth inertial position vector over observation window
    :params: r_moon                 : moon inertial position vector over observation window
    :params: r_mars                 : mars inertial position vector over observation window
    :params: takeImage     : vector of camera control arguments; 1 is take image, 0 is don't

    :return: detectorArrays          : list of image objects
    :return: imgTimes          : list of image times
    :return: imgPos          : list of sc positions during images (passed-through)
    :return: imgMRP          :   sc attitude MRP (passed-through)
    :return: imgBeaconPos    :   beacon inertial positions (passed-through)
    """

    earth = beaconList[0]
    mars = beaconList[1]
    moon = beaconList[2]

    lastTakeImage = 0

    #   Iterate over sampled spacecraft positions
    for i in range(0, len(r_BN)):
        #   Convert from meters to km
        cam.scState = r_BN[i][1:4] / 1000
        earth.state = r_earth[i][1:4] / 1000
        moon.state = r_moon[i][1:4] / 1000
        mars.state = r_mars[i][1:4] / 1000

        #   Set internal dcm to something random for now I guess
        cam.scDCM = np.array([
            [0, 1, 0],
            [-1, 0, 0],
            [0, 0, 1]
        ])

        if i == 100 or i == 101:
            sc2bdy = earth.state - cam.scState
        elif i == 200 or i == 201:
            sc2bdy = moon.state - cam.scState
        else:
            sc2bdy = mars.state - cam.scState

        #sc2bdyNormed = sc2bdy / np.linalg.norm(sc2bdy)
        #RA = np.arctan2(sc2bdyNormed[1], sc2bdyNormed[0])
        #DE = np.arctan2(sc2bdyNormed[2], np.sqrt(sc2bdyNormed[0] ** 2 + sc2bdyNormed[1] ** 2))
        #cam.scDCM = rbk.euler3212C(np.array([RA, -DE, 0]))
        cam.scDCM = rbk.MRP2C(sigma_BN[i,1:4])

        if np.linalg.norm(omega_BN[i,1:4]) < math.degrees(0.1):
            cam.takeImage = 1
            cam.imgTime = r_BN[i][0]
        else:
            cam.takeImage = 0

        cam.updateState()

    detectorArrays = []
    imgTimes = []
    imgPos = []
    imgMRP = []
    imgBeaconPos = []

    for i in range(0, len(cam.images)):
        detectorArrays.append(cam.images[i].detectorArray)
        imgTimes.append(cam.images[i].imgTime)
        imgPos.append(cam.images[i].imgPos)
        imgMRP.append(rbk.C2MRP(cam.images[i].imgDCM))
        imgBeaconPos.append(cam.images[i].imgBeaconPos)

        plt.figure()
        plt.imshow(cam.images[i].detectorArray)

    return detectorArrays, imgTimes, imgPos, imgMRP, imgBeaconPos


def genImgProc(detectorArrays, imgPos, imgMRP, imgBeaconPos, imgTimes, ipParam):
    # required parameters from defineParameters function
    camParamIP = ipParam[0]
    beaconIDs = ipParam[1]
    beaconRadius = ipParam[2]

    imgTimesFound = []
    beaconIDsFound = []
    beaconPLFound = []
    imgMRPFound = []  # will have 'None' entries when not able to detect enough objects
    imgMRPFoundPassThrough = []  # identical attitude as input into the image processing module

    print '\nBeacon ID Check:'
    print beaconIDs

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

    # Generate inputs for navigation modulec
    numNavInputs = len(imgTimesFound)
    imgTimesNav = np.reshape(imgTimesFound, (numNavInputs, 1))
    beaconIDsNav = np.reshape(beaconIDsFound, (numNavInputs, 1))
    beaconPLNav = np.reshape(beaconPLFound, (numNavInputs, 2))
    imgMRPNav = np.reshape(imgMRPFoundPassThrough, (numNavInputs, 3))

    return beaconPLNav, beaconIDsNav, imgMRPNav, imgTimesNav, numNavInputs


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

    # init values for camera that will be set later.
    scState = -1
    scDCM = -1
    takeImage = 0

    cam = camera.camera(
        camSensorSize[0],  # detector width in m
        camSensorSize[1],  # detector height in m
        camFocalLength,  # focal lenght in m
        camResolution[0],  # detector resolution (width direction)
        camResolution[1],  # detector resolution (height direction)
        np.identity(3),  # body2cameraDCM
        1000,  # maximum magnitude (for debugging)
        -1000,  # minimum magnitude (for debugging)
        qe,  # quantum efficiency dictionary
        tc,  # transmission curve dictionary
        lambdaBinSize,  # lambda bin size
        effectiveArea,  # effective area in m^2
        darkCurrent,  # dark current in electrons per second
        readSTD,  # std for read noise in electrons
        binSize,  # bin size
        maxBinDepth,  # max bin depth
        psfSTD,  # std for psf
        simTimeStep,  # simulation timestep
        scState,  # position state of s/c
        scDCM,  # intertal 2 body DCM for s/c
        beacons,  # bodies to track in images
        takeImage,  # takeImage message
        debug={
            'addStars': 0, 'rmOcc': 1, 'addBod': 1, 'psf': 1,
            'raster': 1, 'photon': 1, 'dark': 1, 'read': 1,
            'verbose': 1, 'hotDark': 1},
        db='../dinoModels/SimCode/opnavCamera/db/tycho.db'  # stellar database
    )

    # Camera Module Parameter Creation

    camInputs = cam

    # Image Processing Module Parameter Creation
    ipCamParam = {}
    ipCamParam['resolution'] = (cam.resolutionWidth, cam.resolutionHeight)
    ipCamParam['focal length'] = cam.focalLength
    ipCamParam['sensor size'] = (cam.detectorWidth, cam.detectorHeight)
    ipCamParam['pixel size'] = (
        cam.detectorWidth / cam.resolutionWidth,
        cam.detectorHeight / cam.resolutionHeight)
    ipCamParam['field of view'] = (
        cam.angularWidth,
        cam.angularHeight)
    ipInputs = [ipCamParam, beaconIDs, beaconRadius]

    # Nav Module Parameter Creation

    navInputs = {}

    # SPICE Parameters

    # basic .bsp filename (generic, such as de430, etc)
    navInputs['basic_bsp'] = 'de430.bsp'
    # .bsp filename for mission
    navInputs['mission_bsp'] = 'DINO_kernel.bsp'
    # .tls filename
    navInputs['tls'] = 'naif0011.tls'
    # abcorr for spkzer
    navInputs['abcorr'] = 'NONE'
    # reference frame
    navInputs['ref_frame'] = 'J2000'

    # Force Parameters

    #   Gravity
    # body vector for primary and secondary gravitational bodies
    navInputs['bodies'] = ['SUN', '3', '399']
    # specify primary and secondary indices
    navInputs['primary'] = 0
    navInputs['secondary'] = [1, 2]
    # respective GP vector
    navInputs['mu'] = [1.32712428 * 10 ** 11, 3.986004415 * 10 ** 5, 4.305 * 10 ** 4]
    #   SRP
    # A/M ratio multiplied by solar pressure constant at 1 AU with adjustments
    # Turboprop document Eq (64)
    navInputs['SRP'] = 0.3 ** 2 / 14. * 149597870. ** 2 * 1358. / 299792458. / 1000.
    # coefficient of reflectivity
    navInputs['cR'] = 1.

    # Camera/P&L Parameters

    # Focal Length (mm)
    navInputs['FoL'] = ipCamParam['focal length']
    # default inertial to camera transformation matrices
    navInputs['DCM_BI'] = np.eye(3)
    navInputs['DCM_TVB'] = np.eye(3)
    # Camera resolution (pixels)
    navInputs['resolution'] = [cam.resolutionWidth, cam.resolutionHeight]
    # width and height of pixels in camera. convert from meters to mm
    navInputs['pixel_width'] = ipCamParam['pixel size'][0] * 10 ** 3
    navInputs['pixel_height'] = ipCamParam['pixel size'][1] * 10 ** 3
    # direction coefficient of pixel and line axes
    navInputs['pixel_direction'] = 1.
    navInputs['line_direction'] = 1.

    # Add anomaly detection parameters
    navInputs['anomaly'] = False
    navInputs['anomaly_num'] = 0
    navInputs['anomaly_threshold'] = 4

    # plotting? 'ON' or 'OFF'
    navInputs['nav plots'] = 'ON'

    return camInputs, ipInputs, navInputs


def genNavOutputs(beaconPLNav, beaconIDsNav, imgTimesNav, imgMRPNav, r_BN, v_BN, navParam):
    # Run the Navigation Module
    observationData = {}
    stateValues = {}

    # pull the initial state
    IC = np.append(r_BN[0, 1:4], v_BN[0, 1:4]) / 1000.

    # estimating accelerations - ON or OFF
    navParam['acc_est'] = 'OFF'

    # number of batch iterations per data package
    navParam['iterations'] = 3

    # a priori uncertainty for the referenceStates
    covBar = np.zeros((9, 9))
    covBar[0, 0] = 10000. ** 2
    covBar[1, 1] = 10000. ** 2
    covBar[2, 2] = 10000. ** 2
    covBar[3, 3] = .1 ** 2
    covBar[4, 4] = .1 ** 2
    covBar[5, 5] = .1 ** 2
    covBar[6, 6] = (10. ** (-8)) ** 2
    covBar[7, 7] = (10. ** (-8)) ** 2
    covBar[8, 8] = (10. ** (-8)) ** 2

    if navParam['acc_est'] == 'OFF':
        covBar = covBar[0:6, 0:6]
    else:
        IC = np.append(np.zeros((3,)))

        # Inverse of the observation weighting matrix (W)
    observationUncertainty = np.identity(2)
    observationUncertainty[0, 0] = 2 ** 2
    observationUncertainty[1, 1] = 2 ** 2

    # the initial STM is an identity matrix
    phi0 = np.identity(IC.shape[0])

    # initiate a priori deviation
    stateDevBar = np.zeros(IC.shape)

    observationData['measurements'] = beaconPLNav
    observationData['beaconIDs'] = np.squeeze(beaconIDsNav).tolist()
    for nn, ii in enumerate(observationData['beaconIDs']):
        # if ii == 'Earth':
        #   observationData['beaconIDs'][nn] = '399'
        # if ii == 'Moon':
        #   observationData['beaconIDs'][nn] = '301'
        if ii == 'Mars':
            observationData['beaconIDs'][nn] = '4'

    observationData['observation uncertainty'] = observationUncertainty

    stateValues['IC'] = IC
    stateValues['phi0'] = phi0
    stateValues['covBar'] = covBar
    stateValues['stateDevBar'] = stateDevBar
    stateValues['initial time'] = r_BN[0, 0]

    # convert time from nanoseconds to seconds
    obsTimes = np.squeeze(imgTimesNav / 10 ** 9)

    # convert imgMRPNav from MRP to 3-2-1 euler
    img321Nav = np.zeros(imgMRPNav.shape)
    for ii in xrange(imgMRPNav.shape[0]):
        img321Nav[ii, :] = rbk.MRP2Euler321(imgMRPNav[ii, :])

    filterOutputs = initBatchFnc(stateValues, obsTimes, observationData, \
                                 imgMRPNav, navParam)
    # output of the estimated state for the input times after number of desired iterations
    # (N,d) in shape. Numpy array of floats
    estimatedState = filterOutputs[str(navParam['iterations'] - 1)]['estimatedState']
    return filterOutputs
