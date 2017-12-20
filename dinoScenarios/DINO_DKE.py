import sys, os, inspect
import numpy as np

bskName = 'Basilisk'
bskPath = '../..' + '/' + bskName + '/'
sys.path.append(bskPath + 'modules')
sys.path.append(bskPath + 'PythonModules')
bskSpicePath = bskPath + 'External/EphemerisData/'
sys.path.append('../dinoModels/SimCode/opnavCamera/')

try:
    import macros as mc
    import unitTestSupport as sp
    import sim_model
    import spacecraftPlus
    import gravityEffector
    import simple_nav
    import spice_interface
    import ephemeris_converter
    import radiation_pressure
    import star_tracker
    import imu_sensor

    import simIncludeRW
    import reactionWheelStateEffector
    import rwVoltageInterface
    import simIncludeGravBody
    import ExtForceTorque

    # import message declarations
    import fswMessages
except ImportError:
    from Basilisk import __path__
    import Basilisk.utilities.macros as mc
    import Basilisk.utilities.unitTestSupport as sp
    from Basilisk.utilities import simIncludeRW
    from Basilisk.simulation import sim_model, spacecraftPlus, gravityEffector, simple_nav, spice_interface
    from Basilisk.simulation import ephemeris_converter, radiation_pressure, star_tracker, imu_sensor
    from Basilisk.simulation import reactionWheelStateEffector, rwVoltageInterface

    bskSpicePath = __path__[0] + '/supportData/EphemerisData/'

#import simMessages

import opnavCamera


#   Define the base class for simulation dynamics
class DynamicsClass():
    #   Constructor method; sets basic parameters for the simulation
    def __init__(self, SimBase, updateRate=0.01):
        # Define process name, task name and task time-step
        self.processName = SimBase.DynamicsProcessName
        self.taskName = "DynamicsTask"
        self.taskTimeStep = mc.sec2nano(0.01)
        self.taskTimeStep = mc.sec2nano(updateRate)

        # Create tasks
        SimBase.dynProc.addTask(SimBase.CreateNewTask(self.taskName, self.taskTimeStep))
        #SimBase.dynPyProc.createPythonTask("opnavCameraTask", self.taskTimeStep,True, 30)

        # Instantiate Dyn modules as objects
        self.scObject = spacecraftPlus.SpacecraftPlus()
        self.spiceObject = spice_interface.SpiceInterface()

        self.earthGravBody = gravityEffector.GravBodyData()
        self.moonGravBody = gravityEffector.GravBodyData()
        self.marsGravBody = gravityEffector.GravBodyData()
        self.sunGravBody = gravityEffector.GravBodyData()

        self.simpleNavObject = simple_nav.SimpleNav()
        self.ephemConvert = ephemeris_converter.EphemerisConverter()

        self.extForceTorque = ExtForceTorque.ExtForceTorque()
        self.extForceTorque.ModelTag = "externalForceTorque"
        self.scObject.addDynamicEffector(self.extForceTorque)

        # Lists to store gravity bodies and spice planet data
        self.gravBodyList = []
        self.spicePlanetNames = []

        self.rwFactory = simIncludeRW.rwFactory()
        self.varRWModel = self.rwFactory.BalancedWheels


        # Initialize all modules and write init one-time messages
        self.InitAllDynObjects()


        # Assign initialized modules to tasks
        SimBase.AddModelToTask(self.taskName, self.spiceObject, None, 20)
        SimBase.AddModelToTask(self.taskName, self.ephemConvert, None, 19)
        SimBase.AddModelToTask(self.taskName, self.scObject, None, 10)
        SimBase.AddModelToTask(self.taskName, self.simpleNavObject, None, 9)
        SimBase.AddModelToTask(self.taskName, self.srpDynEffector, None, 18)
        SimBase.AddModelToTask(self.taskName, self.starTracker,None, 8)
        SimBase.AddModelToTask(self.taskName, self.gyroModel,None,7)
        SimBase.AddModelToTask(self.taskName, self.rwStateEffector,None,6)


    # ------------------------------------------------------------------------------------------- #
    # These are module-initialization methods
    def SetSpacecraftObject(self):
        self.scObject.ModelTag = "spacecraftBody"
        # -- Crate a new variable for the sim sc inertia I_sc. Note: this is currently accessed from FSWClass
        self.I_sc = [9., 0., 0.,
                    0., 8., 0.,
                    0., 0., 6.]
        self.scObject.hub.mHub = 10.0  # kg - spacecraft mass
        self.scObject.hub.r_BcB_B = [[0.0], [0.0], [0.0]]  # m - position vector of body-fixed point B relative to CM
        self.scObject.hub.IHubPntBc_B = sp.np2EigenMatrix3d(self.I_sc)
        self.scObject.hub.useTranslation = True
        self.scObject.hub.useRotation = True

    def SetSimpleNavObject(self):
        self.simpleNavObject.ModelTag = "SimpleNavigation"

    def SetGravityBodies(self):

        # clear prior gravitational body and SPICE setup definitions
        self.gravFactory = simIncludeGravBody.gravBodyFactory()

        # setup Earth Gravity Body
        earth = self.gravFactory.createEarth()
        earth.isCentralBody = True  # ensure this is the central gravitational body
        mu = earth.mu
        def AddMoon(self):
            self.moonGravBody.bodyInMsgName = "moon_planet_data"
            self.moonGravBody.outputMsgName = "moon_display_frame_data"
            self.moonGravBody.mu = 4.902799E12  # meters^3/s^2
            self.moonGravBody.radEquator = 1738100.0  # meters
            self.moonGravBody.isCentralBody = False
            self.moonGravBody.useSphericalHarmParams = True
            # Store Moon celestial body in the Gravity and Spice lists
            self.gravBodyList.append(self.moonGravBody)
            self.spicePlanetNames.append(self.moonGravBody.bodyInMsgName[:-12])
        def AddEarth(self):
            self.earthGravBody = self.gravFactory.createEarth()
            self.earthGravBody.isCentralBody = False  # ensure this is the central gravitational body
            #self.earthGravBody.bodyInMsgName = "earth_planet_data"
            #self.earthGravBody.outputMsgName = "earth_display_frame_data"
            #self.earthGravBody.mu = 0.3986004415E+15  # meters^3/s^2
            #self.earthGravBody.radEquator = 6378136.6  # meters
            #self.earthGravBody.isCentralBody = False
            #self.earthGravBody.useSphericalHarmParams = True
            # Store Earth celestial body in the Gravity and Spice lists
            self.gravBodyList.append(self.earthGravBody)
            self.spicePlanetNames.append(self.earthGravBody.bodyInMsgName[:-12])
        def AddMars(self):
            self.marsGravBody.bodyInMsgName = "mars barycenter_planet_data"
            self.marsGravBody.outputMsgName = "mars_barycenter_display_frame_data"
            self.marsGravBody.mu = 4.28283100e13  # meters^3/s^2
            self.marsGravBody.radEquator = 3396190  # meters
            self.marsGravBody.isCentralBody = False
            self.marsGravBody.useSphericalHarmParams = True
            # Store Mars celestial body in the Gravity and Spice lists
            self.gravBodyList.append(self.marsGravBody)
            self.spicePlanetNames.append(self.marsGravBody.bodyInMsgName[:-12])
        def AddSun(self):
            self.sunGravBody.bodyInMsgName = "sun_planet_data"
            self.sunGravBody.outputMsgName = "sun_display_frame_data"
            self.sunGravBody.mu = 1.32712440018E20 # meters^3/s^2
            self.sunGravBody.radEquator = 695508000.0 # meters
            self.sunGravBody.isCentralBody = True
            self.sunGravBody.useSphericalHarmParams = True
            # Store Sun celestial body in the Gravity and Spice lists
            self.gravBodyList.append(self.sunGravBody)
            self.spicePlanetNames.append(self.sunGravBody.bodyInMsgName[:-12])

        AddEarth(self)
        AddMars(self)
        AddSun(self)
        AddMoon(self)
        # Attach gravity model to spaceCraftPlus
        self.scObject.gravField.gravBodies = spacecraftPlus.GravBodyVector(self.gravBodyList)
        print '\n' + 'GRAVITY DATA'
        print 'Gravity cel bodies: ',
        for body in self.gravBodyList:
            print 'body: ', body.bodyInMsgName, 'central: ', body.isCentralBody
            if body.isCentralBody:
                self.mu = body.mu


    def SetSpiceObject(self):
        self.spiceObject.UTCCalInit = "11 Jul 2020 00:00:37.034"
        self.spiceObject.ModelTag = "SpiceInterfaceData"
        self.spiceObject.SPICEDataPath = bskSpicePath
        self.spiceObject.outputBufferCount = 10000
        self.spiceObject.planetNames = spice_interface.StringVector(self.spicePlanetNames)
        print '\n' + 'SPICE DATA'
        print 'Spice cel bodies: ', self.spicePlanetNames
        # By default the SPICE object will use the solar system barycenter as the inertial origin
        # If the spacecraftPlus() output is desired relative to another celestial object, the zeroBase string
        # name of the SPICE object needs to be changed.
        self.spiceObject.zeroBase = 'sun' #'mars barycenter' #'sun'#'mars barycenter'
        print 'Spice zero base: ',  self.spiceObject.zeroBase

    def SetEphemerisConverter(self):
        self.sunConvertName = 'sun_ephemeris_data'
        self.earthConvertName = 'earth_ephemeris_data'
        self.marsConvertName = 'mars_ephemeris_data'
        self.moonConvertName = 'moon_ephemeris_data'

        self.ephemConvert.ModelTag = "ephemerisConverter"
        messageMap = {}
        messageMap[self.sunGravBody.bodyInMsgName] = self.sunConvertName
        messageMap[self.earthGravBody.bodyInMsgName] = self.earthConvertName
        messageMap[self.marsGravBody.bodyInMsgName] = self.marsConvertName
        messageMap[self.moonGravBody.bodyInMsgName] = self.moonConvertName
        self.ephemConvert.messageNameMap = ephemeris_converter.map_string_string(messageMap)

    def SetSRPModel(self):
        self.srpDynEffector =  radiation_pressure.RadiationPressure()
        self.srpDynEffector.ModelTag = "RadiationPressure"
        self.srpDynEffector.setUseCannonballModel(True)
        self.srpDynEffector.area = 4.
        self.srpDynEffector.coefficientReflection = 1.2
        self.srpDynEffector.sunEphmInMsgName = self.sunGravBody.outputMsgName
        self.srpDynEffector.stateInMsgName = self.scObject.scStateOutMsgName
        self.scObject.addDynamicEffector(self.srpDynEffector)


    def SetBeacons(self):
        self.beaconList =[]
        self.smaList = []
        self.eccList = []
        ephemFile = open("observations_log.txt", 'r')

        lines = ephemFile.readlines()
        for ind in range(1, 2):
            self.beaconList.append(spacecraftPlus.SpacecraftPlus())
            self.beaconList[ind-1].ModelTag = 'beaconBody_'+str(ind-1)
            self.beaconList[ind-1].hub.useTranslation = True
            self.beaconList[ind-1].hub.useRotation = False
            self.beaconList[ind-1].scStateOutMsgName = 'beacon_state_output_'+str(ind-1)
            self.beaconList[ind-1].scMassStateOutMsgName = 'beacon_mass_output_'+str(ind-1)
            self.beaconList[ind-1].gravField.gravBodies = self.scObject.gravField.gravBodies
            #print "Grav body list:", len(self.beaconList[ind - 1].gravField.gravBodies)
            self.beaconList[ind-1].hub.mHub = 750.0  # kg - spacecraft mass
            self.beaconList[ind-1].hub.r_BcB_B = [[0.0], [0.0],
                                         [0.0]]  # m - position vector of body-fixed point B relative to CM
            self.beaconList[ind-1].hub.IHubPntBc_B = sp.np2EigenMatrix3d(self.I_sc)
        ephemFile.close()

    def AddStarTracker(self):
        self.starTracker = star_tracker.StarTracker()
        self.starTracker.ModelTag = 'StarTracker'
        self.starTracker.inputStateMessage = self.scObject.scStateOutMsgName
        self.starTracker.outputStateMessage = "star_tracker_state"

        senNoiseStd = 0.01
        PMatrix = [0.0] * 3 * 3
        PMatrix[0 * 3 + 0] = PMatrix[1 * 3 + 1] = PMatrix[2 * 3 + 2] = senNoiseStd
        errorBounds = [1e-15] * 3
        self.starTracker.walkBounds = np.array(errorBounds)
        self.starTracker.PMatrix = np.array(PMatrix).reshape(3,3)

    def AddGyro(self):
        self.gyroModel = imu_sensor.ImuSensor()
        self.gyroModel.ModelTag = "GyroModel"
        self.gyroModel.InputStateMsg = self.scObject.scStateOutMsgName
        self.gyroModel.OutputDataMsg = "gyro_output_data"

        self.gyroModel.sensorPos_B = imu_sensor.DoubleVector([0.0, 0.0, 0.0])
        self.gyroModel.setBodyToPlatformDCM(0.0, 0.0, 0.0)
        try:
            self.gyroModel.accelLSB = 0.0
            self.gyroModel.gyroLSB = 0.0
        except ValueError:
            return
        self.gyroModel.senRotBias = [0.0] * 3
        self.gyroModel.senTransBias = [0.0] * 3
        self.gyroModel.senRotMax = 1.0e6
        self.gyroModel.senTransMax = 1.0e6
        
        senNoiseStd = 0.01
        PMatrixGyro = [0.0] * 3 * 3
        PMatrixGyro[0 * 3 + 0] = PMatrixGyro[1 * 3 + 1] = PMatrixGyro[2 * 3 + 2] = senNoiseStd
        errorBoundsGyro = [1e6] * 3
        
        PMatrixAccel = [0.0] * 3 * 3
        errorBoundsAccel = [1e6] * 3
        
        self.gyroModel.PMatrixGyro = np.array(PMatrixGyro).reshape(3,3)
        self.gyroModel.walkBoundsGyro = np.array(errorBoundsGyro)
        self.gyroModel.PMatrixAccel = np.array(PMatrixAccel).reshape(3,3)
        self.gyroModel.walkBoundsAccel = np.array(errorBoundsAccel)

    def AddRWs(self):
        self.RW1 = self.rwFactory.create('Honeywell_HR16'
                               , [1, 0, 0]
                               , maxMomentum=50.
                               , Omega=0.  # RPM
                               , RWModel=self.varRWModel
                               )
        self.RW2 = self.rwFactory.create('Honeywell_HR16'
                               , [0, 1, 0]
                               , maxMomentum=50.
                               , Omega=0.  # RPM
                               , RWModel=self.varRWModel
                               )

        self.RW3 = self.rwFactory.create('Honeywell_HR16'
                               , [0, 0, 1]
                               , maxMomentum=50.
                               , Omega=0.  # RPM
                               , rWB_B=[0.5, 0.5, 0.5]  # meters
                               , RWModel=self.varRWModel
                               )

        self.numRW = self.rwFactory.getNumOfDevices()

        # create RW object container and tie to spacecraft object
        self.rwStateEffector = reactionWheelStateEffector.ReactionWheelStateEffector()
        self.rwFactory.addToSpacecraft("ReactionWheels", self.rwStateEffector, self.scObject)

    def AddExtForceTorque(self):
        self.extForceTorque.cmdTorqueInMsgName = "Lr_requested"
    def AddOpnavCamera(self):
        self.opnavCamera = opnavCamera.opnavCamera("opnavCamera")
        self.opnavCamera.ModelTag = "opnavCamera"
        self.opnavCamera.InputStateMsg = self.scObject.scStateOutMsgName
        self.opnavCamera.OutputDataMsg = "opnavCameraOutputMsg"


    # Global call to initialize every module
    def InitAllDynObjects(self):
        self.SetSpacecraftObject()
        self.SetSimpleNavObject()
        self.SetGravityBodies()
        self.SetSpiceObject()
        self.SetEphemerisConverter()
        self.SetSRPModel()
        #self.SetBeacons()
        self.AddStarTracker()
        #self.AddOpnavCamera()
        self.AddGyro()
        self.AddRWs()

