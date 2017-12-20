import sys, os, inspect
import numpy as np
import math

bskName = 'Basilisk'
bskPath = '../..' + '/' + bskName + '/'
dinoPath = '../DINObatch/Attitude'
sys.path.append(bskPath + 'modules')
sys.path.append(bskPath + 'PythonModules')
sys.path.append(dinoPath)

try:
    import macros as mc
    # import batch_filter
    import ephem_difference
    import ephem_nav_converter
    # import simulation related support
    import spacecraftPlus
    import simIncludeGravBody
    import simIncludeRW
    import simple_nav
    import reactionWheelStateEffector
    import rwVoltageInterface

    # import FSW Algorithm related support
    import MRP_PD
    import inertial3D
    import attTrackingError
    import rwMotorTorque
    import fswSetupRW
    import rwMotorVoltage
    import simulationArchTypes
    import celestialTwoBodyPoint
    import hillPoint

    # import message declarations
    import fswMessages
except ImportError:
    import Basilisk.utilities.macros as mc
    from Basilisk.utilities import simIncludeRW, simulationArchTypes
    from Basilisk.simulation import sim_model, spacecraftPlus, gravityEffector, simple_nav, spice_interface
    from Basilisk.simulation import ephemeris_converter, radiation_pressure, star_tracker, imu_sensor
    from Basilisk.simulation import reactionWheelStateEffector, rwVoltageInterface
    from Basilisk.fswAlgorithms import ephem_difference, ephem_nav_converter
    from Basilisk.fswAlgorithms import MRP_PD, inertial3D, attTrackingError, rwMotorTorque, celestialTwoBodyPoint
    from Basilisk.fswAlgorithms import hillPoint
import AttitudeFilter as aekf

class FSWClass():
    def __init__(self, SimBase, updateRate=0.1):
        # Define process name, task name and task time-step
        self.processName = SimBase.FSWProcessName
        self.defaultTaskTimeStep = mc.sec2nano(updateRate)

        self.pyTaskName = "pyFswTask"
        self.taskName = "FswTask"

        # Create Python Tasks
        SimBase.fswPyProc.createPythonTask(self.pyTaskName, self.defaultTaskTimeStep,True, 30)

        #   Create normal tasks
        SimBase.dynProc.addTask(SimBase.CreateNewTask(self.taskName, self.defaultTaskTimeStep))

        self.attFilter = aekf.AttitudeFilter("attitudeFilter", True, 100)

        self.mrpControlConfig = MRP_PD.MRP_PDConfig()
        self.mrpControlWrap = SimBase.setModelDataWrap(self.mrpControlConfig)
        self.mrpControlWrap.ModelTag = "MRP_PD"

        self.rwMotorTorqueConfig = rwMotorTorque.rwMotorTorqueConfig()
        self.rwMotorTorqueWrap = SimBase.setModelDataWrap(self.rwMotorTorqueConfig)
        self.rwMotorTorqueWrap.ModelTag = "rwMotorTorque"

        self.attErrorConfig = attTrackingError.attTrackingErrorConfig()
        self.attErrorWrap = SimBase.setModelDataWrap(self.attErrorConfig)
        self.attErrorWrap.ModelTag = "attErrorInertial3D"

        self.attGuideConfig = hillPoint.hillPointConfig()
        self.attGuideWrap = SimBase.setModelDataWrap(self.attGuideConfig)
        self.attGuideWrap.ModelTag = "hillPoint"

        # Initialize all modules
        self.baseEphemeris = SimBase.DynClass.marsConvertName # ephemeris base for planet data and vehicle reference
        self.InitAllFSWObjects(SimBase)

        # Assign initialized modules to tasks
        #SimBase.AddModelToTask("ephemDiffConverterTask", self.ephemDifferenceConvWrap, self.ephemDifferenceConv, 10)
        #SimBase.AddModelToTask("vehicleConverterTask", self.vehicleEphConvWrap, self.vehicleEphConv, 9)
        SimBase.fswPyProc.addModelToTask(self.pyTaskName, self.attFilter)

        SimBase.AddModelToTask(self.taskName, self.mrpControlWrap, self.mrpControlConfig)
        SimBase.AddModelToTask(self.taskName, self.rwMotorTorqueWrap, self.rwMotorTorqueConfig)
        SimBase.AddModelToTask(self.taskName, self.attGuideWrap, self.attGuideConfig)
        SimBase.AddModelToTask(self.taskName, self.attErrorWrap, self.attErrorConfig)

    def SetPdController(self, SimBase):
        self.mrpControlConfig.inputGuidName = self.attErrorConfig.outputDataName
        self.mrpControlConfig.inputGuidName = "attErrorInertial3DMsg"
        self.mrpControlConfig.inputVehicleConfigDataName = "vehicleConfigName"
        self.mrpControlConfig.outputDataName = SimBase.DynClass.extForceTorque.cmdTorqueInMsgName
        self.mrpControlConfig.K = 3.5
        self.mrpControlConfig.P = 30.0
        return

    def SetAttGuidance(self, SimBase):
        self.attGuideConfig.outputDataName = "attGuideInertial"
        self.attGuideConfig.inputNavDataName = SimBase.DynClass.simpleNavObject.outputTransName
        self.attGuideConfig.inputCelMessName = SimBase.DynClass.earthGravBody.bodyInMsgName
        return


    def SetAttError(self, SimBase):
        self.attErrorConfig.outputDataName = "attErrorInertial3DMsg"
        self.attErrorConfig.inputRefName = self.attGuideConfig.outputDataName
        self.attErrorConfig.inputNavName = SimBase.DynClass.simpleNavObject.outputAttName #self.attFilter.outputMsgName
        return

    def SetAttitudeFilter(self):
        self.attFilter.ModelTag = "attitudeFilter"
        self.attFilter.inputIMUMsgName =  "gyro_output_data"
        self.attFilter.inputStMsgName = "star_tracker_state"
        self.attFilter.outputMsgName = "aekf_output_data"
        self.attFilter.filterMsgName = "aekf_filter_data"
        self.attFilter.dt = mc.NANO2SEC*self.defaultTaskTimeStep
        return

    def SetMotorTorqueConv(self, SimBase):
        self.rwMotorTorqueConfig.outputDataName = SimBase.DynClass.rwStateEffector.InputCmds
        self.rwMotorTorqueConfig.inputVehControlName = self.mrpControlConfig.outputDataName
        self.rwMotorTorqueConfig.rwParamsInMsgName = "rwa_config_data_parsed"
        # Make the RW control all three body axes
        controlAxes_B = [
            1, 0, 0
            , 0, 1, 0
            , 0, 0, 1
        ]
        self.rwMotorTorqueConfig.controlAxes_B = controlAxes_B

    def InitAllFSWObjects(self, SimBase):
        self.SetAttGuidance(SimBase)
        self.SetAttError(SimBase)
        self.SetPdController(SimBase)
        self.SetMotorTorqueConv(SimBase)
        self.SetAttitudeFilter()
        return