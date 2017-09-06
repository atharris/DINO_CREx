import sys, os, inspect

bskName = 'Basilisk'
bskPath = '../..' + '/' + bskName + '/'
sys.path.append(bskPath + 'modules')
sys.path.append(bskPath + 'PythonModules')

import macros as mc
import batch_filter
import ephem_difference
import ephem_nav_converter

class FSWClass():
    def __init__(self, SimBase):
        # Define process name, task name and task time-step
        self.processName = SimBase.FSWProcessName
        self.defaultTaskTimeStep = mc.sec2nano(0.1)

        # Create Tasks
        SimBase.fswProc.addTask(SimBase.CreateNewTask("vehicleConverterTask", self.defaultTaskTimeStep), 12)
        SimBase.fswProc.addTask(SimBase.CreateNewTask("ephemDiffConverterTask", self.defaultTaskTimeStep), 11)
        SimBase.fswProc.addTask(SimBase.CreateNewTask("batchFilterTask", self.defaultTaskTimeStep), 10)

        # Create module data and module wraps
        self.batchFilterData = batch_filter.BatchConfig()
        self.batchFilterWrap = SimBase.setModelDataWrap(self.batchFilterData)
        self.batchFilterWrap.ModelTag = "batchFilter"

        self.ephemDifferenceConv = ephem_difference.EphemDifferenceData()
        self.ephemDifferenceConvWrap = SimBase.setModelDataWrap(self.ephemDifferenceConv)
        self.ephemDifferenceConvWrap.ModelTag = "ephemerisDifferenceConverter"

        self.vehicleEphConv = ephem_nav_converter.EphemNavConverterData()
        self.vehicleEphConvWrap = SimBase.setModelDataWrap(self.vehicleEphConv)
        self.vehicleEphConvWrap.ModelTag = "vehicleEphemerisConverter"

        # Initialize all modules
        self.baseEphemeris = SimBase.DynClass.marsConvertName # ephemeris base for planet data and vehicle reference
        self.InitAllFSWObjects(SimBase)

        # Assign initialized modules to tasks
        SimBase.AddModelToTask("ephemDiffConverterTask", self.ephemDifferenceConvWrap, self.ephemDifferenceConv, 10)
        #SimBase.AddModelToTask("vehicleConverterTask", self.vehicleEphConvWrap, self.vehicleEphConv, 9)

    def SetEphemDifferenceConverter(self, SimBase):
        self.ephemDifferenceConv.ephBaseInMsgName = self.baseEphemeris
        self.ephemDifferenceConv.ephBdyCount = 4
        self.ephemDifferenceConv.baseScale = 1.0

        self.ephemDifferenceConv.changeBodies[0].ephInMsgName = SimBase.DynClass.marsConvertName
        self.ephemDifferenceConv.changeBodies[0].ephOutMsgName = 'mars_state_corr_eph'
        self.ephemDifferenceConv.changeBodies[1].ephInMsgName = SimBase.DynClass.earthConvertName
        self.ephemDifferenceConv.changeBodies[1].ephOutMsgName = 'earth_state_corr_eph'
        self.ephemDifferenceConv.changeBodies[2].ephInMsgName = SimBase.DynClass.sunConvertName
        self.ephemDifferenceConv.changeBodies[2].ephOutMsgName = 'sun_state_corr_eph'
        self.ephemDifferenceConv.changeBodies[3].ephInMsgName = SimBase.DynClass.moonConvertName
        self.ephemDifferenceConv.changeBodies[3].ephOutMsgName = 'moon_state_corr_eph'

        #self.ephemDifferenceConv.changeBodies[4].ephInMsgName = "vehicle_state_eph"
        #self.ephemDifferenceConv.changeBodies[4].ephOutMsgName = 'vehicle_state_corr_eph'

    def SetVehicleEphemNavConverter(self):
        self.vehicleEphConv.ephInMsgName = self.baseEphemeris
        self.vehicleEphConv.stateOutMsgName = "vehicle_eph_state_est"

    def SetBatchODFilter(self):
        self.batchFilterData.navStateOutMsgName = "trans_filter_nav_state"
        self.batchFilterData.filtDataOutMsgName = "trans_filter_data"

        self.batchFilterData.X0_star = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        self.batchFilterData.x0_bar = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

        self.batchFilterData.P0_bar = [
            1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 1.0
        ]

    def InitAllFSWObjects(self, SimBase):
        self.SetEphemDifferenceConverter(SimBase)
        #self.SetVehicleEphemNavConverter()