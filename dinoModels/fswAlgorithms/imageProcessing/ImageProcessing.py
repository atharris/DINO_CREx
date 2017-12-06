import pytest
import sys, os, inspect
import matplotlib
import numpy as np
import ctypes
import math
import csv
import logging

#imu_sensor.h type file in         dinoModels/messages
import NavMsg       #msg struct     navTransIntMsg
import CamMsg       #msg struct
import AttdeNavMsg  #msg struct     navAttIntMsg
import spiceBeacon  #msg struct     beacon positions
import ImageProcessingExecutive as IP

# @cond DOXYGEN_IGNORE
filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))

bskName = 'Basilisk'
bskPath = '../..' + '/' + bskName + '/'
sys.path.append(bskPath + 'modules')
sys.path.append(bskPath + 'PythonModules')

# if this script is run from a custom folder outside of the Basilisk folder, then uncomment the
# following line and specify the absolute bath to the Basilisk folder
# bskPath = '/Users/hp/Documents/Research/' + bskName + '/'
# @endcond

class ImageProcessing(simulationArchTypes.PythonModelClass):
    def __init__(self, modelName, modelActive=True, modelPriority=-1, cameraParameters, beaconRadius):
        super(ImageProcessing, self).__init__(modelName, modelActive, modelPriority)

        ## Input message names
        self.inputCamMsgName = "Cam Msg"
        self.inputBeaconMsgName = []        # individual beacon messages added in addBeacon() function
        self.inputNavMsgName = "Nav Msg"
        self.inputAttdeNavMsgName = "Nav Attde Msg"

        ## Output image processing message name
        self.outputMsgName = "IP Msg"

        ## Input message ID (initialized to -1 to break messaging if unset)
        self.inputCamMsgID = -1
        self.inputBeaconMsgID = []                   # check if this is right ...
        self.inputNavMsgID = -1
        self.inputAttdeNavMsgID = -1            # check if this is right ...

        ## Output message ID (initialized to -1 to break messaging if unset)
        self.outputMsgID = -1

        ## Output message structure instantiation.  Note that this creates the whole class for interation
        # with the messaging system
        # self.outputMsgData = other_module.OutputMsgStruct()
        self.outputMsgData = ImageProcessing.OutputMsgStruct()

        ## Input message structure instantiation.  Note that this creates the whole class for iteration
        # with the messaging system
        # self.inputMsgData = other_module.OutputMsgStruct()
        self.inputCamMsgData = CameraModule.OutputMsgStruct()
        self.inputNavMsgData = NavModule.OutputMsgStruct()
        self.inputAttdeNavMsgData = AttitudeNavModule.OutputMsgStruct()
        self.inputBeaconMsgData = []

        # python dict with camera parameters
        self.cameraParameters = cameraParameters

        # python dict with beacon radius
        self.beaconRadius = beaconRadius

    ## The selfInit method is used to initialize all of the output messages of a class.
    # It is important that ALL outputs are initialized here so that other models can
    # subscribe to these messages in their crossInit method.
    def selfInit(self):

        # ??? What are simulationArchTypes?
        self.outputMsgID = simulationArchTypes.CreateNewMessage(self.outputMsgName, self.outputMsgData,
                                                                 self.moduleID)
        return


    ## The addBeacon method is used to intialize input messages for beacon positions.
    #  This subscription assumes that all of the other models present in a given simulation
    #  instance have initialized their messages during the selfInit step.
    def addBeacon(self, beaconID):

        self.inputBeaconMsgID.append(-1)
        self.inputBeaconMsgName.append('Beacon ID '+beaconID+' Msg')
        # need to find appropriate Beacon Module
        self.inputBeaconMsgData.append(BeaconModule.OutputMsgStruct())

        return


    ## The crossInit method is used to initialize all of the input messages of a class.
    #  This subscription assumes that all of the other models present in a given simulation
    #  instance have initialized their messages during the selfInit step.
    def crossInit(self):

        # original command
        # self.inputMsgID = simulationArchTypes.SubscribeToMessage(
        #   self.inputMsgName, self.inputMsgData, self.moduleID)

        # self.moduleID need updating?
        self.inputCamMsgID = simulationArchTypes.SubscribeToMessage(
            self.inputCamMsgName, self.inputCamMsgData, self.moduleID)

        self.inputAttdeNavID = simulationArchTypes.SubscribeToMessage(
            self.inputAttdeNavMsgName, self.inputAttdeNavMsgData, self.moduleID)

        self.inputNavMsgID = simulationArchTypes.SubscribeToMessage(
            self.inputNavMsgName, self.inputNavMsgData, self.moduleID)

        numBeacons = len(self.inputBeaconMsgID)
        for indBeacon in range(numBeacons):
            self.inputBeaconMsgID[indBeacon] = simulationArchType.SubscribeToMessage(
                self.inputBeaconMsgName, self.inputBeaconMsgData, self.moduleID)

        return

    ## The reset method is used to clear out any persistent variables that need to get changed
    #  when a task is restarted.  This method is typically only called once after selfInit/crossInit,
    #  but it should be written to allow the user to call it multiple times if necessary.
    def reset(self, currentTime):
        return

    ## The updateState method is the cyclical worker method for a given Basilisk class.  It
    # will get called periodically at the rate specified in the Python task that the model is
    # attached to.  It persists and anything can be done inside of it.  If you have realtime
    # requirements though, be careful about how much processing you put into a Python updateState
    # method.  You could easily detonate your sim's ability to run in realtime.
    def updateState(self, currentTime):

        #   First, read messages we've subscribed to:
        # simulationArchTypes.ReadMessage(self.inputMsgID, self.inputMsgData, self.moduleID)
        simulationArchTypes.ReadMessage(self.inputCamMsgID, self.inputCamMsgData, self.moduleID)
        simulationArchTypes.ReadMessage(self.inputNavMsgID, self.inputNavMsgData, self.moduleID)
        simulationArchTypes.ReadMessage(self.inputAttdeNavID, self.inputAttdeNavData, self.moduleID)

        numBeacons = len(self.inputBeaconMsgID)
        for indBeacon in range(numBeacons):
            simulationArchTypes.ReadMessage(self.inputBeaconMsgID[indBeacon],
                                            self.inputBeaconMsgData[indBeacon],
                                            self.moduleID)

        #   Next, implement your routines or functions to process the input data and store it:
        currentImageMap = self.inputCamMsgData.imageMap

        r_N_currentBeacons = []
        for indBeacon in range(numBeacons):
            r_N_currentBeacons.append(self.inputBeaconMsgData[indBeacon].r_N)   # verify beacon message contents

        r_N_sc_estimate = self.inputNavMsgData.r_N                              # verify nav message contents
        sigma_BN_estimate = self.inputAttdeNavMsgData.sigma_BN                  # verify attde nav message contents

        # image processing executive function
        idOutput, pixelLineBeaconOutput, sigmaOutput = IP.imageProcessing(
            currentImageMap, self.cameraParameters, r_N_sc_estimate,
            sigma_BN_estimate, r_N_currentBeacons, self.beaconRadius)

        # self.outputMsgData = localOutputData
        self.outputMsgData.beaconID = idOutput
        self.outputMsgData.timeImage = self.inputCamMsgData.timeImage
        self.outputMsgData.beaconPL = pixelLineBeaconOutput
        self.outputMsgData.sigma_BN = sigma_BN_estimate

        #   Finally, write the output message types:
        simulationArchTypes.WriteMessage(self.outputMsgID, currentTime, self.outputMsgData, self.moduleID)

        return
