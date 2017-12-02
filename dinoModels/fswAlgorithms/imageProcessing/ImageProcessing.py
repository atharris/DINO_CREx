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
    def __init__(self, modelName, modelActive=True, modelPriority=-1):
        super(ImageProcessing, self).__init__(modelName, modelActive, modelPriority)

        ## Input guidance structure message name
        self.inputCamMsgName = "Cam Msg"
        self.inputBeaconMsgName = []        # create separate 'AddBeacon' function to add variable number of beacons

        ## Output body torque message name
        self.outputMsgName = ""

        ## Input message ID (initialized to -1 to break messaging if unset)
        self.inputCamMsgID = -1
        self.inputNavID = -1
        self.inputAttdeNavID = -1

        ## Output message ID (initialized to -1 to break messaging if unset)
        self.outputMsgID = -1

        ## Output Lr torque structure instantiation.  Note that this creates the whole class for interation
        # with the messaging system
        self.outputMsgData = other_module.OutputMsgStruct()

        ## Input vehicle configuration structure instantiation.  Note that this creates the whole class for interation
        # with the messaging system
        self.inputMsgData = other_module.OutputMsgStruct()

    ## The selfInit method is used to initialze all of the output messages of a class.
    # It is important that ALL outputs are initialized here so that other models can
    # subscribe to these messages in their crossInit method.
    def selfInit(self):
        self.outputMsgID = simulationArchTypes.CreateNewMessage(self.outputMsgName, self.outputMsgData,
                                                                 self.moduleID)
        return

    ## The crossInit method is used to initialize all of the input messages of a class.
    #  This subscription assumes that all of the other models present in a given simulation
    #  instance have initialized their messages during the selfInit step.
    def crossInit(self):
        self.inputMsgID = simulationArchTypes.SubscribeToMessage(self.inputMsgName, self.inputMsgData, self.moduleID)
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
        simulationArchTypes.ReadMessage(self.inputMsgID, self.inputMsgData, self.moduleID)
        #   Next, implement your routines or functions to process the input data and store it:
        localInputData = self.inputMsgData
        localOutputData = foo(localInputData)
        self.outputMsgData = localOutputData
        #   Finally, write the output message types:
        simulationArchTypes.WriteMessage(self.outputMsgID, currentTime, self.outputMsgData, self.moduleID)
