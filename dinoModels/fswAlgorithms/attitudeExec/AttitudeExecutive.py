import pytest
import sys, os, inspect
import matplotlib
import numpy as np
import ctypes
import math
import csv
import logging

import AttExec as ae

# @cond DOXYGEN_IGNORE
filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))

bskName = 'Basilisk'
bskPath = '../../../..' + '/' + bskName + '/'
sys.path.append(bskPath + 'modules')
sys.path.append(bskPath + 'PythonModules')

import spice_interface
import simple_nav


# if this script is run from a custom folder outside of the Basilisk folder, then uncomment the
# following line and specify the absolute bath to the Basilisk folder
#bskPath = '/Users/hp/Documents/Research/' + bskName + '/'
# @endcond
class PythonMRPPD(simulationArchTypes.PythonModelClass):
    def __init__(self, modelName, modelActive=True, modelPriority=-1):
        super(PythonMRPPD, self).__init__(modelName, modelActive, modelPriority)

        ## Proportional gain term used in control
        self.K = 0
        ## Derivative gain term used in control
        self.P = 0
        ##  Convergence Counter
        self.convCounter = 0
        self.convTarget = 10
        self.obsCounter = 0
        self.obsTime = 30
        self.targInd = 0
        self.targList = []
        self.sysMode = 0 #  0 for observation, 1 if in safe mode. Needs to be reset.

        ## Input guidance structure message name
        self.inputGuidName = ""
        ## Input vehicle configuration structure message name
        self.inputVehicleConfigDataName = ""
        ##  Input target SPICE data msg names
        self.inputTargetMessageNames = []
        ## Output body torque message name
        self.outputDataName = ""

        ## Input message ID (initialized to -1 to break messaging if unset)
        self.inputGuidID = -1
        ## Input message ID (initialized to -1 to break messaging if unset)
        self.inputVehConfigID = -1
        self.inputTargetIDs = [-1]
        ## Output message ID (initialized to -1 to break messaging if unset)
        self.outputDataID = -1


        ## Output Lr torque structure instantiation.  Note that this creates the whole class for interation
        # with the messaging system
        self.outputMessageData = MRP_PD.CmdTorqueBodyIntMsg()
        ## Input guidance error structure instantiation.  Note that this creates the whole class for interation
        # with the messaging system
        self.inputGuidMsg = MRP_PD.AttGuidFswMsg()
        ## Input vehicle configuration structure instantiation.  Note that this creates the whole class for interation
        # with the messaging system
        self.inputConfigMsg = MRP_PD.VehicleConfigFswMsg()

    ## The selfInit method is used to initialze all of the output messages of a class.
    # It is important that ALL outputs are initialized here so that other models can
    # subscribe to these messages in their crossInit method.
    def selfInit(self):
        self.outputDataID = simulationArchTypes.CreateNewMessage(self.outputDataName, self.outputMessageData,
                                                                 self.moduleID)
        return

    ## The crossInit method is used to initialize all of the input messages of a class.
    #  This subscription assumes that all of the other models present in a given simulation
    #  instance have initialized their messages during the selfInit step.
    def crossInit(self):

        self.inputGuidID = simulationArchTypes.SubscribeToMessage(self.inputGuidName, self.inputGuidMsg, self.moduleID)
        self.inputVehConfigID = simulationArchTypes.SubscribeToMessage(self.inputVehicleConfigDataName,
                                                                       self.inputConfigMsg, self.moduleID)
        return

    ## The reset method is used to clear out any persistent variables that need to get changed
    #  when a task is restarted.  This method is typically only called once after selfInit/crossInit,
    #  but it should be written to allow the user to call it multiple times if necessary.
    def reset(self, currentTime):
        self.sysMode = 0
        return

    def addTarget(self, targetMsgName):
        self.targetList.append(self.targetList[-1]+1)
        self.inputTargetMessageNames.append(targetMsgName)

    ## The updateState method is the cyclical worker method for a given Basilisk class.  It
    # will get called periodically at the rate specified in the Python task that the model is
    # attached to.  It persists and anything can be done inside of it.  If you have realtime
    # requirements though, be careful about how much processing you put into a Python updateState
    # method.  You could easily detonate your sim's ability to run in realtime.
    def updateState(self, currentTime):
        simulationArchTypes.ReadMessage(self.inputGuidID, self.inputGuidMsg, self.moduleID)
        lrCmd = np.array(self.inputGuidMsg.sigma_BR) * self.K + np.array(self.inputGuidMsg.omega_BR_B) * self.P
        self.outputMessageData.torqueRequestBody = (-lrCmd).tolist()
        simulationArchTypes.WriteMessage(self.outputDataID, currentTime, self.outputMessageData, self.moduleID)

        return

    def findTargetMrp(self):

        return

    def checkIfObserving(self):
        if self.sigmaError < self.sigmaThreshold and self.rateError < self.rateThreshold:
            self.convCounter = self.convCounter+1
        else:
            convCounter = 0

        if self.convCounter > self.convTarget and self.obsCounter < self.obsTime:
            self.takeImage = 1
        else:
            self.takeImage = 0

        if self.obsCounter > self.obsTime:
            self.obsCounter = 0
            self.convCounter = 0
            try:
                self.targInd = self.targInd+1
            except:
                self.sysMode = 1
                print "Entering safe mode!"
