import pytest
import sys, os, inspect
import matplotlib
import numpy as np
import ctypes
import math
import csv
import logging

# @cond DOXYGEN_IGNORE
filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))

bskName = 'Basilisk'
bskPath = '../../..' + '/' + bskName + '/'
sys.path.append(bskPath + 'modules')
sys.path.append(bskPath + 'PythonModules')
import import simple_nav

navPath = '../../../DINObatch/P\&L/common_functions'
import batchFilter
import posVel

# if this script is run from a custom folder outside of the Basilisk folder, then uncomment the
# following line and specify the absolute bath to the Basilisk folder
#bskPath = '/Users/hp/Documents/Research/' + bskName + '/'
# @endcond

class FswStateProp(simulationArchTypes.PythonModelClass):
    def __init__(self, modelName, modelActive=True, modelPriority=-1):
        super(PythonModule, self).__init__(modelName, modelActive, modelPriority)

        ## Output body torque message name
        self.outputMsgName = "fswPropStateMsg"
        ## Output message ID (initialized to -1 to break messaging if unset)
        self.outputMsgID = -1

        self.internalState = np.zeros([6,])
        self.t0 = 0.0   #   Time offset in nanoseconds.
        self.currPropTime = 0.0 #   Current time.

        ## Output Lr torque structure instantiation.  Note that this creates the whole class for interation
        # with the messaging system
        self.outputMsgData = simple_nav.NavTransIntMsg()


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

        self.currPropTime = self.t0 + currentTime
        propagatorInput = (self.internalState, self.internalSTM, self.currPropTime-self.prevPropTime, self.extras)
        self.internalState = batchFilter.runRef(propagatorInput)

        self.outputMsgData.r_BN_N = self.internalState
        self.prevPropTime = self.currPropTime
        #   Finally, write the output message types:
        simulationArchTypes.WriteMessage(self.outputMsgID, currentTime, self.outputMsgData, self.moduleID)
