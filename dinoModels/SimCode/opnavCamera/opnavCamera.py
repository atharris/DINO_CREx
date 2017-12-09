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
bskPath = '../..' + '/' + bskName + '/'
sys.path.append(bskPath + 'modules')
sys.path.append(bskPath + 'PythonModules')
sys.path.append('../dinoModels/SimCode/opnavCamera/dependencies')

try:
    import simulationArchTypes
    import camera
    import spacecraftPlus
except ImportError:
    import Basilisk.utilities.simulationArchTypes as simulationArchTypes
    import Basilisk.simulation.spacecraftPlus as spacecraftPlus

import bodies as bod
from constants import au
import numpy as np

# if this script is run from a custom folder outside of the Basilisk folder, then uncomment the
# following line and specify the absolute bath to the Basilisk folder
#bskPath = '/Users/hp/Documents/Research/' + bskName + '/'
# @endcond

class opnavCamera(simulationArchTypes.PythonModelClass):
    def __init__(self, modelName, modelActive=True, modelPriority=-1):
        super(opnavCamera, self).__init__(modelName, modelActive, modelPriority)
        ## Input guidance structure message name
        self.scStateInputMsgName = "scStateInputMsgName"
        ## Output body torque message name
        self.outputMsgName = "opnavCameraOutputMsg"
        ## Input message ID (initialized to -1 to break messaging if unset)
        self.scStateInputMsgID = -1
        ## Output message ID (initialized to -1 to break messaging if unset)
        self.outputMsgID = -1
        import pdb
        pdb.set_trace()
        # Output Lr torque structure instantiation.  Note that this creates the whole class for interation
        # with the messaging system
        #self.outputMsgData = spice_interface.OutputMsgStruct()

        # Input vehicle configuration structure instantiation.  Note that this creates the whole class for interation
        # with the messaging system
        self.SCstate = spacecraftPlus.SCPlusStatesSimMsg()

        #self.inputMsgData = other_module.OutputMsgStruct()

        #load tranmission curve for Canon 20D
        _20D = np.load('../dinoModels/SimCode/opnavCamera/tc/20D.npz')
        tc = {}
        tc['lambda'] = _20D['x']
        tc['throughput'] = _20D['y']

        #load QE curve for Hubble Space Telecope Advanced Camera for Surveys SITe CCD
        ACS = np.load('../dinoModels/SimCode/opnavCamera/qe/ACS.npz')
        qe = {}
        qe['lambda'] = ACS['x']
        qe['throughput'] = ACS['y']

        scState = np.array([au/1000-100000,0,0,0,0,0])
        scDCM = np.identity(3)
        bod.earth.state = np.array([au/1000, 0,0,0,0,0])
        bod.luna.state = np.array([au/1000, 250000,0,0,0,0])
        takeImage = 0
        bodies = [bod.earth, bod.luna]

        self.cam = camera.camera(
            2,                  #detectorHeight
            2,                  #detectorWidth
            5.0,                #focalLength
            512,                #resolutionHeight
            512,                #resolutionWidth
            np.identity(3),     #body2cameraDCM
            1000,               #maximum magnitude
            -1000,              #minimum magnitude (for debugging)
            qe,                 #quantum efficiency dictionary
            tc,                 #transmission curve dictionary
            1,                  #wavelength bin size in nm
            0.01**2,            #effective area in m^2
            100,                #dark current in electrons per second
            100,                #std for read noise in electrons
            100,                #bin size
            2**32,              #max bin depth
            1,                  #sigma for gaussian psf
            0.001,               #simulation timestep
            scState,            #position state of s/c
            scDCM,              #intertal 2 body DCM for s/c
            bodies,             #bodies to track in images
            takeImage,          #takeImage message
            db='../dinoModels/SimCode/opnavCamera/db/tycho.db'
            )

    ## The selfInit method is used to initialze all of the output messages of a class.
    # It is important that ALL outputs are initialized here so that other models can
    # subscribe to these messages in their crossInit method.
    def selfInit(self):
        #self.outputMsgID = simulationArchTypes.CreateNewMessage(self.outputMsgName, self.outputMsgData,
        #                                                        self.moduleID)
        return

    ## The crossInit method is used to initialize all of the input messages of a class.
    #  This subscription assumes that all of the other models present in a given simulation
    #  instance have initialized their messages during the selfInit step.
    def crossInit(self):
        import pdb
        print('############################################################')
        pdb.set_trace()
        self.scStateInputMsgID = simulationArchTypes.SubscribeToMessage(
            self.scStateInputMsgName, self.SCstate, self.moduleID)

        # self.inputMsgID = simulationArchTypes.SubscribeToMessage(
        #     self.inputMsgName, self.SCstate, self.moduleID)
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

        print('img')

        simulationArchTypes.ReadMessage(
            self.scStateInputMsgID, self.SCstate, self.moduleID)
        print self.SCstate.r_BN_N
        import pdb
        pdb.set_trace()
        # self.msg['takeImage']=1
        # self.cam.updateState()
        # self.msg['takeImage']=0
        # self.cam.updateState()

        #   First, read messages we've subscribed to:
        # simulationArchTypes.ReadMessage(self.inputMsgID, self.inputMsgData, self.moduleID)
        # #   Next, implement your routines or functions to process the input data and store it:
        # localInputData = self.inputMsgData
        # localOutputData = foo(localInputData)
        # self.outputMsgData = localOutputData
        # #   Finally, write the output message types:
        # simulationArchTypes.WriteMessage(self.outputMsgID, currentTime, self.outputMsgData, self.moduleID)
