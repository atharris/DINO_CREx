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

import ekfFunctions as ekf
try:
    import simulationArchTypes
    import imu_sensor
    import spacecraftPlus
    import star_tracker
    import simple_nav
    import RigidBodyKinematics as rbk
except ImportError:
    from Basilisk.utilities import simulationArchTypes
    from Basilisk.fswAlgorithms import fswMessages
    import Basilisk.utilities.RigidBodyKinematics as rbk
    from Basilisk.simulation import spacecraftPlus, imu_sensor, star_tracker, simple_nav


# if this script is run from a custom folder outside of the Basilisk folder, then uncomment the
# following line and specify the absolute bath to the Basilisk folder
#bskPath = '/Users/hp/Documents/Research/' + bskName + '/'
# @endcond

##  AttitudeFilter Class
#   Inputs:
#   IMU_data - IMUmessage - used for gyro reading.
#   st_data - StarTracker message - used for star tracker reading.

#   Outputs:
#   nav_aggregate - provides estimated attitude, angular rates.

class AttitudeFilter(simulationArchTypes.PythonModelClass):
    def __init__(self, modelName, modelActive=True, modelPriority=-1):
        super(AttitudeFilter, self).__init__(modelName, modelActive, modelPriority)

        ## Input gyro, star tracker message names
        self.inputStMsgName = "star_tracker_state"
        self.inputIMUMsgName = "won_t_work"

        ## Output body torque message name
        self.outputMsgName = "aekf_output_data"
        self.filterMsgName = "aekf_filter_data"

        ## Input message ID (initialized to -1 to break messaging if unset)
        self.inputStMsgID = -1
        self.inputIMUMsgID = -1

        ## Output message ID (initialized to -1 to break messaging if unset)
        self.outputMsgID = -1
        self.filterMsgID = -1

        ## Input IMU, Star Tracker structures
        self.inputIMUMsgData = imu_sensor.IMUSensorIntMsg()
        self.inputStMsgData = star_tracker.STSensorIntMsg()

        ## Output navigation estimate structure.
        self.outputMsgData = simple_nav.NavAttIntMsg()
        self.filterMsgData = simple_nav.NavAttIntMsg()

        ##  Define Estimate variables
        self.stateEst = np.array([0.1,0.,0.,0.2,0.,0.]) #    state estimate is 3 delta MRPs, 3 bias states
        self.outEst = np.zeros(6) #      Output is 3 BN MRPs, 3 angular rates
        self.estCov = np.identity(6) #     Estimated state covariance

        ##  Define noise variables
        self.stateNoise = np.identity(6)
        self.stateNoise[0:3,0:3] = 0.00001**2*np.identity(3)
        self.stateNoise[3:6,3:6] = 0.01**2*np.identity(3)
        self.measNoise = 0.001*np.identity(3)

        self.postFitRes = 0.

        ##  Define system models
        self.H = np.zeros([3, 6])
        self.H[:3, :3] = np.identity(3)

        self.dt = 0.1

        A_ls, G_ls = ekf.linearizeSystem(self.stateEst[0:3], self.stateEst[3:6])

        self.predOptions = ekf.biasReplacementOptions(np.zeros([3,]), A_ls, G_ls, self.stateNoise)

    ## The selfInit method is used to initialze all of the output messages of a class.
    # It is important that ALL outputs are initialized here so that other models can
    # subscribe to these messages in their crossInit method.

    def selfInit(self):
        print "selfing:"
        print self.outputMsgName
        print self.moduleID
        self.outputMsgID = simulationArchTypes.CreateNewMessage(self.outputMsgName, self.outputMsgData,
                                                                 self.moduleID)
        self.filterMsgID = simulationArchTypes.CreateNewMessage(self.filterMsgName, self.filterMsgData,
                                                                 self.moduleID)
        print "Output AEKF ID:", self.outputMsgID
        print "Output AEKF Filter ID:", self.filterMsgID
        return

    ## The crossInit method is used to initialize all of the input messages of a class.
    #  This subscription assumes that all of the other models present in a given simulation
    #  instance have initialized their messages during the selfInit step.
    def crossInit(self):
        self.inputStMsgID = simulationArchTypes.SubscribeToMessage(self.inputStMsgName, self.inputStMsgData, self.moduleID)
        self.inputIMUMsgID = simulationArchTypes.SubscribeToMessage(self.inputIMUMsgName, self.inputIMUMsgData, self.moduleID)

        print "Input ST ID:", self.inputStMsgID
        print "Input IMU ID:", self.inputIMUMsgID
        return

    ## The reset method is used to clear out any persistent variables that need to get changed
    #  when a task is restarted.  This method is typically only called once after selfInit/crossInit,
    #  but it should be written to allow the user to call it multiple times if necessary.
    def reset(self, currentTime):
        return

    def readMsg(self):
        simulationArchTypes.ReadMessage(self.inputMsgID, self.inputMsgData, self.moduleID)
        simulationArchTypes.ReadMessage(self.inputMsgID, self.inputMsgData, self.moduleID)

    def kalmanStep(self):
        #   Prediction Step: predict the new mean and covariance based on the assumed model
        self.predOptions.omega_bn_meas = self.gyroMeas
        #print "State Estimate:", self.stateEst
        F, G = ekf.linearizeSystem(self.stateEst[0:3], self.gyroMeas - self.stateEst[3:6])
        self.predOptions.F = F
        self.predOptions.G = G
        self.predOptions.Q = self.stateNoise

        propState = np.zeros([42,])
        propState[0:6] = self.stateEst
        propState[6:] = np.reshape(self.estCov,[36,])

        predState = ekf.rk4(ekf.biasReplacementEOM, 0.0, self.dt, propState, self.predOptions)
        predCov = np.resize(predState[6:], (6,6))
        #   Correction Step: If new star tracker measurements are available, attempt to correct the attitude and bias estimates
        predState = predState[0:6]
        #   Compute the new Kalman gain
        temp = self.H.dot(predCov.dot(np.transpose(self.H))) + self.measNoise
        tempInv = np.linalg.inv(temp)
        kalmanGain = predCov.dot(np.transpose(self.H).dot(tempInv))

        #   Compute measurement innovation
        innov = self.stMeas - np.dot(self.H, predState)
        if np.linalg.norm(innov) > (1.0 / 3.0):
            innov_s = rbk.MRPswitch(self.stMeas, 1.0) - np.dot(self.H, predState)
            if np.linalg.norm(innov_s) < np.linalg.norm(innov):
                innov = innov_s


        newStateEst = predState + kalmanGain.dot(innov)

        KH = np.dot(kalmanGain, self.H)
        I = np.identity(KH.shape[0])
        M = I - KH
        estCov = np.dot(M, np.dot(predCov, M.T)) + np.dot(kalmanGain, np.dot(self.measNoise, kalmanGain.T))

        self.stateEst = newStateEst
        self.estCov = estCov

    ## The updateState method is the cyclical worker method for a given Basilisk class.  It
    # will get called periodically at the rate specified in the Python task that the model is
    # attached to.  It persists and anything can be done inside of it.  If you have realtime
    # requirements though, be careful about how much processing you put into a Python updateState
    # method.  You could easily detonate your sim's ability to run in realtime.
    def updateState(self, currentTime):
        #   First, read messages we've subscribed to:
        simulationArchTypes.ReadMessage(self.inputStMsgID, self.inputStMsgData, self.moduleID)
        simulationArchTypes.ReadMessage(self.inputIMUMsgID, self.inputIMUMsgData, self.moduleID)

        self.stMeas = rbk.EP2MRP(self.inputStMsgData.qInrtl2Case)
        self.gyroMeas =  self.inputIMUMsgData.AngVelPlatform
        #   Next, implement your routines or functions to process the input data and store it:
        self.kalmanStep()

        localEstAtt = self.stateEst[0:3]
        localRateEst = self.gyroMeas - self.stateEst[3:]

        self.outputMsgData.timeTag = currentTime
        self.outputMsgData.sigma_BN = localEstAtt
        self.outputMsgData.omega_BN_B = localRateEst

        #   Finally, write the output message types:
        simulationArchTypes.WriteMessage(self.outputMsgID, currentTime, self.outputMsgData, self.moduleID)

        # Write variable for post fit residuals
        # This is a hack because of the messaging format currently and the different versions of BSK
        # The filter message can't quite be used, so the filter data is being put in another msg
        self.postFitRes = self.stMeas - np.dot(self.H, self.stateEst)
        covarDiag1 = np.array([self.estCov[0,0],self.estCov[1,1],self.estCov[2,2]])
        covarDiag2 = np.array([self.estCov[3,3],self.estCov[4,4],self.estCov[5,5]])

        self.filterMsgData.timeTag = currentTime
        self.filterMsgData.sigma_BN = covarDiag1
        self.filterMsgData.omega_BN_B = self.postFitRes
        self.filterMsgData.vehSunPntBdy = covarDiag2
        simulationArchTypes.WriteMessage(self.filterMsgID, currentTime, self.filterMsgData, self.moduleID)








