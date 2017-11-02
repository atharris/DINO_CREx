import attEkfAlgs as ekf
import numpy as np
from numpy import linalg as la
import scipy as sci
import attKalmanAlgs as aka
import RigidBodyKinematics as rbk

class attitudeEKF:
    def __init__(self, initState, initCovar, procNoise, measNoise, measPartials, measTransmission, dt):
        #! Function spec:
        #! Inputs
        #! initState - 6x1 - initial state vector; 1st 3 are MRP components, second 3 are body/inertial angular rates
        #! initCovar - 6x6 - initial covariance matrix, used as a tuning parameter.
        #! procNoise - 6x6 - process noise matrix; used as a tuning paramter.
        #! measNoise - 6x6 - measurement noise matrix; used as a tuning paramter.
        #! statePartials - 6x6 - state transition matrix F used to propagate dynamics within the filter.
        #! directTranmission - 6x3 - G matrix used to transform from the actuator space into the state space. Can be set to zero if no actuator effects are considered.
        #! measPartials - 6x6 - partials of measurements with respect to the states. User-spec-able.
        #! measTransmission - 6x3 - partials of measurements with respect to the control inputs. Typically will be zero.

        self.estState = initState
        self.estCovar = initCovar

        self.Q = procNoise
        self.R = measNoise

        self.H = measPartials
        self.J = measTransmission

        self.dt = dt

        self.newMeasurements = False
        self.newGyro = False
        self.newSt = False

    def newGyroMeas(self, gyroMeas):
        self.rateObs = gyroMeas
        self.newGyro = True

    def newStMeas(self, stMeas):
        self.stObs = stMeas
        self.newSt = True

    def UpdateState(self):
        self.newMeasurements = self.newSt or self.newGyro
        if self.newMeasurements:
            self.filterUpdate(self.newSt, self.newGyro)
            self.newMeasurements = False
            self.newSt = False
            self.newGyro = False
        else:
            self.filterProp()

    def filterProp(self):
        stm = aka.linearizeSystem(self.estState[0:3],self.estState[4:6])
        self.estState = stm * self.estState
        self.estCovar = np.transpose(stm).dot(self.estCovar).dot(stm)

    def filterUpdate(self, newSt, newGyro):
        # KF PREDICTION
        omega_ref = self.rateObs - self.estState[3:6]
        (A_LS, G_LS) = aka.linearizeSystem(self.estState[0:3], omega_ref)

        # Propagate KF State Estimate:
        params_KF = (self.rateObs, A_LS, G_LS, self.Q)
        predState = aka.RK4(self.estState, ekf.propagate_dynKF, params_KF, self.dt)
        aka.check_ShadowTransform(predState)

        # KF CORRECTION
        if newSt and newGyro:
            tmpMeas = self.stMeas
            tmpH = self.H
            tmpR = self.R
            (self.estState, self.estCovar) = aka.computeKalmanCorrection(predState, tmpMeas, tmpH, tmpR)
        elif newGyro:
            tmpMeas = self.stMeas
            self.estState[0:3] = predState[0:3]
        else:
            print "Error: This shouldn't have happened."

        (self.estState, self.estCovar) = aka.computeKalmanCorrection(predState, tmpMeas, tmpH, tmpR)
        rbk.MRPShadow(self.estState, 1.0)

if __name__ == "__main__":
    (k_f, dt, k_meas) = aka.define_numericData()
    print 't_sim = ', k_f*dt

    (cov_ARW, cov_RRW) = aka.define_noiseIntensities()
    print 'cov_ARW = ', cov_ARW
    print 'cov_RRW = ', cov_RRW

    (sigma0_true, omega_true, bias_true) = aka.initialize_groundTruth()
    X0_ground = np.append(sigma0_true, bias_true)
    print 'sigma0_true = ', sigma0_true
    print 'omega_true = ', omega_true
    print 'bias_true = ', bias_true
    W = ekf.define_PSD_Q()
    X0_est = ekf.define_initialGuess_KF()
                    #self, initState, initCovar, procNoise, measNoise, statePartials, directTransmission, measPartials, measTransmission, dt
    ekfClass = attitudeEKF(X0_est, np.identity(6), cov_ARW, cov_RRW,