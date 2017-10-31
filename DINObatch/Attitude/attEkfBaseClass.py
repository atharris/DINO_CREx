import attEkfAlgs as ekf
import numpy as np
import scipy as sci
import attKalmanAlgs as aka

class attitudeEKF:
    def __init__(self, initState, initCovar, procNoise, measNoise, statePartials, directTransmission, measPartials, measTransmission, dt):
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

        self.F = statePartials
        self.G = directTransmission
        self.H = measPartials
        self.J = measTransmission

        self.dt = dt

        self.newMeasurements = False

    def newGyroMeas(self, gyroMeas):
        self.rateObs = gyroMeas
        self.newMeasurements = True

    def newStMeas(self, stMeas):
        self.self.stObs = stMeas
        self.newMeasurements = True

    def UpdateState(self):
        if self.newMeasurements:
            self.filterCorrect()
            self.newMeasurements = False
            self.currentMeas = np.zeros([6,])
        else:
            self.filterProp()

    def filterProp(self):
        stm = aka.linearizeSystem(self.estState[0:3],self.estState[4:6])
        self.estState = stm * self.estState
        self.estCovar = np.transpose(stm).dot(self.estCovar).dot(stm)

    def filterUpdate(self):
        # KF PREDICTION

        omega_ref = self.rateObs - self.estState[3:6]
        (A_LS, G_LS) = aka.linearizeSystem(self.estState[0:3], omega_ref)
        # Propagate KF State Estimate:
        params_KF = (self.rateObs, A_LS, G_LS, self.Q)
        predState = aka.RK4(self.estState, ekf.propagate_dynKF, params_KF, dt)
        aka.check_ShadowTransform(predState)

        # KF CORRECTION
        (self.estState, e_y) = aka.computeKalmanCorrection(predState, Y_meas_vec[:, k], H, R)
        E_y_vec[k] = e_y
        E_y_ergo_vec[k] = float(1.0/(k+1)) * np.sum(E_y_vec[:k])

        check_ShadowTransform(self.estState)


