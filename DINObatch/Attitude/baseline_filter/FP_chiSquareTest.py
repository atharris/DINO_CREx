import numpy as np
import matplotlib.pyplot as plt
import KF_commonAlgs as kf
import FP_extendedKalman as ekf

import sys, os
userPath = os.path.expanduser('~')
sys.path.append(userPath + '/Desktop/astroPy/')
import RigidBodyKinematics as rbk


def plot_MC_run():
    return
def create_N_monteCarloRuns(N_MC, P_groundTruth, Xm_groundTruth, X0_est, k_f, dt, k_meas, omega_true):
    print 'N_MC = ', N_MC
    np.random.seed(100)
    X_truthModel_samp = np.random.multivariate_normal(Xm_groundTruth, P_groundTruth, size=N_MC)
    for i in range(N_MC):
        X0_ground = X_truthModel_samp[i]
        X0_ground[0:3] = rbk.sigmaUnitSwitchCheck(X0_ground[0:3])
        print 'i_MC = ', i
        print 'X0_ground = ', X0_ground
        (k_vec, X_state_vec, X_std_vec, X_error_vec) = ekf.create_monteCarlo(X0_ground, X0_est, k_f, dt)
    kf.checkArePlotsSaved('EKF', 'MC_dynamics', 'MC State Dynamics')
    plt.show()

def create_chiErrorSum(N_MC, P_groundTruth, Xm_groundTruth, X0_est, k_f, dt, k_meas, omega_true):
    print 'N_MC = ', N_MC
    np.random.seed(100)
    X_truthModel_samp = np.random.multivariate_normal(Xm_groundTruth, P_groundTruth, size=N_MC)
    (r1, r2) = kf.chi2inv(N_MC)
    E_x_sum = np.zeros(k_f)
    for i in range(N_MC):
        X0_ground = X_truthModel_samp[i]
        X0_ground[0:3] = rbk.sigmaUnitSwitchCheck(X0_ground[0:3])
        print 'i_MC = ', i
        print 'X0_ground = ', X0_ground
        (k_vec, X_std_vec, X_error_vec, X_ground_vec, Y_meas_vec, omega_obs_vec, E_x_vec, E_y_ergo_vec) = \
            ekf.create_EKF_estimate(X0_ground, X0_est, k_f, dt, k_meas, omega_true)
        #ekf.plot_EKF_estimates(k_vec, X_std_vec, X_error_vec, dt)
        #plt.show()
        E_x_sum += E_x_vec
    E_x = (1.0/N_MC) * E_x_sum
    kf.plotChiResults(k_vec, r1, r2, E_x, N_MC, dt)
    kf.checkArePlotsSaved('EKF', 'chi_test_N_'+ str(N_MC), 'NEES Estimation Results')
    kf.checkArePlotsSaved('EKF', 'MC_dynamics', 'MC State Dynamics')
    plt.show()

if __name__ == "__main__":
    (k_f, dt, k_meas) = ekf.define_numericData()
    X0_est = ekf.define_initialGuess_KF()
    #X0_est = ekf.define_wrongInitialStateGuess_KF()
    #X0_est = ekf.define_largeInitialCovar_KF()
    (sigma0_true, omega_true, bias_true) = ekf.initialize_groundTruth()
    print 'sigma0_true = ', sigma0_true
    print 'omega_true = ', omega_true
    print 'bias_true = ', bias_true

    Xm_groundTruth = np.append(sigma0_true, bias_true)
    P_groundTruth = np.identity(6)
    P_groundTruth[0:3, 0:3] *= 0.0125
    P_groundTruth[3:6, 3:6] *= 0.#000005
    N_MC = 10


    create_chiErrorSum(N_MC, P_groundTruth, Xm_groundTruth, X0_est, k_f, dt, k_meas, omega_true)
    #create_chiRandom(P_groundTruth, Xm_groundTruth, X0_est, k_f, dt, k_meas, omega_true)





