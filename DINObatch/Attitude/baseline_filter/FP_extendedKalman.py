import numpy as np
from scipy import linalg as la
import KF_commonAlgs as kf
import matplotlib.pyplot as plt

import sys, os
userPath = os.path.expanduser('~')
sys.path.append(userPath + '/Desktop/astroPy/')
import RigidBodyKinematics as rbk

# ------------------------------------------ #
# This function plots error of mean state and std deviation 2*sigma bounds
def plotSigmaEstimate(t, m1_error, s1_vec, nfig, sigmaLim):
    sigmaLabel = 'MRP Attitude Set'
    plt.figure(nfig)
    kf.plotMeanErrorAndCovar(t, m1_error[0, :], s1_vec[0, :], kf.color_x, 'b', sigmaLim, sigmaLabel)
    kf.checkArePlotsSaved('EKF', 'sigma_plus_1', 'Estimate of state $\sigma_1$')
    plt.figure(nfig + 1)
    kf.plotMeanErrorAndCovar(t, m1_error[1, :], s1_vec[1, :], kf.color_y, 'g', sigmaLim, sigmaLabel)
    kf.checkArePlotsSaved('EKF', 'sigma_plus_2', 'Estimate of state $\sigma_2$')
    plt.figure(nfig + 2)
    kf.plotMeanErrorAndCovar(t, m1_error[2, :], s1_vec[2, :], kf.color_z, 'r', sigmaLim, sigmaLabel)
    kf.checkArePlotsSaved('EKF', 'sigma_plus_3', 'Estimate of state $\sigma_3$')

def plotOmegaBiasEstimate(t, m2_error, s2_vec, nfig, omegaLim):
    omegaLabel = 'Angular Rate [rad / s]'
    plt.figure(nfig + 3)
    kf.plotMeanErrorAndCovar(t, m2_error[0, :], s2_vec[0, :], kf.color_x, 'b', omegaLim, omegaLabel)
    kf.checkArePlotsSaved('EKF', 'bias_plus_1', 'Estimate of state $\omega_{b1}$')
    plt.figure(nfig + 4)
    kf.plotMeanErrorAndCovar(t, m2_error[1, :], s2_vec[1, :], kf.color_y, 'g', omegaLim, omegaLabel)
    kf.checkArePlotsSaved('EKF', 'bias_plus_2', 'Estimate of state $\omega_{b1}$')
    plt.figure(nfig + 5)
    kf.plotMeanErrorAndCovar(t, m2_error[2, :], s2_vec[2, :], kf.color_z, 'r', omegaLim, omegaLabel)
    kf.checkArePlotsSaved('EKF', 'bias_plus_3', 'Estimate of state $\omega_{b1}$')
# --------------------------------------------
def plot_groundTruth(X_ground_vec, k_vec, dt):
    sigma_lim = 1
    kf.plotSigma(k_vec, X_ground_vec[:3, :], 0, dt, sigma_lim)
    kf.checkArePlotsSaved('EKF', 'sigma_true', 'True MRP')
    omega_lim = -1
    kf.plotOmega(k_vec, X_ground_vec[3:, :], 1, dt, omega_lim)
    kf.checkArePlotsSaved('EKF', 'bias_true', 'True Angular Bias')
def plot_STmeas(Y_meas_vec, k_vec, dt):
    sigma_lim = 1
    kf.plotSigma(k_vec, Y_meas_vec, 2, dt, sigma_lim)
    kf.checkArePlotsSaved('EKF', 'sigma_meas', 'Measured MRP')

def plot_gyroObs(omega_obs_vec, k_vec, dt):
    omega_lim = -1
    kf.plotOmega(k_vec, omega_obs_vec, 3, dt, omega_lim)
    kf.checkArePlotsSaved('EKF', 'gyro_obs', 'Observed Angular Rates')

def plot_stateMeanPropagationError(k_vec, X_error_vec):
    sigma_lim = 0.01
    kf.plotSigma(k_vec, X_error_vec[:3, :], 4, dt, sigma_lim)
    kf.checkArePlotsSaved('EKF', 'sigma_propagation_error', 'MRP State Propagation Error')
    omega_lim = -1
    kf.plotOmega(k_vec, X_error_vec[3:, :], 5, dt, omega_lim)
    kf.checkArePlotsSaved('EKF', 'bias_propagation_error', 'Bias State Propagation Error')

def plot_statePropagation(k_vec, X_state_vec, X_std_vec, X_error_vec, dt):
    plot_stateMeanPropagationError(k_vec, X_error_vec)
    kf.plot_statePropagation(k_vec, X_state_vec[:3, :], X_std_vec[:3, :], X_state_vec[3:, :], X_state_vec[3:, :], dt, 6)

def plot_EKF_estimates(k_vec, X_std_vec, X_error_vec, dt):
    t = k_vec * dt
    sigmaLim = 0.05 #* 10.0
    plotSigmaEstimate(t, X_error_vec[0:3, :], X_std_vec[0:3, :], 12, sigmaLim)
    omegaLim = 0.001 #* 10.0 * 10.0
    plotOmegaBiasEstimate(t, X_error_vec[3:6, :], X_std_vec[3:6, :], 15, omegaLim)

def plot_chiSquare_NEES(k_vec, E_x_vec, dt):
    (r1, r2) = kf.chi2inv(1)
    kf.plotChiResults(k_vec, r1, r2, E_x_vec, 21, dt)
    kf.checkArePlotsSaved('EKF', 'NEES_test', 'NIS Estimation Results')

def plot_prefitResiduals(k_meas, E_y_ergo_vec, k_vec, dt):
    def chi2inv_NIS():
        # The values r1_y and r2_Y have been computed in Matlab
        # for k = 2000, p = 3
        r1_y = 2.8485
        r2_y = 3.1560
        return (r1_y, r2_y)

    (r1_y, r2_y) = chi2inv_NIS()
    idx_y = k_meas -1
    E_yy_ergo = E_y_ergo_vec[idx_y::idx_y]
    t_yy = dt * k_vec[idx_y::idx_y]
    plt.figure(22)
    plt.plot(t_yy, E_yy_ergo, 'bo')
    plt.axhline(r1_y / k_meas, color='dodgerblue')
    plt.axhline(r2_y/ k_meas, color='dodgerblue')
    plt.xlabel('Time, sec')
    plt.ylabel('NIS statistic, $\\bar{\epsilon}_y$')
    plt.legend(['NIS$_K$','$r_1$', '$r_2$'])
    kf.checkArePlotsSaved('EKF', 'NIS_test', 'NIS Estimation Results')


# --------------------------------------------
def F_NL(X, omega_obs):
    dsigma = kf.sigmaDot(X[0:3], omega_obs - X[3:6])
    dX = np.full_like(X, 0.0)
    dX[0:3] = dsigma
    return dX

def G_NL(X, eta):
    dsigma = - kf.sigmaDot(X[0:3], eta[0:3])
    dbias = eta[3:6]
    dX = np.full_like(X, 0.0)
    dX[0:3] = dsigma
    dX[3:6] = dbias
    return dX

def propagate_dynNL(X, params):
    dX = F_NL(X, params[0]) + G_NL(X, params[1])
    return dX

def propagate_dynKF(X, params_KF):
    dX_extended = np.full_like(X, 0.0)

    omega_obs = params_KF[0]
    A_LS = params_KF[1]
    G_LS = params_KF[2]
    Q = params_KF[3]
    # Propagate State: dX = F(X)
    dX_state =  F_NL(X[0:6], omega_obs)
    # Propagate Covariance: differential Ricatti equation
    P = np.resize(X[6:], (6,6))
    dF = np.dot(A_LS, P) + np.dot(P, A_LS.T)
    dG = np.dot(G_LS.T, np.dot(Q, G_LS))
    dP_mat = dF + dG
    n = dP_mat.shape[0] * dP_mat.shape[1]
    dP = np.resize(dP_mat, n)
    dX_extended[0:6] = dX_state
    dX_extended[6:] = dP
    return dX_extended
# --------------------------------------------
def store_KF_estimates(X_KF, X_state_vec, X_std_vec, k):
    X_state_vec[:,k] = X_KF[0:6]
    P = np.resize(X_KF[6:], (6,6))
    X_std_vec[:, k] = np.sqrt(np.diag(P))
    return

def check_ShadowTransform(X_KF):
    sigma = X_KF[0:3]
    s = la.norm(sigma)
    if s > 1.0:
        Lambda = np.identity(6)
        s_dot = np.dot(sigma, sigma)
        Lambda[0:3, 0:3] = 2.0* np.power(s_dot, -4) * np.outer(sigma, sigma) - np.power(s_dot, -2) * np.identity(3)
        P = np.resize(X_KF[6:], (6, 6))
        P_shadow = np.dot(Lambda.T, np.dot(P, Lambda))
        X_KF[6:] = np.resize(P_shadow, 6*6)
        X_KF[0:3] = kf.sigmaUnitSwitch(sigma)
    return X_KF


def computeInnovationError(P_minus, H, R, innov):
    P_yy = np.dot(H, np.dot(P_minus, (H.T))) + R
    P_yy_inv = la.inv(P_yy)
    e_y = np.dot(innov, np.dot(P_yy_inv, innov))
    return e_y

def computeKalmanCorrection(X_KF_minus, y, H, R):
    m_minus = X_KF_minus[0:6]
    P_minus = np.resize(X_KF_minus[6:], (6,6))

    K = kf.compute_KalmanGain(P_minus, H, R)
    innov = kf.compute_innovationFactor(y, H, m_minus)
    m_plus = m_minus + np.dot(K, innov)

    KH = np.dot(K, H)
    I = np.identity(KH.shape[0])
    M = I - KH
    P_plus = np.dot(M, np.dot(P_minus, M.T)) + np.dot(K, np.dot(R, K.T))

    P_plus_vec = np.resize(P_plus, 6*6)
    X_KF_plus = np.append(m_plus, P_plus_vec)

    e_y = computeInnovationError(P_minus, H, R, innov)
    return (X_KF_plus, e_y)
# --------------------------------------------

def initialize_groundTruth():
    def define_LKF():
        sigma0_true = np.array([0.2, 0.6, 0.1])
        omega_true = np.array([0., 0., 0.]) * rbk.D2R
        return (sigma0_true, omega_true)
    def define_EKF():
        sigma0_true = np.array([0., 0., 0.])
        omega_true = np.array([0.8, 0.2, 0.4]) * rbk.D2R
        return (sigma0_true, omega_true)
    (sigma0_true, omega_true) = define_LKF()
    #(sigma0_true, omega_true) = define_EKF()
    sigma0_true = rbk.MRPswitch(sigma0_true,1.0)
    bias_true = np.array([0.0001, 0.0002, 0.0004])
    #bias_true = np.array([0.001, 0.002, 0.004])
    return(sigma0_true, omega_true, bias_true)

def define_numericData():
    k_f = 2000
    dt = 0.5
    k_meas = 10
    return(k_f, dt, k_meas)

def define_noiseIntensities():
    def define_increasedNoise():
        (cov_ARW, cov_RRW) = define_nominalNoise()
        factor = 1E1
        return (cov_ARW*factor, cov_RRW)

    def define_nominalNoise():
        H2SEC = 3600.0
        std_ARW = 0.55 # deg / h^(1/2)
        std_RRW = 7.0 # deg / h^(3/2)
        cov_ARW = std_ARW*std_ARW * np.power(rbk.D2R, 2) / H2SEC
        cov_RRW = std_RRW*std_RRW * np.power(rbk.D2R, 2) / np.power(H2SEC, 3)
        return (cov_ARW, cov_RRW)


    #(cov_ARW, cov_RRW) = define_nominalNoise()
    (cov_ARW, cov_RRW) = define_increasedNoise()
    return (cov_ARW, cov_RRW)
# --------------------------------------------
def define_initialGuess_KF():
    (sigma0_true, omega_true, bias_true) = initialize_groundTruth()
    X0 = np.append(sigma0_true, bias_true)

    P0_sigma = np.identity(3) * 0.175
    P0_bias = np.identity(3) * 0.005  # (rad/s)^2
    P0 = la.block_diag(P0_sigma, P0_bias)
    n = P0.shape[0] * P0.shape[1]
    P0_vec = np.resize(P0, n)

    X0_KF = np.append(X0, P0_vec)
    return X0_KF

def define_largeInitialCovar_KF():
    (sigma0_true, omega_true, bias_true) = initialize_groundTruth()
    X0 = np.append(sigma0_true, bias_true)
    P0 = np.identity(6) * 100
    n = P0.shape[0] * P0.shape[1]
    P0_vec = np.resize(P0, n)

    X0_KF = np.append(X0, P0_vec)
    return X0_KF

def define_wrongInitialStateGuess_KF():
    sigma0 = np.array([0.5, 0.5, 0.5])
    bias0 = np.zeros(3)
    X0 = np.append(sigma0, bias0)
    P0_sigma = np.identity(3) * 0.175
    P0_bias = np.identity(3) * 0.005  # (rad/s)^2
    P0 = la.block_diag(P0_sigma, P0_bias)
    n = P0.shape[0] * P0.shape[1]
    P0_vec = np.resize(P0, n)
    X0_KF = np.append(X0, P0_vec)
    return X0_KF

# --------------------------------------------
def define_PSD_Q():
    (cov_ARW, cov_RRW) = define_noiseIntensities()
    P_ARW = cov_ARW * np.identity(3)
    P_RRW = cov_RRW * np.identity(3)
    Q = la.block_diag(P_ARW, P_RRW)
    return Q

def define_covarMat_R():
    R = np.identity(3) * 0.0004
    return R

def create_measNoiseSamples(k_f):
    R = define_covarMat_R()
    np.random.seed(100)
    v_samp = np.random.multivariate_normal(np.zeros(3), R, size=k_f)
    return v_samp

def create_processNoiseSamples(k_f):
    Q = define_PSD_Q()
    np.random.seed(100)
    w_samp = np.random.multivariate_normal(np.zeros(6), Q, size=k_f)
    return w_samp

# --------------------------------------------
def create_noiseSamples(cov_ARW, cov_RRW, R, k_f):
    np.random.seed(100)
    P_ARW = cov_ARW * np.identity(3)
    P_RRW = cov_RRW * np.identity(3)
    P_w = la.block_diag(P_ARW, P_RRW)
    w_samp = np.random.multivariate_normal(np.zeros(6), P_w, size=k_f)
    v_samp = np.random.multivariate_normal(np.zeros(3), R, size=k_f)
    return (w_samp, v_samp)
# --------------------------------------------

def create_noisyMeasurements(X, H, v_samp, k):
    sigma_meas = np.dot(H, X) + v_samp[k, :]
    return sigma_meas



def create_gyroObservations(B0, omega_true, cov_ARW, cov_RRW, dt, k, omega_obs_vec):
    N_v = np.random.multivariate_normal(np.zeros(3), np.identity(3))
    N_u = np.random.multivariate_normal(np.zeros(3), np.identity(3))
    B = B0 + np.sqrt(cov_RRW*dt) * N_u
    disc = cov_ARW/dt + (1.0/12.0)*cov_RRW*dt
    omega_obs = omega_true + 0.5*(B + B0) + N_v * np.sqrt(disc)
    #print 'w_add = ', 0.5*(B + B0) + N_v * np.sqrt(disc)
    omega_obs_vec[:,k] = omega_obs
    return(omega_obs, B)

def create_groundTruth(X0, k_f, dt):
    k_vec = np.zeros(k_f)
    X_ground_vec = np.zeros([6, k_f])
    X_old = X0.copy()

    Y_meas_vec = np.zeros([3, k_f])
    H = kf.define_H_mat()
    v_samp = create_measNoiseSamples(k_f)
    w_samp = create_processNoiseSamples(k_f)

    omega_obs_vec = np.zeros([3, k_f])
    (cov_ARW, cov_RRW) = define_noiseIntensities()
    B0 = X0[3:6].copy()

    for k in range(0, k_f):
        k_vec[k] = k
        X_ground_vec[:, k] = X_old
        Y_meas_vec[:, k] = create_noisyMeasurements(X_old, H, v_samp, k)
        (omega_obs, B) = create_gyroObservations(B0, omega_true, cov_ARW, cov_RRW, dt, k, omega_obs_vec)
        params_NL = (omega_obs, w_samp[k, :])
        X = kf.RK4(X_old, propagate_dynNL, params_NL, dt)
        X[0:3] = rbk.sigmaUnitSwitchCheck(X[0:3])
        X_old = X
        B0 = B
    return (X_ground_vec, Y_meas_vec, omega_obs_vec, k_vec)


def create_monteCarlo(X0_ground, X_est, k_f, dt):
    # GROUND TRUTH
    (sigma0_true, omega_true, bias_true) = initialize_groundTruth()
    k_vec = np.zeros(k_f)
    X_ground_vec = np.zeros([6, k_f])
    X0_GT = X0_ground.copy()

    Y_meas_vec = np.zeros([3, k_f])
    H = kf.define_H_mat()
    v_samp = create_measNoiseSamples(k_f)
    w_samp = create_processNoiseSamples(k_f)

    omega_obs_vec = np.zeros([3, k_f])
    (cov_ARW, cov_RRW) = define_noiseIntensities()
    B0 = X0_ground[3:6].copy()

    # KF ESTIMATES
    X_state_vec = np.zeros([6, k_f])
    X_std_vec = np.zeros([6, k_f])
    X0_KF = X_est.copy()
    Q = define_PSD_Q()

    # MEAN ERROR
    X_error_vec = np.zeros([6, k_f])

    for k in range(0, k_f):
        # GROUND TRUTH
        k_vec[k] = k
        X_ground_vec[:, k] = X0_GT
        Y_meas_vec[:, k] = create_noisyMeasurements(X0_GT, H, v_samp, k)
        (omega_obs, B) = create_gyroObservations(B0, omega_true, cov_ARW, cov_RRW, dt, k, omega_obs_vec)
        params_NL = (omega_obs, w_samp[k, :])
        # Propagate Ground Truth
        X_GT = kf.RK4(X0_GT, propagate_dynNL, params_NL, dt)
        X_GT[0:3] = rbk.MRPswitch(X_GT[0:3],1)
        # Update Current Ground Truth
        X0_GT = X_GT
        B0 = B

        # KF PREDICTION
        store_KF_estimates(X0_KF, X_state_vec, X_std_vec, k)
        omega_ref = omega_obs - X0_KF[3:6]
        (A_LS, G_LS) = kf.linearizeSystem(X0_KF[0:3], omega_ref)
        # Propagate KF State Estimate:
        params_KF = (omega_obs, A_LS, G_LS, Q)
        X_KF = kf.RK4(X0_KF, propagate_dynKF, params_KF, dt)
        check_ShadowTransform(X_KF)
        # Update KF State Estimate:
        X0_KF = X_KF

        # COMPUTE ERROR
        sigma_error = rbk.MRPswitch(X_GT[0:3] - X_KF[0:3],1.0)
        bias_error = X_GT[3:6] - X_KF[3:6]
        X_error_vec[:,k] = np.append(sigma_error, bias_error)


    return(k_vec, X_state_vec, X_std_vec, X_error_vec)

def Q_sigma_corrupted():
    Q = define_PSD_Q()
    Q_corr = Q.copy()
    Q_corr[0:3, 0:3] = np.identity(3) * 1E-4
    return Q_corr

def Q_bias_corrupted():
    Q = define_PSD_Q()
    Q_corr = Q.copy()
    Q_corr[3:6, 3:6] = np.identity(3) * 1E-4
    return Q_corr

def Q_corrupted():
    Q = define_PSD_Q()
    Q_corr = Q.copy()
    Q_corr *= 1E-2
    #Q_corr = np.zeros(6)
    #Q_corr = np.identity(6) * 1E-4
    return Q_corr

def compute_chiSquareError(e_xk, X_KF):
    P_plus_k = np.resize(X_KF[6:], (6,6))
    E_xk = np.dot(e_xk, np.dot(la.inv(P_plus_k), e_xk))
    return E_xk

def create_EKF_estimate(X0_ground, X_est, k_f, dt, k_meas, omega_true):
    # GROUND TRUTH
    k_vec = np.zeros(k_f)
    X_ground_vec = np.zeros([6, k_f])
    X0_GT = X0_ground.copy()

    Y_meas_vec = np.zeros([3, k_f])
    H = kf.define_H_mat()
    v_samp = create_measNoiseSamples(k_f)
    w_samp = create_processNoiseSamples(k_f)

    omega_obs_vec = np.zeros([3, k_f])
    (cov_ARW, cov_RRW) = define_noiseIntensities()
    B0 = X0_ground[3:6].copy()

    # KF ESTIMATES
    X_state_vec = np.zeros([6, k_f])
    X_std_vec = np.zeros([6, k_f])
    X0_KF = X_est.copy()
    #Q = define_PSD_Q() #Q = Q_sigma_corrupted() #Q = Q_bias_corrupted()
    Q = Q_corrupted()
    R = define_covarMat_R()

    # MEAN ERROR
    X_error_vec = np.zeros([6, k_f])

    # CHI SQUARE TEST
    E_x_vec = np.zeros(k_f)
    E_y_vec = np.zeros(k_f)
    E_y_ergo_vec = np.zeros(k_f)

    k_counter = 0

    #(sigma0_true, omega_true, bias_true) = initialize_groundTruth()
    for k in range(0, k_f):
        # GROUND TRUTH
        k_vec[k] = k
        X_ground_vec[:, k] = X0_GT
        Y_meas_vec[:, k] = create_noisyMeasurements(X0_GT, H, v_samp, k)
        (omega_obs, B) = create_gyroObservations(B0, omega_true, cov_ARW, cov_RRW, dt, k, omega_obs_vec)
        params_NL = (omega_obs, w_samp[k, :])
        # Propagate Ground Truth
        X_GT = kf.RK4(X0_GT, propagate_dynNL, params_NL, dt)
        X_GT[0:3] = rbk.MRPSwitch(X_GT[0:3],1.0)
        # Update Current Ground Truth
        X0_GT = X_GT
        B0 = B

        # KF PREDICTION
        store_KF_estimates(X0_KF, X_state_vec, X_std_vec, k)
        omega_ref = omega_obs - X0_KF[3:6]
        (A_LS, G_LS) = kf.linearizeSystem(X0_KF[0:3], omega_ref)
        (A_LS, G_LS) = kf.linearizeSystem(X0_ground[0:3], omega_true)
        # Propagate KF State Estimate:
        params_KF = (omega_obs, A_LS, G_LS, Q)
        X_KF_minus = kf.RK4(X0_KF, propagate_dynKF, params_KF, dt)
        check_ShadowTransform(X_KF_minus)

        # KF CORRECTION
        if k_counter == (k_meas - 1):
            (X_KF_plus, e_y) = computeKalmanCorrection(X_KF_minus, Y_meas_vec[:, k], H, R)
            E_y_vec[k] = e_y
            E_y_ergo_vec[k] = float(1.0/(k+1)) * np.sum(E_y_vec[:k])

            check_ShadowTransform(X_KF_plus)
            X0_KF = X_KF_plus
            k_counter = 0

        else:
            X0_KF = X_KF_minus

        k_counter += 1

        # COMPUTE ERROR
        sigma_error = rbk.MRPSwitch(X_GT[0:3] - X0_KF[0:3],1.0)
        bias_error = X_GT[3:6] - X0_KF[3:6]
        e_x = np.append(sigma_error, bias_error)
        X_error_vec[:,k] = e_x

        # CHI SQUARE TEST
        E_xk = compute_chiSquareError(e_x, X0_KF)
        E_x_vec[k] = E_xk
    return (k_vec, X_std_vec, X_error_vec, X_ground_vec, Y_meas_vec, omega_obs_vec, E_x_vec, E_y_ergo_vec)


def run_groundTruthSim(X0, k_f, dt):
    (X_ground_vec, Y_meas_vec, omega_obs_vec, k_vec) = create_groundTruth(X0, k_f, dt)
    plot_groundTruth(X_ground_vec, k_vec, dt)
    plot_STmeas(Y_meas_vec, k_vec, dt)
    plot_gyroObs(omega_obs_vec, k_vec, dt)
    plt.show()

def run_monteCarloSim(X0_ground, X0_est, k_f, dt):
    (k_vec, X_state_vec, X_std_vec, X_error_vec) = create_monteCarlo(X0_ground, X0_est, k_f, dt)
    plot_statePropagation(k_vec, X_state_vec, X_std_vec, X_error_vec, dt)
    plt.show()


def run_extended_KF(X0_ground, X0_est, k_f, dt, k_meas, omega_true):
    (k_vec, X_std_vec, X_error_vec, X_ground_vec, Y_meas_vec, omega_obs_vec, E_x_vec, E_y_ergo_vec) = \
        create_EKF_estimate(X0_ground, X0_est, k_f, dt, k_meas, omega_true)
    #plot_groundTruth(X_ground_vec, k_vec, dt)
    #plot_STmeas(Y_meas_vec, k_vec, dt)
    #plot_gyroObs(omega_obs_vec, k_vec, dt)
    plot_EKF_estimates(k_vec, X_std_vec, X_error_vec, dt)
    plot_prefitResiduals(k_meas, E_y_ergo_vec, k_vec, dt)
    plot_chiSquare_NEES(k_vec, E_x_vec, dt)
    plt.show()



if __name__ == "__main__":
    (k_f, dt, k_meas) = define_numericData()
    print 't_sim = ', k_f*dt

    (cov_ARW, cov_RRW) = define_noiseIntensities()
    print 'cov_ARW = ', cov_ARW
    print 'cov_RRW = ', cov_RRW

    (sigma0_true, omega_true, bias_true) = initialize_groundTruth()
    X0_ground = np.append(sigma0_true, bias_true)
    print 'sigma0_true = ', sigma0_true
    print 'omega_true = ', omega_true
    print 'bias_true = ', bias_true

    #run_groundTruthSim(X0_ground, k_f, dt)
    X0_est = define_initialGuess_KF()
    X0_est = define_wrongInitialStateGuess_KF() #X0_est = define_largeInitialCovar_KF()
    print 'X0_est = ', X0_est[0:6]
    run_monteCarloSim(X0_ground, X0_est, k_f, dt)
    run_extended_KF(X0_ground, X0_est, k_f, dt, k_meas, omega_true)

    #(A_LS, G_LS) = kf.linearizeSystem(sigma0_true, omega_true)
    W = define_PSD_Q()
    (F_LS, Q_LS) = kf.compute_DTSystem(G_LS, W, A_LS, dt)
    kf.printM6(A_LS, 'A_LS')
    kf.printM6(Q_LS, 'Q_LS')
    kf.printM6(F_LS, 'F_LS')








