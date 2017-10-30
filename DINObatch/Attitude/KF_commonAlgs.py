import matplotlib.pyplot as plt
import numpy as np
from scipy import linalg as la

import sys, os
userPath = os.path.expanduser('~')
sys.path.append(userPath + '/Desktop/astroPy/')
import RigidBodyKinematics as rbk


#  ------------------- PLOTS DIRECTORY ------------------- #
# Optionally specify the full path of the directory where you want figures to be saved
#plotsPath = '/Users/marcolsmargenet/Desktop/FALL_TERM_2016/5044/Final_Project/FiguresEKF/'
plotsPath = '/Users/marcolsmargenet/Desktop/FALL_TERM_2016/5044/Final_Report/FiguresMC/'
arePlotsSaved = True
TheSim = 'EKF'
if arePlotsSaved:
    factor = 1#.5
    plt.rcParams['figure.figsize'] = 1.75*factor, 1.5*factor #3.5, 3.5 #
    plt.rcParams.update({'font.size': 6}) # 8
else:
    plt.rcParams['figure.figsize'] = 5.2, 4. #1.75, 1.5
    plt.rcParams.update({'font.size': 11})

color_x = 'dodgerblue'
color_y = 'limegreen'
color_z = 'salmon'

# ------------------- COMMON METHODS ------------------- #
def printM6(M, M_str):
    print M_str + '_11 = \n', M[:3, :3]
    print M_str + '_12 = \n', M[:3, 3:]
    print M_str + '_22 = \n', M[3:, 3:]
    print M_str + '_21 = \n', M[3:, :3]
    return

def printEigenM(M, M_str):
    (w_A, v_A) = la.eig(M)
    print 'M = ', M_str
    for i in range(0, len(w_A)):
        print 'w_M = ', w_A[i]
        print 'v_M = \n', v_A[:,i]
    return

def printEigen(A, F):
    print 'Pole Locations:'
    (w_A, v_A) = la.eig(A)
    (w_F, v_F) = la.eig(F)
    for i in range(0, len(w_A)):
        print 'w_A = ', w_A[i]
        print 'v_A = \n', v_A[:,i]
    for i in range(0, len(w_F)):
        print 'w_F = ', w_F[i]
        print 'v_F = \n', v_F[:,i]
    return

def checkArePlotsSaved(TheSim, pltName, pltTitle):
    if arePlotsSaved:
        plt.savefig(plotsPath + pltName +".pdf", bbox_inches='tight')
    else:
        plt.title(TheSim + ': ' + pltTitle)

# ------------------------------------------ #

def chi2inv(N_MC):
    # The values of r1 and r2 have been computed in Matlab
    # For n = 6 and alpha = 0.005
    if N_MC == 1:
        r1 = 0.5266
        r2 =  20.2494
    elif N_MC == 5:
        r1 = 2.5529
        r2 = 11.2665
    elif N_MC == 10:
        r1 = 3.3791
        r2 = 9.5344
    elif N_MC == 15:
        r1 = 3.7928
        r2 = 8.8170
    elif N_MC == 25:
        r1 = 4.2377
        r2 = 8.1286
    elif N_MC == 50:
        r1 = 4.7163
        r2 = 7.4670
    else:
        raise ValueError('Chi 2 inv bounds not defined for N_MC = ' + str(N_MC))
    return (r1, r2)

def plotChiResults(k_vec, r1, r2, E_x_vec, nfig, dt):
    t_vec = k_vec * dt
    plt.figure(nfig)
    plt.plot(t_vec, E_x_vec, 'ro')
    plt.axhline(r1, color='lightsalmon')
    plt.axhline(r2, color='lightsalmon')
    lim_offset = 2
    plt.ylim([ int(r1) - lim_offset*2 - 4, int(r2) + lim_offset*2 + 4])
    plt.xlabel('Time, sec')
    plt.ylabel('NEES statistic, $\\bar{\epsilon}_x$')
    plt.legend(['NEES$_K$','$r_1$', '$r_2$'])


# ------------------------------------------ #
# This function plots error of mean state and std deviation 2*sigma bounds
def plotSigmaEstimate(t, m1_error, s1_vec, nfig):
    sigmaLim = 0.005
    sigmaLabel = 'MRP Attitude Set'
    plt.figure(nfig)
    plotMeanErrorAndCovar(t, m1_error[0, :], s1_vec[0, :], color_x, 'b', sigmaLim, sigmaLabel)
    checkArePlotsSaved(TheSim, 'sigma_plus_1', 'Estimate of state $\sigma_1$')
    plt.figure(nfig + 1)
    plotMeanErrorAndCovar(t, m1_error[1, :], s1_vec[1, :], color_y, 'g', sigmaLim, sigmaLabel)
    checkArePlotsSaved(TheSim, 'sigma_plus_2', 'Estimate of state $\sigma_2$')
    plt.figure(nfig + 2)
    plotMeanErrorAndCovar(t, m1_error[2, :], s1_vec[2, :], color_z, 'r', sigmaLim, sigmaLabel)
    checkArePlotsSaved(TheSim, 'sigma_plus_3', 'Estimate of state $\sigma_3$')

def plotOmegaBiasEstimate(t, m2_error, s2_vec, nfig):
    omegaLim = -1
    omegaLabel = 'Angular Rate [rad / s]'
    plt.figure(nfig + 3)
    plotMeanErrorAndCovar(t, m2_error[0, :], s2_vec[0, :], color_x, 'b', omegaLim, omegaLabel)
    checkArePlotsSaved(TheSim, 'bias_plus_1', 'Estimate of state $\omega_{b1}$')
    plt.figure(nfig + 4)
    plotMeanErrorAndCovar(t, m2_error[1, :], s2_vec[1, :], color_y, 'g', omegaLim, omegaLabel)
    checkArePlotsSaved(TheSim, 'bias_plus_2', 'Estimate of state $\omega_{b1}$')
    plt.figure(nfig + 5)
    plotMeanErrorAndCovar(t, m2_error[2, :], s2_vec[2, :], color_z, 'r', omegaLim, omegaLabel)
    checkArePlotsSaved(TheSim, 'bias_plus_3', 'Estimate of state $\omega_{b1}$')

def plot_KalmanStateEstimates(k_vec, m1_vec, sigma_true, s1_vec, m2_vec, bias_true, s2_vec, dt, nfig):
    m1_error = m1_vec - sigma_true
    m2_error = m2_vec - bias_true
    t = k_vec * dt
    plotSigmaEstimate(t, m1_error, s1_vec, nfig)
    plotOmegaBiasEstimate(t, m2_error, s2_vec, nfig)

    return


# ------------------------------------------ #
def plotMeanAndCovar(t, m, s, c_m, c_s, lim, yLabel):
    upperBound = m + 2 * s
    lowerBound = m - 2 * s
    plt.plot(t, m, c_m)
    plt.plot(t, upperBound, c_s, marker='_')
    plt.plot(t, lowerBound, c_s, marker='_')
    plt.legend(['Mean $\mu$', '2$\sigma$ bounds'])
    plt.xlabel('Time, sec')
    plt.ylabel(yLabel)
    if lim > 0:
        plt.ylim([-lim, lim])

def plotMeanErrorAndCovar(t, m_error, s, c_m, c_s, lim, yLabel):
    upperBound = 2 * s
    lowerBound = -2 * s
    plt.plot(t, m_error, c_m)
    plt.plot(t, upperBound, c_s, marker='_')
    plt.plot(t, lowerBound, c_s, marker='_')
    plt.legend(['Mean $\mu$', '2$\sigma$ bounds'])
    plt.xlabel('Time, sec')
    plt.ylabel(yLabel)
    if lim > 0:
        plt.ylim([-lim, lim])


def plot_statePropagation(k_vec, m1_vec, s1_vec, m2_vec, s2_vec, dt, nfig):
    def plotSigmaState(t, m1_vec, s1_vec, nfig):
        sigmaLim = -0.1
        sigmaLabel = 'MRP Attitude Set'
        plt.figure(nfig)
        plotMeanAndCovar(t, m1_vec[0, :], s1_vec[0, :], color_x, 'b', sigmaLim, sigmaLabel)
        checkArePlotsSaved(TheSim, 'sigma_minus_1', 'Propagation of state $\sigma_1$')
        plt.figure(nfig+1)
        plotMeanAndCovar(t, m1_vec[1, :], s1_vec[1, :], color_y, 'g', sigmaLim, sigmaLabel)
        checkArePlotsSaved(TheSim, 'sigma_minus_2', 'Propagation of state $\sigma_2$')
        plt.figure(nfig+2)
        plotMeanAndCovar(t, m1_vec[2, :], s1_vec[2, :], color_z, 'r', sigmaLim, sigmaLabel)
        checkArePlotsSaved(TheSim, 'sigma_minus_3', 'Propagation of state $\sigma_3$')

    def plotOmegaBiasState(t, m2_vec, s2_vec, nfig):
        omegaLim = -s2_vec[0, 0] * 3
        omegaLim = -0.01
        omegaLabel = 'Angular Rate [rad / s]'
        plt.figure(nfig+3)
        plotMeanAndCovar(t, m2_vec[0, :], s2_vec[0, :], color_x, 'b', omegaLim, omegaLabel)
        checkArePlotsSaved(TheSim, 'bias_minus_1', 'Propagation of state $\omega_{b1}$')
        plt.figure(nfig+4)
        plotMeanAndCovar(t, m2_vec[1, :], s2_vec[1, :], color_y, 'g', omegaLim, omegaLabel)
        checkArePlotsSaved(TheSim, 'bias_minus_2', 'Propagation of state $\omega_{b1}$')
        plt.figure(nfig+5)
        plotMeanAndCovar(t, m2_vec[2, :], s2_vec[2, :], color_z, 'r', omegaLim, omegaLabel)
        checkArePlotsSaved(TheSim, 'bias_minus_3', 'Propagation of state $\omega_{b1}$')

    t = k_vec * dt
    plotSigmaState(t, m1_vec, s1_vec, nfig)
    plotOmegaBiasState(t, m2_vec, s2_vec, nfig)
# ------------------------------------------ #

def plot_meanErrors(k_vec, m1_vec, sigma_true_vec, m2_vec, omegaBias_true_vec, dt):
    def plotSigmaError(k_vec, m1_vec, sigma_true_vec, dt):
        sigmaLim = 1
        m1_error = m1_vec - sigma_true_vec
        for i in range(0, m1_vec.shape[1]):
            m1_error[:, i] = rbk.sigmaUnitSwitchCheck(m1_error[:, i])
        plotSigma(k_vec, m1_error, 9, dt,  sigmaLim)
        checkArePlotsSaved(TheSim, 'error_sigma', 'Error for state $\sigma$')

    def plotOmegaBiasError(k_vec, m2_vec, omegaBias_true_vec, dt):
        omegaLim = -1
        m2_error = m2_vec - omegaBias_true_vec
        plotOmega(k_vec, m2_error, 10, dt,  omegaLim)
        checkArePlotsSaved(TheSim, 'error_bias', 'Error for state $\omega_{b}$')

    plotSigmaError(k_vec, m1_vec, sigma_true_vec, dt)
    #plotOmegaBiasError(t, m2_error)

# ------------------------------------------ #

def plot_meanPropagationOnly(k_vec, m1_vec, sigma_true_vec, m2_vec, omegaBias_true_vec, dt):
    def plotSigmaMean(k_vec, m1_vec, dt):
        sigmaLim = 1
        plotSigma(k_vec, m1_vec, 11, dt,  sigmaLim)
        checkArePlotsSaved(TheSim, 'sigma_minus', 'Mean Propagation for state $\sigma$')

    def plotOmegaBiasMeans(k_vec, m2_vec, dt):
        omegaLim = -1
        plotOmega(k_vec, m2_vec, 12, dt,  omegaLim)
        checkArePlotsSaved(TheSim, 'bias_minus', 'Mean Propagation for state $\omega_{b}$')

    plotSigmaMean(k_vec, m1_vec, dt)
    plotOmegaBiasMeans(k_vec, m2_vec, dt)

# ------------------------------------------ #

def plot_noisyMeasurements(k_vec, u_vec, y_vec, dt):
    sigmaLim = 1
    plotSigma(k_vec, y_vec, 13, dt,  sigmaLim)
    checkArePlotsSaved(TheSim, 'sigma_meas', 'Star Tracker MRP Measurements $\\tilde{\sigma}$')

    omegaLim = -1
    plotOmega(k_vec, u_vec, 14, dt,  omegaLim)
    checkArePlotsSaved(TheSim, 'omega_meas', 'Gyro Input Rate $\\tilde{\omega}$')

# ------------------------------------------ #
def plotCovarBounds(t, sigma_vec, std_vec, nfig):
    plt.figure(nfig)
    plt.plot(t, sigma_vec[0, :] + 2*std_vec[0, :], color_x, marker='_')
    plt.plot(t, sigma_vec[0, :] - 2*std_vec[0, :], color_x, marker='_')
    plt.plot(t, sigma_vec[1, :] + 2*std_vec[1, :], color_y, marker='_')
    plt.plot(t, sigma_vec[1, :] - 2*std_vec[1, :], color_y, marker='_')
    plt.plot(t, sigma_vec[2, :] + 2*std_vec[2, :], color_z, marker='_')
    plt.plot(t, sigma_vec[2, :] - 2*std_vec[2, :], color_z, marker='_')

def plotSigma_meanAndCovar(k_vec, sigma_vec, std_vec, nfig, dt, sigma_lim):
    t = k_vec * dt
    plt.figure(nfig)
    plt.plot(t, sigma_vec[0, :], color_x)
    plt.plot(t, sigma_vec[1, :], color_y)
    plt.plot(t, sigma_vec[2, :], color_z)
    plotCovarBounds(t, sigma_vec, std_vec, nfig)
    plt.legend(['$\sigma_1$', '$\sigma_2$', '$\sigma_3$'])
    plt.xlabel('Time, sec')
    plt.ylabel('MRP Attitude Set')
    if sigma_lim > 0:
        plt.ylim([-sigma_lim, sigma_lim])

# ------------------------------------------ #
def plotBiasObs(k_vec, omega_i_obs, omega_true, i, nfig, dt, omega_lim):
    plt.figure(nfig)
    omega_i_true = np.full_like(omega_i_obs, omega_true[i])
    t_vec = k_vec * dt
    plt.plot(t_vec, omega_i_true, 'b')
    plt.plot(t_vec, omega_i_obs, color_x)
    plt.legend(['$\omega_{true}$','$\omega_{meas}$'])
    plt.xlabel('Time, sec')
    plt.ylabel('Angular Rate [rad/s]')
    if omega_lim > 0:
        plt.ylim([-omega_lim, omega_lim])

def plotSigmaMeas(k_vec, k_meas, sigma_i_meas, nfig, dt, sigma_lim):
    plt.figure(nfig)
    t_vec = k_vec * dt
    plt.plot(t_vec, sigma_i_meas, color_x)
    #plt.scatter(t_vec[0::k_meas], sigma_i_meas[0::k_meas], color='b')
    #plt.plot(t_vec[0::k_meas], sigma_i_meas[0::k_meas], 'bo')
    plt.plot(t_vec[0::k_meas], sigma_i_meas[0::k_meas], c='b', marker='.')
    plt.legend(['$\sigma$','$\sigma_{meas}$'])
    plt.xlabel('Time, sec')
    plt.ylabel('MRP Attitude Set')
    if sigma_lim > 0:
        plt.ylim([-sigma_lim, sigma_lim])
# ------------------------------------------ #

# This function plots MRP states
def plotSigma(k_vec, sigma_vec, nfig, dt, sigma_lim):
    t = k_vec * dt
    plt.figure(nfig)
    plt.plot(t, sigma_vec[0, :], color_x)
    plt.plot(t, sigma_vec[1, :], color_y)
    plt.plot(t, sigma_vec[2, :], color_z)
    plt.legend(['$\sigma_1$', '$\sigma_2$', '$\sigma_3$'])
    plt.xlabel('Time, sec')
    plt.ylabel('MRP Attitude Set')
    if sigma_lim > 0:
        plt.ylim([-sigma_lim, sigma_lim])

# This function plots angular velocity states
def plotOmega(k_vec, omega_vec, nfig, dt, omega_lim):
    t = k_vec * dt
    plt.figure(nfig)
    plt.plot(t, omega_vec[0, :], color_x)
    plt.plot(t, omega_vec[1, :], color_y)
    plt.plot(t, omega_vec[2, :], color_z)
    plt.legend(['$\omega_1$', '$\omega_2$', '$\omega_3$'])
    plt.xlabel('Time, sec')
    plt.ylabel('Angular velocity, rad/s')
    if omega_lim > 0:
        plt.ylim([-omega_lim, omega_lim])

def plot_groundTruth(k_vec, sigma_vec, omegaBias_vec, omega_vec, dt):
    nfig = 0
    sigma_lim = 1
    plotSigma(k_vec, sigma_vec, nfig, dt, sigma_lim)
    checkArePlotsSaved(TheSim, 'sigma_true', 'True MRP')

    nfig = 1
    omega_lim = -1
    plotOmega(k_vec, omegaBias_vec, nfig, dt, omega_lim)
    checkArePlotsSaved(TheSim, 'omegaBias_true', 'True Angular Bias')

    nfig = 2
    omega_lim = -1
    plotOmega(k_vec, omega_vec, nfig, dt, omega_lim)
    checkArePlotsSaved(TheSim, 'omega_true', 'True Angular Velocity')

def storeKalmanOutputVars(X_m0, P0, m1_vec, m2_vec, s1_vec, s2_vec, k_vec, k):
    m1_vec[:, k] = X_m0[0:3]
    m2_vec[:, k] = X_m0[3:6]
    std = np.diag(P0)
    s1_vec[:, k] = std[0:3]
    s2_vec[:, k] = std[3:6]
    k_vec[k] = k
    return
# ------------------------------------------ #
# This is a Runge-Kutta 4 integration scheme
def RK4(y, g, params, dt):
    k1 = g(y, params) * dt
    k2 = g(y + 0.5 * k1, params) * dt
    k3 = g(y + 0.5 * k2, params) * dt
    k4 = g(y + k3, params) * dt
    y_new = y + 1./6. * (k1 + 2*k2 + 2*k3 + k4)
    return y_new
# ------------------------------------------ #
# This function returns the time-derivative of sigma (MRP set)
def sigmaDot(sigma, omega): # dSigma = 1 / 4. * B[sigma] * omega
    B_mat = rbk.BmatMRP(sigma)
    dsigma = 1 / 4. * np.dot(B_mat, omega)
    return dsigma
# ------------------------------------------ #
# These functions compute the time derivatives of the state X = [sigma, omega_bias].T
# return: dX = [dsigma, domega_bias]

# Reference State: omega = omega_true, X[3:6] = omega_bias_true
def F_ref(X, params):
    sigma = X[0:3]
    omega = params # true rate
    dsigma = sigmaDot(sigma, omega)
    dX = np.full_like(X, 0.0)
    dX[0:3] = dsigma
    dX[3:6] = np.array([0.0, 0.0, 0.0]) # true bias is constant
    return dX

# State propagated with omega = omega_meas - X[3:6] - eta_w
def F_vec(X, params):
    sigma = X[0:3]
    omega_meas = params[0] # gyro measurement
    eta_w = params[1] # omega_noise
    eta_wb = params[2] # bias_noise
    omega = omega_meas - X[3:6] - eta_w
    dsigma = sigmaDot(sigma, omega)
    domega_bias = eta_wb
    dX = np.full_like(X, 0.0)
    dX[0:3] = dsigma
    dX[3:6] = domega_bias
    return dX
# ------------------------------------------ #

def linearizeSystem(sigma_ref, omega_ref):
    # Compute F11 submatrix
    swT = np.outer(sigma_ref, omega_ref)
    wsT =  np.outer(omega_ref, sigma_ref)
    w_tilde = rbk.tildeMatrix(omega_ref)
    swI = np.dot(sigma_ref, omega_ref) * np.identity(3)
    F_11 = 0.5 * (swT - wsT - w_tilde + swI)
    # Compute F12 submatrix
    B = rbk.BmatMRP(sigma_ref)
    F_12 = -0.25 * B
    # Compute F matrix
    F = np.zeros([6, 6])
    F[0:3, 0:3] = F_11
    F[0:3, 3:6] = F_12
    # Compute G matrix
    G = np.identity(6)
    G[0:3, 0:3] = -0.25 * B
    return (F, G)


def compute_DTSystem(Tau, W, A, dt):
    TWT = np.dot(Tau, np.dot(W, Tau.T))
    n = A.shape[0]
    M = np.zeros([2 * n, 2 * n])
    M[0:n, 0:n] = -A
    M[n: 2 * n, n: 2 * n] = A.T
    M[0:n, n:2 * n] = TWT
    M = dt * M
    e_M = la.expm(M)
    F_T = e_M[n: 2 * n, n: 2 * n]
    FinvQ = e_M[0:n, n:2 * n]
    Q = np.dot(F_T.T, FinvQ)
    F = F_T.T
    return (F, Q)
# ------------------------------------------ #
# This function propagates the covariance matrix P
# using Euler integration of the differential Ricatti Equation
def propagate_covariance(P_old, F, G, Q, dt):
    dP = np.dot(F, P_old) + np.dot(P_old, F.T) + np.dot(G, np.dot(Q, G.T))
    P_new = P_old + dP * dt
    return P_new
# ------------------------------------------ #
# These functions are used to avoid MRP singularities.
# They switch to the shadow set if applicable,
def sigmaUnitSwitch(sigma):
    sigma_m = la.norm(sigma)
    sigma = - sigma / (sigma_m * sigma_m)
    return sigma

def check_ShadowTransform(X, P):
    sigma = X[0:3]
    s = la.norm(X[0:3])
    if s > 1.0:
        S = 2.0/(s*s*s*s) * np.outer(sigma, sigma) - 1.0/(s*s) * np.identity(3)
        P_ss = P[0:3, 0:3]
        P_sw = P[0:3, 3:6]
        P[0:3, 0:3] = np.dot(S, np.dot(P_ss, S.T))
        P[0:3, 3:6] = P[3:6, 0:3] = np.dot(S, P_sw)
        X[0:3] = sigmaUnitSwitch(sigma)
    return (X, P)

def apply_ShadowTransform(sigma, P):
    s = la.norm(sigma)
    S = 2.0/(s*s*s*s) * np.outer(sigma, sigma) - 1.0/(s*s) * np.identity(3)
    P_ss = P[0:3, 0:3]
    P_sw = P[0:3, 3:6]
    P_shadow = P.copy()
    P_shadow[0:3, 0:3] = np.dot(S, np.dot(P_ss, S.T))
    P_shadow[0:3, 3:6] = P_shadow[3:6, 0:3] = np.dot(S, P_sw)
    sigma_shadow = rbk.sigmaUnitSwitchCheck(sigma)
    return (sigma_shadow, P_shadow)
# ------------------------------------------ #
# This is a basic Kalman Filter algorithm.
# Note: it doesn't handle the scalar case for measurements y
def Kalman_filter_original(m_old, P_old, u, y, F, G, Q, H, R):
    def compute_KalmanGain(P_minus, H, R):
        Temp = np.dot(H, np.dot(P_minus, (H.T))) + R
        Inv = la.inv(Temp)
        K = np.dot(P_minus, np.dot(H.T, Inv))
        return K

    m_minus = np.dot(F, m_old) + np.dot(G, u)
    P_minus = np.dot(F, np.dot(P_old, F.T)) + Q

    K = compute_KalmanGain(P_minus, H, R)
    innov = y - np.dot(H, m_minus)
    m_plus = m_minus + np.dot(K, innov)

    KH = np.dot(K, H)
    I = np.identity(KH.shape[0])
    M = I - KH
    P_plus = np.dot(M, np.dot(P_minus, M.T)) + np.dot(K, np.dot(R, K.T))
    return (m_plus, P_plus)


# This function computes the Kalman Gain K
def compute_KalmanGain(P_minus, H, R):
    Temp = np.dot(H, np.dot(P_minus, (H.T))) + R
    Inv = la.inv(Temp)
    K = np.dot(P_minus, np.dot(H.T, Inv))
    return K

# This function computes the measurement residual (i.e. innovation factor)
def compute_innovationFactor(y, H, m_minus):
    innov = y - np.dot(H, m_minus)
    if la.norm(innov) > (1.0/3.0):
        innov_s = sigmaUnitSwitch(y) - np.dot(H, m_minus)
        if la.norm(innov_s) < la.norm(innov):
            innov = innov_s
    return innov

# This is the correction step of a conventional Kalman filter algorithm
def computeKalmanCorrection(m_minus, P_minus, y, H, R):
    K = compute_KalmanGain(P_minus, H, R)
    innov = compute_innovationFactor(y, H, m_minus)
    m_plus = m_minus + np.dot(K, innov)
    KH = np.dot(K, H)
    I = np.identity(KH.shape[0])
    M = I - KH
    P_plus = np.dot(M, np.dot(P_minus, M.T)) + np.dot(K, np.dot(R, K.T))
    return (m_plus, P_plus, innov)

# ------------------------------------------ #
def generate_noiseSamples(k_f, R, W_w, W_wb):
    np.random.seed(100)
    m_awgn = np.array([0.0, 0.0, 0.0]) # mean of Additive White Gaussian Noise
    v_samp = np.random.multivariate_normal(m_awgn, R, size=k_f)
    w_nw_samp = np.random.multivariate_normal(m_awgn, W_w, size=k_f)
    w_nwb_samp = np.random.multivariate_normal(m_awgn, W_wb, size=k_f)
    W = la.block_diag(W_w, W_wb)
    return (v_samp, w_nw_samp, w_nwb_samp, W)

# ------------------------------------------ #
def define_H_mat():
    H = np.zeros([3, 6])
    H[:3, :3] = np.identity(3)
    return H
# ------------------- NUMERICAL DATA ------------------- #

def numericData_Karlgaard():
    dt = 0.1
    k_meas = 10
    w_omega = 1E-13 # [rad^2 / s]
    w_bias = 1E-15 # [rad^2 / s^3]
    v = 7.16*1E-5 # [rad^2]
    P0_sigma = 0.0122 * np.identity(3)
    P0_bias = (2.35*1E-9) * np.identity(3)

    W = la.block_diag(w_omega*np.identity(3), w_bias*np.identity(3))
    R = v * np.identity(3)

    P0 = la.block_diag(P0_sigma, P0_bias)
    X0 = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    return(X0, P0, R, W, dt, k_meas)


def numericData_oKeefe():
    sigma0_true = np.array([0.3, 0.1, -0.5])
    omega_true = np.array([-0.2, 0.2, -0.1]) * rbk.D2R
    omegaBias_true = np.array([-1.0, 2.0, 3.0]) * rbk.D2R / 3600.0
    f_sigma = 0.2 # Hz
    f_gyro = 2.0 # Hz
    R = 0.0004 * np.identity(3)

    eta_w = np.sqrt(10) * 1E-7 # rad / s^(1/2)
    eta_wb = np.sqrt(10) * 1E-10 # rad / s^(3/2)

    W_w = eta_w*eta_w * np.identity(3)
    W_wb = eta_wb*eta_wb * np.identity(3)

    P0_sigma = 0.175 * np.identity(3)  # rad^2
    P0_bias = 0.005 * np.identity(3)  # (rad/s)^2
    #P0_bias = 0.1 * np.identity(3)  # (rad/s)^2
    P0 = la.block_diag(P0_sigma, P0_bias)

    tf = 60.0 * 60.0
    t0 = 0.0
    dt = 1.0/f_gyro
    k_f = int(tf/dt)
    k_meas = int(f_gyro / f_sigma)
    span = (tf - t0) / dt + 1
    t_vec = np.linspace(t0, tf, span)

    def printData():
        print 'sigma0_true = ', sigma0_true
        print 'omega_true = ', omega_true
        print 'omegaBias_true = ', omegaBias_true
        print 'P0 = \n', P0
        print 'W_w = \n', W_w
        print 'W_wb = \n', W_wb
        print 'R = \n', R
        print 't_vec = ', t_vec
        print 'dt = ', dt
        print 'k_meas = ', k_meas
        print 'k_f = ', k_f

    printData()
    return (sigma0_true, omega_true, omegaBias_true, P0, W_w, W_wb, R, t_vec, dt, k_meas, k_f)

