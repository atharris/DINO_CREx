import numpy as np
from scipy import linalg as la
import matplotlib.pyplot as plt
import KF_commonAlgs as kf
import FP_vanillaKalman as vk
import sys, os
userPath = os.path.expanduser('~')
sys.path.append(userPath + '/Desktop/astroPy/')
import RigidBodyKinematics as rbk



def X_truthModel():
    X_trueModel = np.zeros(6)
    return X_trueModel

def create_gyroMeas(omega_true, bias_true, eta_w):
    omega_meas = omega_true + bias_true + eta_w
    return omega_meas

def create_stMeas(sigma_true):
    v_samp = np.random.multivariate_normal(np.zeros(3), vk.R, size=1)
    sigma_meas = sigma_true + v_samp[0,:]
    return sigma_meas

def computeLTI():
    F = vk.F
    Q = vk.Q
    H = vk.H
    R = vk.R
    return (F, Q, H, R)

def X_KF_estimate(X_old, P_old, k_counter, k_max, sigma_true, omega_true, bias_true):
        eta_w = np.random.multivariate_normal(np.zeros(3), vk.W_nw, size=1)[0,:]
        eta_b = np.random.multivariate_normal(np.zeros(3), vk.W_nwb, size=1)[0,:]

        # PROPAGATION STEP
        omega_meas = create_gyroMeas(omega_true, bias_true, eta_w)
        params = (omega_meas, eta_w, eta_b)
        X_minus = kf.RK4(X_old, kf.F_vec, params, dt)
        P_minus = np.dot(vk.F, np.dot(P_old, vk.F.T)) + vk.Q
        kf.check_ShadowTransform(X_minus, P_minus)
        # CORRECTION STEP
        if k_counter == (k_max - 1):
            sigma_meas = create_stMeas(sigma_true)
            (X_plus, P_plus, innov) = kf.computeKalmanCorrection(X_minus, P_minus, sigma_meas, vk.H, vk.R)
            kf.check_ShadowTransform(X_plus, P_plus)
            k_counter = 0
            return(X_plus, P_plus, k_counter)
        else:
            k_counter += 1
            return(X_minus, P_minus, k_counter)

(f_gyro, dt, k_f) = vk.define_numericData()
(sigma0_true, omega_true, omega_bias_true) = vk.initialConditions()

# ---- Initial Distribution ---- #
X_m0 = np.zeros(6)
x0 = np.array([0.0, 0.0, 0.0])
P0_sigma = np.identity(3) #* 0.175
P0_bias = np.identity(3) * 0.005  # (rad/s)^2
P0 = la.block_diag(P0_sigma, P0_bias)
kf.check_ShadowTransform(X_m0, P0)


# ---- Begin ---- #

def vanilla_KF(X_m0, P0):
    m1_vec = np.zeros([3, k_f])
    m2_vec = np.zeros([3, k_f])
    s1_vec = np.zeros([3, k_f])
    s2_vec = np.zeros([3, k_f])
    k_vec = np.zeros(k_f)
    k_max = 10
    k_counter = 0
    X_old = X_m0.copy()
    P_old = P0.copy()
    for i in range(0, k_f):
        kf.storeKalmanOutputVars(X_old, P_old, m1_vec, m2_vec, s1_vec, s2_vec, k_vec, i)
        (X, P, k_counter) = X_KF_estimate(X_old, P_old, k_counter, k_max, sigma0_true, omega_true, omega_bias_true)
        (X_old, P_old) = (X, P)
    print 'X_f = ', X_old
    print 'std_f = ', np.sqrt(np.diag(P_old))


    (k_vec, sigma_vec, bias_vec, omega_vec) = vk.referenceTrajectory(k_f, dt)
    m1_error = m1_vec - sigma_vec
    m2_error = m2_vec - bias_vec
    sigmaLim = 0.1
    omegaLim = 0.005
    vk.plotMeanError(dt, k_vec, m1_error, s1_vec, '\sigma', '', 0, sigmaLim)
    vk.plotMeanError(dt, k_vec, m2_error, s2_vec, '\omega', 'rad/s', 4, omegaLim)
    vk.plotSigma(k_vec, m1_error, 7, dt, 'Mean MRP Error', 0.05)


def chi2inv(N_MC):
    # For n = 6 and alpha = 0.005
    if N_MC == 1:
        r1 = 0.5266
        r2 =  20.2494
    elif N_MC == 5:
        r1 = 2.5529
        r2 = 11.2665
    else:
        raise ValueError('Chi 2 inv bounds not defined for N_MC = ' + str(N_MC))
    return (r1, r2)

def ChiSquareTest_old(X_m0, P0):
    k_max = 100
    k_counter = 0
    X_old = X_m0.copy()
    P_old = P0.copy()

    n_MC = 1
    X_true = np.append(sigma0_true, omega_bias_true)
    P_true = la.block_diag(np.identity(3), np.zeros([3,3]))
    x_true_samp = np.random.multivariate_normal(X_true, P_true, size=n_MC)
    (r1, r2) = chi2inv(n_MC)

    plt.figure()

    E_sum_vec = np.zeros(k_f)
    k_vec = np.zeros(k_f)
    for i in range(0, n_MC):
        print 'i_MC = ', i
        X_true_MC = x_true_samp[i]
        X_true_MC = X_true
        print 'X_true_MC = ', X_true_MC
        for k in range(0, k_f):
            (X, P, k_counter) = X_KF_estimate(X_old, P_old, k_counter, k_max, X_true_MC[:3], omega_true, X_true_MC[3:])
            (X_old, P_old) = (X, P)
            e = X_true_MC - X_old
            E = np.dot(e, np.dot(la.inv(P), e))
            #print 'E = ', E
            E_sum_vec[k] += E
            k_vec[k] = k
        print 'E_sum_vec (i_MC) = ', E_sum_vec
    E_sum_vec *= (1.0/n_MC)
    for k in range(0, k_f):
        if (k % k_max == 0):
            print 'E = ', E_sum_vec[k]
            plt.scatter(k, E_sum_vec[k])
    plt.plot(k_vec, E_sum_vec, 'g', marker='o')
    plt.axhline(r1, color='lime', marker = '_')
    plt.axhline(r2, color='lime',marker = '_')
    plt.show()
    return

vanilla_KF(X_m0, P0)
ChiSquareTest_old(X_m0, P0)
plt.show()