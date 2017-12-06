import numpy as np
from scipy import linalg as la
import matplotlib.pyplot as plt
import KF_commonAlgs as kf


import sys, os
userPath = os.path.expanduser('~')
sys.path.append(userPath + '/Documents/SoftwareDev/Basilisk/PythonModules')
import RigidBodyKinematics as rbk

# -------------------------------------- SUPPORT FUNCTIONS -------------------------------------- #

def define_noiseIntensities():
    w_nw = 1E-1
    w_nwb = 1E-10
    v = 4.0*1E-6
    # AWGN Spectral Densities:
    W_nw = w_nw * np.identity(3)
    W_nwb = w_nwb * np.identity(3)
    V = v * np.identity(3)
    return (W_nw, W_nwb, V)

def define_numericData():
    f_gyro = 10.0 # Hz
    dt = 1.0/ f_gyro
    tf = 800
    k_f = int(tf / dt)
    return (f_gyro, dt, k_f)

def initialConditions():
    def init_refConstant():
        sigma0 = np.array([0.2, 0.4, 0.8])
        omega_true = np.array([0., 0., 0.])  # deg/s 0.8, 0.4, 0.6
        return (sigma0, omega_true)
    (sigma0_true, omega_true) = init_refConstant()
    omega_true *= rbk.D2R
    omega_bias_true = np.array([0.001, 0.0002, 0.0004])  # rad/s
    return (sigma0_true, omega_true, omega_bias_true)

def referenceTrajectory(k_f, dt):
    (sigma0_true, omega_true, omega_bias_true) = initialConditions()
    X0 = np.zeros(6)
    X0[0:3] = sigma0_true.copy()
    X0[3:6] = omega_bias_true.copy()
    sigma_vec = np.zeros([3, k_f])
    omega_b_vec = np.zeros([3, k_f])
    omega_vec = np.zeros([3, k_f])
    k_vec = np.array([])
    params = (omega_true)
    for k in range(0, k_f):
        X0[0:3] = rbk.sigmaUnitSwitchCheck(X0[0:3])
        sigma_vec[:, k] = X0[0:3]
        omega_b_vec[:, k] = X0[3:6]
        omega_vec[:, k] = omega_true
        k_vec = np.append(k_vec, k)
        X = kf.RK4(X0, kf.F_ref, params, dt)
        X0 = X
    return(k_vec, sigma_vec, omega_b_vec, omega_vec)

# -------------------------------------- MAIN MC FUNCTION -------------------------------------- #

def MC_statePropagation(sigma_vec, omega_vec, omega_b_vec, w_nw_samp, w_nwb_samp, W, dt):
    def plotLegend():
        plt.legend(['MC$_1$', 'MC$_2$', 'MC$_3$', 'MC$_4$', 'MC$_5$',
                    'MC$_6$', 'MC$_7$', 'MC$_8$', 'MC$_9$', 'MC$_{10}$',
                    'MC mean',
                    #'MC $2\sigma$','MC $2\sigma$',
                    'Pred mean',
                    #'Pred $2\sigma$','Pred $2\sigma$'
                    ])
    def mcPlotting(color, mark, k_vec, m1_vec, var_str, var_num):
        plt.figure(var_num*100 + 1)
        plt.title('sigma_1')
        plt.plot(k_vec, m1_vec[0, :], color, marker=mark)
        plotLegend()
        plt.figure(var_num*100 + 2)
        plt.title('sigma_2')
        plt.plot(k_vec, m1_vec[1, :], color, marker=mark)
        plotLegend()
        plt.figure(var_num*100 + 3)
        plt.title('sigma_3')
        plt.plot(k_vec, m1_vec[2, :], color, marker=mark)
        plotLegend()

    def MeanCovarPlotting(color, k_vec, m1_vec, s1_vec, var_str, var_num):
        plt.figure(var_num*100 + 1)
        plt.title(var_str + '1')
        plt.plot(k_vec, m1_vec[0, :], color)
        plt.plot(k_vec, m1_vec[0, :] + 2*np.sqrt(s1_vec[0, :]), color, marker='.')
        plt.plot(k_vec, m1_vec[0, :] - 2*np.sqrt(s1_vec[0, :]), color, marker='.')
        plotLegend()
        plt.figure(var_num*100 + 2)
        plt.title(var_str + '2')
        plt.plot(k_vec, m1_vec[1, :], color)
        plt.plot(k_vec, m1_vec[1, :] + 2*np.sqrt(s1_vec[1, :]), color, marker='.')
        plt.plot(k_vec, m1_vec[1, :] - 2*np.sqrt(s1_vec[1, :]), color, marker='.')
        plotLegend()
        plt.figure(var_num*100 + 3)
        plt.title(var_str + '3')
        plt.plot(k_vec, m1_vec[2, :], color)
        plt.plot(k_vec, m1_vec[2, :] + 2*np.sqrt(s1_vec[2, :]), color, marker='.')
        plt.plot(k_vec, m1_vec[2, :] - 2*np.sqrt(s1_vec[2, :]), color, marker='.')
        plotLegend()

    # ------------------------------------------ #

    def statePropagation(x0, P0, m1_sum, s1_sum, m2_sum, s2_sum):
        m1_vec = np.zeros([3, k_f])
        m2_vec = np.zeros([3, k_f])
        s1_vec = np.zeros([3, k_f])
        s2_vec = np.zeros([3, k_f])
        for k in range(0, k_f):
            if la.norm(x0[0:3]) > 1.0:
                (x0[0:3], P0) = kf.apply_ShadowTransform(x0[0:3], P0)
            m1_vec[:, k] = x0[0:3]
            m2_vec[:, k] = x0[3:6]
            std = np.diag(P0)
            s1_vec[:, k] = std[0:3]
            s2_vec[:, k] = std[3:6]
            k_vec[k] = k

            omega_meas = omega_vec[:, k] + omega_b_vec[:, k] + w_nw_samp[k, :]
            params = (omega_meas, w_nw_samp[k, :], w_nwb_samp[k, :])
            x = kf.RK4(x0, kf.F_vec, params, dt)
            P = np.dot(F, np.dot(P0, F.T)) + Q
            (x0, P0) = (x, P)
        m1_sum += m1_vec
        s1_sum += s1_vec
        m2_sum += m2_vec
        s2_sum += s2_vec
        return (m1_vec, m2_vec, s1_vec, s2_vec)

    N_mc = 10
    X_init = np.zeros(6)
    P_init = np.identity(6)
    P_init[:3, :3] *= 0.0125
    P_init[3:,3:] *= 0.0005
    X0_samp = np.random.multivariate_normal(X_init, P_init, size = N_mc)
    k_f = 400 * 2

    k_vec = np.zeros(k_f)

    m1_sum = np.zeros([3, k_f])
    s1_sum = np.zeros([3, k_f])

    m2_sum = np.zeros([3, k_f])
    s2_sum = np.zeros([3, k_f])

    (A, G) = kf.linearizeSystem(sigma_vec[:, 0], omega_vec[:, 0])
    (F, Q) = kf.compute_DTSystem(G, W, A, dt)

    color_vec = np.array(['b', 'r', 'g', 'm', 'k', 'y', 'pink','dodgerblue', 'lightgreen', 'lightsalmon'])

    x_counter = 0

    sigmaNum = 1
    biasNum = 2
    for x0 in X0_samp:
        print 'x0 = ', x0
        P0 = P_init.copy()
        (m1_vec, m2_vec, s1_vec, s2_vec) = statePropagation(x0, P0, m1_sum, s1_sum, m2_sum, s2_sum)
        mcPlotting(color_vec[x_counter], '_', k_vec, m1_vec, 'sigma', sigmaNum)
        mcPlotting(color_vec[x_counter], '_', k_vec, m2_vec, 'bias', biasNum)
        x_counter += 1

    m1_sum *= 1.0/x_counter
    #MeanCovarPlotting('purple', k_vec, m1_sum, s1_sum, 'sigma', sigmaNum)
    #MeanCovarPlotting('purple', k_vec, m2_sum, s2_sum, 'bias', biasNum)
    mcPlotting('purple', '.', k_vec, m1_sum, 'sigma', sigmaNum)
    mcPlotting('purple', '.', k_vec, m2_sum, 'bias', biasNum)

    X0 = X_init.copy()
    print 'x_init = ', X0
    (m1_vec, m2_vec, s1_vec, s2_vec) = statePropagation(X0, P_init, m1_sum, s1_sum, m2_sum, s2_sum)
    #MeanCovarPlotting('violet', k_vec, m1_vec, s1_vec, 'sigma', sigmaNum)
    #MeanCovarPlotting('violet', k_vec, m2_vec, s2_vec, 'bias', biasNum)
    mcPlotting('violet', '.', k_vec, m1_vec, 'sigma', sigmaNum)
    mcPlotting('violet', '.', k_vec, m2_vec, 'bias', biasNum)


(f_gyro, dt, k_f) = define_numericData()
(k_vec, sigma_vec, omega_b_vec, omega_vec) = referenceTrajectory(k_f, dt)
(W_nw, W_nwb, V) = define_noiseIntensities()
np.random.seed(100)
R = 1.0/dt * V
v_samp = np.random.multivariate_normal(np.zeros(3), R, size=k_f)
w_nw_samp = np.random.multivariate_normal(np.zeros(3), W_nw, size=k_f)
w_nwb_samp = np.random.multivariate_normal(np.zeros(3), W_nwb, size=k_f)
W = la.block_diag(W_nw, W_nwb)


MC_statePropagation(sigma_vec, omega_vec, omega_b_vec, w_nw_samp, w_nwb_samp, W, dt)
plt.show()