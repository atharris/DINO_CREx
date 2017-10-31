import numpy as np
from scipy import linalg as la
import matplotlib.pyplot as plt
import KF_commonAlgs as kf
import sys, os
userPath = os.path.expanduser('~')
sys.path.append(userPath + '/Desktop/astroPy/')
import RigidBodyKinematics as rbk

# -------------------------------------- PLOTS DIRECTORY -------------------------------------- #
plotsPath = '/Users/marcolsmargenet/Desktop/FALL_TERM_2016/5044/Final_Project/Figures'
arePlotsSaved = False
if arePlotsSaved:
    plt.rcParams['figure.figsize'] = 2., 2.
    plt.rcParams.update({'font.size': 7})
color_x = 'dodgerblue'
color_y = 'lightgreen'
color_z = 'r'
# -------------------------------------- PLOTTING FUNCTIONS -------------------------------------- #
# This function plots error of mean state and std deviation 2*sigma bounds
def plotMeanError(dt, k_vec, m_vec, s_vec, var_str, units_str, nfig, lim):
    def commonLabels(var_str, dt, var_num):
        plt.legend(['Filter Error $\mu_{'+ var_str + '_'+ str(var_num) +'}$', 'Filter 2$\sigma$ bounds'])
        plt.xlabel('Time, sec')
        plt.ylabel(units_str)
        title = 'Filter error for state $' + var_str + '_'+ str(var_num) +'$ and $\Delta t$ = ' + str(dt)
        if lim > 0:
            plt.ylim([-lim , lim ])
        return title
    varNum = 0
    plt.figure(nfig + varNum)
    t = k_vec * dt
    plt.plot(t, m_vec[0, :], color_x)
    plt.plot(t, 2*s_vec[0, :], 'b', marker='_')
    plt.plot(t, -2*s_vec[0, :], 'b', marker='_')
    title = commonLabels(var_str, dt, varNum + 1)
    savePlot(nfig + varNum, title)

    varNum = 1
    plt.figure(nfig + varNum)
    plt.plot(t, m_vec[1, :], color_y)
    plt.plot(t, 2*s_vec[1, :], 'g', marker='_')
    plt.plot(t, -2*s_vec[1, :], 'g', marker='_')
    title = commonLabels(var_str, dt, varNum + 1)
    savePlot(nfig + varNum, title)

    varNum = 2
    plt.figure(nfig + varNum)
    plt.plot(t, m_vec[2, :], color_z)
    plt.plot(t, 2 * s_vec[2, :], 'lightsalmon', marker='_')
    plt.plot(t, -2 * s_vec[2, :], 'lightsalmon', marker='_')
    title = commonLabels(var_str, dt, varNum + 1)
    savePlot(nfig + varNum, title)

# This function is used to save plots
def savePlot(n, title):
    if arePlotsSaved:
        plt.savefig(plotsPath + "/x" + str(n) + ".pdf", bbox_inches='tight')
    else:
        plt.title(title)
# This function plots MRP states
def plotSigma(k_vec, sigma_vec, nfig, dt, title, lim):
    t = k_vec * dt
    plt.figure(nfig)
    plt.plot(t, sigma_vec[0, :], color_x)
    plt.plot(t, sigma_vec[1, :], color_y),
    plt.plot(t, sigma_vec[2, :], color_z)
    plt.legend(['$\sigma_1$', '$\sigma_2$', '$\sigma_3$'])
    plt.xlabel('Time, sec')
    plt.ylabel('MRP Attitude Set')
    title = title + ' for $\Delta t$ = ' + str(dt)
    if lim > 0:
        plt.ylim([-lim, lim])
    else:
        plt.ylim([-1, 1])
    savePlot(nfig, title)

# This function plots angular velocity states
def plotOmega(k_vec, omega_vec, nfig, dt, title):
    t = k_vec * dt
    plt.figure(nfig)
    plt.plot(t, omega_vec[0, :], color_x)
    plt.plot(t, omega_vec[1, :], color_y)
    plt.plot(t, omega_vec[2, :], color_z)
    plt.legend(['$\omega_1$', '$\omega_2$', '$\omega_3$'])
    plt.xlabel('Time, sec')
    plt.ylabel('Angular velocity, rad/s')
    title = title + ' for $\Delta t$ = ' + str(dt)
    savePlot(nfig, title)

# -------------------------------------- PRINTING FUNCTIONS -------------------------------------- #
def computeTimeInvariantMatrices(W, dt, sigma_true, omega_true):
    (A, G) = kf.linearizeSystem(sigma_true, omega_true)
    (F, Q) = kf.compute_DTSystem(G, W, A, dt)  # G = Tau
    def printMatrices():
        kf.printM6(A, 'A')
        kf.printM6(A, 'F')
        kf.printEigen(A, F)
    return (F, Q)

# -------------------------------------- INTEGRATION FUNCTIONS -------------------------------------- #
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
# This functions compute the time derivatives of the state X = [sigma, omega_bias].T
# returns: dX = [dsigma, domega_bias]

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
#   params = (omega_meas, eta) where P_eta = Q
def F_vec_Q(X, params):
    sigma = X[0:3]
    omega_meas = params[0] # gyro measurement
    eta = params[1]
    omega = omega_meas - X[3:6] - eta[:3]
    dsigma = sigmaDot(sigma, omega)
    domega_bias = eta[3:]
    dX = np.full_like(X, 0.0)
    dX[0:3] = dsigma
    dX[3:6] = domega_bias
    return dX

# State propagated with omega = omega_meas - X[3:6] - eta_w
#   params = (omega_meas, eta_w, eta_wb) where P_eta = [Ww, Wb]
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

# State propagated with pre-computed omega_meas
#   params = (omega_meas, eta_wb)
#   X_m = RK4(X_m0, F_meas, params, dt)
def F_meas(X, params):
    sigma = X[0:3]
    omega_meas = params[0] # gyro measurement
    eta_wb = params[1] # bias_noise
    dsigma = sigmaDot(sigma, omega_meas)
    domega_bias = eta_wb
    dX = np.full_like(X, 0.0)
    dX[0:3] = dsigma
    dX[3:6] = domega_bias
    return dX


# -------------------------------------- KF FUNCTIONS -------------------------------------- #

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
def propagateState(X0, F, G, eta):
    X = np.dot(F, X0) + np.dot(G, eta)
    return X

# -------------------------------------- KF INIT FUNCTIONS -------------------------------------- #
def define_numericData():
    f_gyro = 10.0 # Hz
    dt = 1.0/ f_gyro
    tf = 800
    k_f = int(tf / dt)
    return (f_gyro, dt, k_f)

def define_noiseIntensities():
    w_nw = 1E-1
    w_nwb = 1E-10 # rad / s^(3/2)
    v = 4.0*1E-6 # rad^2
    # AWGN Spectral Densities:
    W_nw = w_nw * np.identity(3)
    W_nwb = w_nwb * np.identity(3)
    V = v * np.identity(3)
    return (W_nw, W_nwb, V)

def initialConditions():
    def init_refVariant():
        sigma0 = np.array([0.0, 0.0, 0.0])
        omega_true = np.array([0.0, 0.8, 0.])  # deg/s 0.8, 0.4, 0.6
        return (sigma0, omega_true)
    def init_refConstant():
        sigma0 = np.array([0.2, 0.4, 0.8])
        #sigma0 = np.array([0.2, 0.2, 0.2])
        #sigma0 = np.array([0.1, 0.2, 0.3])
        #sigma0 = np.array([0.4, 0.4, 0.4])
        omega_true = np.array([0., 0., 0.])  # deg/s 0.8, 0.4, 0.6
        return (sigma0, omega_true)
    def define_nominalBias():
        omega_bias_true = np.array([-1.0, 2.0, -3.0])  # deg/hr
        omega_bias_true *= rbk.D2R * (1.0 / 3600.0)
        return omega_bias_true

    (sigma0_true, omega_true) = init_refConstant()
    omega_true *= rbk.D2R
    omega_bias_true = np.array([0.001, 0.0002, 0.0004])  # rad/s
    return (sigma0_true, omega_true, omega_bias_true)


# -------------------------------------- KF MAIN FUNCTIONS -------------------------------------- #
def referenceTrajectory(k_f, dt):
    def plot_trueStates(k_vec, sigma_vec, omega_b_vec, omega_vec):
        plotSigma(k_vec, sigma_vec, 100, dt, 'True MRP Set', -1)
        plotOmega(k_vec, omega_b_vec, 101, dt, 'True Gyro Bias')
        #plotOmega(k_vec, omega_vec, 102, dt, 'True Angular Velocity')
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
        X = RK4(X0, F_ref, params, dt)
        X0 = X
    #plot_trueStates(k_vec, sigma_vec, omega_b_vec, omega_vec)
    return(k_vec, sigma_vec, omega_b_vec, omega_vec)



# ------------------------------------------ #

def sim_vanillaKalman(X_m0, P0, k_f, sigma_vec, omega_vec, omega_b_vec, w_nw_samp, w_nwb_samp, W, H, R, dt):
    m1_vec = np.zeros([3, k_f])
    m2_vec = np.zeros([3, k_f])
    s1_vec = np.zeros([3, k_f])
    s2_vec = np.zeros([3, k_f])
    k_max = 10
    k_counter = 0
    # Shadow switch check
    kf.check_ShadowTransform(X_m0, P0)
    for i in range(0, k_f):
        # Store states for post-process
        kf.storeKalmanOutputVars(X_m0, P0, m1_vec, m2_vec, s1_vec, s2_vec, k_vec, i)

        # PROPAGATION STEP
        # Create gyro input: omega_meas = omega_true + omega_bias_true + omega_noise
        omega_meas = omega_vec[:, i] + omega_b_vec[:, i] + w_nw_samp[i, :]
        # Propagate non-linear dynamics integrating EOM numerically
        params = (omega_meas, w_nw_samp[i, :], w_nwb_samp[i, :])
        X_minus = RK4(X_m0, F_vec, params, dt)

        # Propagate covariance P
        (A, G) = kf.linearizeSystem(sigma_vec[:, i], omega_vec[:, i])
        (F, Q) = kf.compute_DTSystem(G, W, A, dt) # G = Tau
        P_minus = np.dot(F, np.dot(P0, F.T)) + Q
        kf.check_ShadowTransform(X_minus, P_minus)

        # CORRECTION STEP
        if k_counter == (k_max - 1):
            # Create MRP measurement: sigma_meas = sigma_true + sigma_noise
            sigma_meas = sigma_vec[:, i] + v_samp[i, :]
            # Correct KF estimates
            (X_plus, P_plus, innov) = kf.computeKalmanCorrection(X_minus, P_minus, sigma_meas, H, R)
            kf.check_ShadowTransform(X_plus, P_plus)
            (X_m0, P0) = (X_plus, P_plus)
            k_counter = 0
        else:
            (X_m0, P0) = (X_minus, P_minus)
            k_counter += 1
    print 'KF Final Estimate:'
    print 'X_f = ', X_m0
    print 'std_f = ', np.sqrt(np.diag(P0))
    m1_error = m1_vec - sigma_vec
    m2_error = m2_vec - omega_b_vec
    sigmaLim = 0.1
    omegaLim = 0.005
    plotMeanError(dt, k_vec, m1_error, s1_vec, '\sigma', '', 0, sigmaLim)
    plotMeanError(dt, k_vec, m2_error, s2_vec, '\omega', 'rad/s', 4, omegaLim )
    plotSigma(k_vec, m1_error, 7, dt, 'Mean MRP Error', 0.05)



# ---- Numeric Data ---- #
(f_gyro, dt, k_f) = define_numericData()

# ---- Compute Reference ---- #
(k_vec, sigma_vec, omega_b_vec, omega_vec) = referenceTrajectory(k_f, dt)

# ---- CT Process Noise Parameters ---- #
(W_nw, W_nwb, V) = define_noiseIntensities()
R = 1.0/dt * V
np.random.seed(100)
m_awgn = np.array([0.0, 0.0, 0.0]) # mean of Additive White Gaussian Noise
v_samp = np.random.multivariate_normal(m_awgn, R, size=k_f)
w_nw_samp = np.random.multivariate_normal(m_awgn, W_nw, size=k_f)
w_nwb_samp = np.random.multivariate_normal(m_awgn, W_nwb, size=k_f)
W = la.block_diag(W_nw, W_nwb)

# ---- DT Process Noise Parameters ---- #
(F, Q) = computeTimeInvariantMatrices(W, dt, sigma_vec[:,0], omega_vec[:,0])
eta_samp = np.random.multivariate_normal(np.zeros(6), Q, size=k_f)

# ---- Measurement Matrix ---- #
H = kf.define_H_mat()

# ---- Initial Distribution ---- #
X_m0 = np.zeros(6)
x0 = np.array([0.0, 0.0, 0.0])
P0_sigma = np.identity(3) #* 0.175
P0_bias = np.identity(3) * 0.005  # (rad/s)^2
P0 = la.block_diag(P0_sigma, P0_bias)


if __name__ == "__main__":
    print 'w_nw = ', W_nw[0,0]
    print 'w_nwb = ', W_nwb[0,0]
    print 'R = \n', R
    print 'H = \n', H
    print '\n'
    print 'True Values:'
    print 'sigma0_true = ', sigma_vec[:, 0]
    print 'omegaBias_true = ', omega_b_vec[:, 0]
    print '\n'
    print 'Initial Guess'
    print 'X0 = ', X_m0
    print 'std0 = ', np.sqrt(np.diag(P0))
    print '\n'

    sim_vanillaKalman(X_m0, P0, k_f, sigma_vec, omega_vec, omega_b_vec, w_nw_samp, w_nwb_samp, W, H, R, dt)
    plt.show()