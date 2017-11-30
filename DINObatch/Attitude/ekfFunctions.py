import numpy as np
import scipy as sci
import math
import RigidBodyKinematics as rbk

def F_NL(X, options):
    dsigma = 0.25 * rbk.BmatMRP(X[0:3]).dot(options.omega_bn - X[3:6])
    dX = np.full_like(X, 0.0)
    dX[0:3] = dsigma
    return dX

class biasReplacementOptions:
    def __init__(self, omega_bn, F, G, Q):
        self.omega_bn = omega_bn
        self.F = F
        self.G = G
        self.Q = Q

def biasReplacementEOM(t, x, options):
    dx = np.zeros([42,])
    dx[0:3] = 0.25 * rbk.BmatMRP(x[0:3]).dot(options.omega_bn - x[3:6])

    P = np.resize(x[6:], (6,6))
    dF = np.dot(options.F, P) + np.dot(P, options.F.T)
    dG = np.dot(options.G.T, np.dot(options.Q, options.G))
    dP_mat = dF + dG
    n = dP_mat.shape[0] * dP_mat.shape[1]
    dP = np.resize(dP_mat, n)
    dx[6:] = dP

    return dx

def rk4(fun,t0,dt,y0, funOptions):
    #INPUTS:
    # fun: function to be integrated, defined as ynew = fun(t,y)
    # t0: time at y0
    # dt: designated step size (also ref'd as 'h')
    # y0: initial conditions
    k1 = fun(t0,y0,funOptions)
    k2 = fun(  t0+dt/2.0,  y0 + dt/2.0 * k1,funOptions)
    k3 = fun(t0+dt/2.0, y0 + dt/2.0 * k2,funOptions)
    k4 = fun(t0 + dt, y0 + dt*k3,funOptions)
    y1 = y0 + dt/6.0 * (k1 + 2.0*k2 + 2.0*k3 + k4)

    return y1

def skew(vec):
    skewMat = np.array([[ 0, -vec[2], vec[1]], [vec[2], 0, -vec[0]], [-vec[1], vec[0], 0]])
    return skewMat

def linearizeSystem(sigma_ref, omega_ref):
    # Compute F11 submatrix
    swT = np.outer(sigma_ref, omega_ref)
    wsT =  np.outer(omega_ref, sigma_ref)
    w_tilde = skew(omega_ref)
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
    return F,G

