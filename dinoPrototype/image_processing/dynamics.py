import numpy as np
import math
import sympy as sym

##################################################
##################################################

def eulerDCM_321(theta1,theta2,theta3):
    ## Computes a Euler 3-2-1 rotation matrix
    #  @param theta1 3rd axis rotation angle [rad]
    #  @param theta2 2nd axis rotation angle [rad]
    #  @param theta3 1st axis rotation angle [rad]
    c1 = math.cos(theta1)
    s1 = math.sin(theta1)
    c2 = math.cos(theta2)
    s2 = math.sin(theta2)
    c3 = math.cos(theta3)
    s3 = math.sin(theta3)
    C = np.array([[c2*c1, c2*s1, -s2],
                  [s3*s2*c1-c3*s1, s3*s2*s1+c3*c1, s3*c2],
                  [c3*s2*c1+s3*s1, c3*s2*s1-s3*c1, c3*c2]])
    return C


##################################################
##################################################

#theta in radians
def eulerDCM_313(theta1,theta2,theta3):
    'Input theta1, theta2, theta3 (rad) to get DCM for Euler 3-1-3 Rotation'
    c1 = math.cos(theta1)
    s1 = math.sin(theta1)
    c2 = math.cos(theta2)
    s2 = math.sin(theta2)
    c3 = math.cos(theta3)
    s3 = math.sin(theta3)
    C = np.array([[c3*c1-s3*c2*s1, c3*s1+s3*c2*c1, s3*s2],
                 [-s3*c1-c3*c2*s1, -s3*s1+c3*c2*c1, c3*s2],
                 [s2*s1, -s2*c1, c2]])
    return C


##################################################
##################################################
# Converts right ascension and declination to unit vector

# Inputs:   radec   tuple of right ascension and declination [deg]
# Outputs   ehat    1x3 array of unit vector

def radec_to_unitvector(radec):

    ehat = np.array([   math.cos(math.radians(radec[1])) * math.cos(math.radians(radec[0])),
                        math.cos(math.radians(radec[1])) * math.sin(math.radians(radec[0])),
                        math.sin(math.radians(radec[1]))
                        ])

    return ehat


##################################################
##################################################
# Converts right ascension and declination to unit vector

# Inputs    ehat    1x3 array of unit vector
# Outputs:  radec   tuple of right ascension and declination [deg]
#                   right ascensions is 0 < RA < 360

def unitvector_to_radec(ehat):

    e1 = ehat[0]
    e2 = ehat[1]
    e3 = ehat[2]

    hplane = math.sqrt(e1**2 + e2**2)

    dec = math.degrees(math.atan2(e3, hplane))
    ra = math.degrees(math.atan2(e2, e1))

    if ra < 0:
        ra = ra + 360

    radec = (ra, dec)

    return radec

def dcm2mrp(C):

    zeta = math.sqrt((C[0,0]+C[1,1]+C[2,2])+1)

    if zeta != 0:
        zeta2 = 1/(zeta*(zeta+2))

        sig1 = zeta2*(C[1,2] - C[2,1])
        sig2 = zeta2*(C[2,0] - C[0,2])
        sig3 = zeta2*(C[0,1] - C[1,0])
        mrp = np.array([sig1, sig2, sig3])
    else:
        phi, ehat = dcm2prv(C)
        beta = prv2ep(phi, ehat)
        mrp = ep2mrp(beta)

    return mrp


def mrp2dcm(mrp):

    if np.linalg.norm(mrp) > 1:
        mrp = mrpshadow(mrp)

    mrp1 = mrp[0]
    mrp2 = mrp[1]
    mrp3 = mrp[2]
    mrpsq = mrp1**2 + mrp2**2 + mrp3**2
    a = 1/(1+mrpsq)**2

    c11 = 4*(mrp1**2-mrp2**2-mrp3**2)+(1-mrpsq)**2
    c12 = 8*mrp1*mrp2+4*mrp3*(1-mrpsq)
    c13 = 8*mrp1*mrp3-4*mrp2*(1-mrpsq)
    c21 = 8*mrp2*mrp1-4*mrp3*(1-mrpsq)
    c22 = 4*(-mrp1**2+mrp2**2-mrp3**2)+(1-mrpsq)**2
    c23 = 8*mrp2*mrp3+4*mrp1*(1-mrpsq)
    c31 = 8*mrp3*mrp1+4*mrp2*(1-mrpsq)
    c32 = 8*mrp3*mrp2-4*mrp1*(1-mrpsq)
    c33 = 4*(-mrp1**2-mrp2**2+mrp3**2)+(1-mrpsq)**2

    C = a*np.array([
        [c11, c12, c13],
        [c21, c22, c23],
        [c31, c32, c33]
    ])

    return C

def crp2dcm(crp):

    q1 = crp[0]
    q2 = crp[1]
    q3 = crp[2]

    sf = 1/(1+q1**2+q2**2+q3**2)
    c11 = sf*(1+q1**2-q2**2-q3**2)
    c12 = sf*2*(q2*q1+q3)
    c13 = sf*2*(q1*q3-q2)
    c21 = sf*2*(q2*q1-q3)
    c22 = sf*(1-q1**2+q2**2-q3**2)
    c23 = sf*2*(q2*q3+q1)
    c31 = sf*2*(q3*q1+q2)
    c32 = sf*2*(q3*q2-q1)
    c33 = sf*(1-q1**2-q2**2+q3**2)

    dcm = np.array([[c11,c12,c13],
                    [c21,c22,c23],
                    [c31,c32,c33]])

    return dcm


def mrp2dcm(mrp):

    if np.linalg.norm(mrp) > 1:
        mrp = mrpshadow(mrp)

    mrp1 = mrp[0]
    mrp2 = mrp[1]
    mrp3 = mrp[2]
    mrpsq = mrp1**2 + mrp2**2 + mrp3**2
    a = 1/(1+mrpsq)**2

    c11 = 4*(mrp1**2-mrp2**2-mrp3**2)+(1-mrpsq)**2
    c12 = 8*mrp1*mrp2+4*mrp3*(1-mrpsq)
    c13 = 8*mrp1*mrp3-4*mrp2*(1-mrpsq)
    c21 = 8*mrp2*mrp1-4*mrp3*(1-mrpsq)
    c22 = 4*(-mrp1**2+mrp2**2-mrp3**2)+(1-mrpsq)**2
    c23 = 8*mrp2*mrp3+4*mrp1*(1-mrpsq)
    c31 = 8*mrp3*mrp1+4*mrp2*(1-mrpsq)
    c32 = 8*mrp3*mrp2-4*mrp1*(1-mrpsq)
    c33 = 4*(-mrp1**2-mrp2**2+mrp3**2)+(1-mrpsq)**2

    C = a*np.array([
        [c11, c12, c13],
        [c21, c22, c23],
        [c31, c32, c33]
    ])

    return C


def dcm2mrp(C):

    zeta = math.sqrt((C[0,0]+C[1,1]+C[2,2])+1)

    if zeta != 0:
        zeta2 = 1/(zeta*(zeta+2))

        sig1 = zeta2*(C[1,2] - C[2,1])
        sig2 = zeta2*(C[2,0] - C[0,2])
        sig3 = zeta2*(C[0,1] - C[1,0])
        mrp = np.array([sig1, sig2, sig3])
    else:
        phi, ehat = dcm2prv(C)
        beta = prv2ep(phi, ehat)
        mrp = ep2mrp(beta)

    return mrp


##################################################
##################################################

def quest(B_ehat, N_ehat, *w):
    ## Attitude estimation using the QUEST method. Variable number of vectors allowed.
    # @param B_ehat List of (1x3) numpy arrays for each unit vector in body coord. frame
    # @param N_ehat List of (1x3) numpy arrays for each unit vector in inertial coord. frame
    # @param w List of weights to be used
    # @return DCM BN direction cosine matrix of attitude solution

    n_obj = len(B_ehat)

    if len(N_ehat) != n_obj:
        print '\nERROR: Uneven number of unit vectors in Quest Method'

    # If weight values are not provided, initialize to all 1
    if not w:

        w = []
        for ind in range(n_obj):
            w.append(1)

    else: w = w[0]


    # Wahba's problem gain function g(Beta) = Beta^T * [K] * Beta
    # Solve for [K] matrix in Wahba's problem gain function

    B = np.zeros((3, 3))

    for indEhat in range(n_obj):

        Btemp = w[indEhat] * np.outer(B_ehat[indEhat], N_ehat[indEhat])
        B = B+Btemp

    S = B+B.T

    sig = np.trace(B)

    Ztrans = np.array([B[1,2] - B[2,1], B[2,0] - B[0,2], B[0,1] - B[1,0]])

    Kup = np.append(sig, Ztrans)
    Klowright = S - sig * np.identity(3)

    Z = Ztrans.reshape(-1, 1)
    Klow = np.concatenate((Z, Klowright), axis=1)

    K = np.vstack((Kup, Klow))

    # characteristic equation for eigenvalues of K
    s = sym.Symbol('s')
    fmat = sym.Matrix(K - s * np.identity(4))
    f = fmat.det()

    # first eigenvalue set to sum of weights
    lam0 = sum(w)

    # solve characteristic equation for other eigenvalues with Newton-Rhapson method
    fdiff = sym.diff(f, s)
    lam1 = lam0 - f.subs(s, lam0) / fdiff.subs(s, lam0)
    lam2 = lam1 - f.subs(s, lam1) / fdiff.subs(s, lam1)
    lam3 = lam2 - f.subs(s, lam2) / fdiff.subs(s, lam2)
    lam4 = lam3 - f.subs(s, lam3) / fdiff.subs(s, lam3)

    # solve for CRP
    mat1 = (lam4 + sig) * np.identity(3)
    mat1 = np.array(mat1, dtype='float')
    mat1 = mat1 - S
    mat1 = np.array(mat1, dtype='float')
    matinv = np.linalg.inv(mat1)
    crp = np.matmul(matinv, Ztrans)

    # convert CROP to a DCM
    DCM = crp2dcm(crp)

    return DCM