
# This function is called to solve the MRP kinematic differential equations
# found through equation 3.154 in Schaub and jjunkins


import numpy as np
import pdb


# Inertia Matrix of the Spacecraft - Body
Ibody = np.array([[10, 0, 0], [0, 5, 0],  [0, 0, 7.5]])

# Initializing Arrays
Q = np.zeros([3, 3])
ds = np.zeros(6)


def scfun(yjj, sbr, wbr, K, P):

    # Body-Inertial MRPs
    sig = yjj[0:3]

    # Body-Inertial Angular Rate
    om = yjj[3:6]

    # Body-Target MRPs
    sbr = sbr

    # Body-Target Angular Rate
    wbr = wbr

    # Breaking down the Body-Inertial MRPs for simplicity for the equations below
    s1 = sig[0]
    s2 = sig[1]
    s3 = sig[2]
    ss = np.dot(sig, sig)

    Q[0, 0] = 1 - ss + 2 * s1**2
    Q[0, 1] = 2 * (s1 * s2 - s3)
    Q[0, 2] = 2 * (s1 * s3 + s2)

    Q[1, 0] = 2 * (s1 * s2 + s3)
    Q[1, 1] = 1 - ss + 2 * s2**2
    Q[1, 2] = 2 * (s3 * s2 - s1)

    Q[2, 0] = 2 * (s1 * s3 - s2)
    Q[2, 1] = 2 * (s3 * s2 + s1)
    Q[2, 2] = 1 - ss + 2 * s3**2

    # MRP rate of change
    hh = 0.25*np.dot(Q, om)
    ds[0] = hh[0]
    ds[1] = hh[1]
    ds[2] = hh[2]

    # Control Law
    uu = -K * sbr - P * wbr

    # Angular Acceleration
    ds[3] = (-(Ibody[2, 2] - Ibody[1, 1]) * om[1] * om[2] + uu[0]) / Ibody[0, 0]
    ds[4] = (-(Ibody[0, 0] - Ibody[2, 2]) * om[2] * om[0] + uu[1]) / Ibody[1, 1]
    ds[5] = (-(Ibody[1, 1] - Ibody[0, 0]) * om[0] * om[1] + uu[2]) / Ibody[2, 2]

    # Output these values to the Runge-Kutta for integration

    return ds


