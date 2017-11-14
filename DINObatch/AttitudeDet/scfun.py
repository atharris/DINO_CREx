
# This function is called to solve the MRP kinematic differential equations
# found through equation 3.154 in Schaub and Junkins


import numpy as np
import pdb
K = 10 
P = 10 

Ibody = np.array([[10, 0, 0], [0, 5, 0],  [0, 0, 7.5]])

Q = np.zeros([3, 3])
ds = np.zeros([1, 6])


def scfun(yj, sbr, wbr):
    pdb.set_trace()
    sig = yj[0:3, 0]
    om = yj[3:6, 0]
    sbr = sbr
    wbr = wbr

    s1 = sig(1)
    s2 = sig(2)
    s3 = sig(3)
    ss = sig.transpose()*sig

    Q[1, 1] = 1 - ss + 2 * s1 ^ 2
    Q[1, 2] = 2 * (s1 * s2 - s3)
    Q[1, 3] = 2 * (s1 * s3 + s2)

    Q[2, 1] = 2 * (s1 * s2 + s3)
    Q[2, 2] = 1 - ss + 2 * s2 ^ 2
    Q[2, 3] = 2 * (s3 * s2 - s1)

    Q[3, 1] = 2 * (s1 * s3 - s2)
    Q[3, 2] = 2 * (s3 * s2 + s1)
    Q[3, 3] = 1 - ss + 2 * s3 ^ 2

    ds[1:3] = 0.25*Q*om

    u = -K * sbr - P * wbr

    ds[4] = (-(Ibody[9] - Ibody[5]) * om[2] * om[3] + u[1]) / Ibody[1]
    ds[5] = (-(Ibody[1] - Ibody[9]) * om[3] * om[1] + u[2]) / Ibody[5]
    ds[6] = (-(Ibody[5] - Ibody[1]) * om[1] * om[2] + u[3]) / Ibody[9]

    return ds


