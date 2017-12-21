

# Unit test for the Spacecraft Attitude Control Function for the Attitude Executive Branch


import numpy as np
import pdb


# Inertia Matrix of the Spacecraft - Body
Ibody = np.array([[10, 0, 0], [0, 5, 0],  [0, 0, 7.5]])
# Gains for control law
K = 10
P = 10

# Body-Inertial MRPs
sigma_bi = np.array([0.3, -0.4, 0.2])

# Body-Inertial Angular Rate
omega_bi = np.array([0.4, -0.3, 0.2])

# Body-Target MRPs
sigma_bt = np.array([0.2, -0.3, 0.1])

# Body-Target Angular Rate
omega_bt = np.array([0.1, 0.1, 0.])


# Initializing Arrays
Q = np.zeros([3, 3])
ds = np.zeros(6)

# Breaking down the Body-Inertial MRPs for simplicity for the equations below
s1 = sigma_bi[0]
s2 = sigma_bi[1]
s3 = sigma_bi[2]
ss = np.dot(sigma_bi, sigma_bi)

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
hh = 0.25*np.dot(Q, omega_bi)
ds[0] = hh[0]
ds[1] = hh[1]
ds[2] = hh[2]

# Control Law
uu = -K * sigma_bt - P * omega_bt

# Angular Acceleration
ds[3] = (-(Ibody[2, 2] - Ibody[1, 1]) * omega_bi[1] * omega_bi[2] + uu[0]) / Ibody[0, 0]
ds[4] = (-(Ibody[0, 0] - Ibody[2, 2]) * omega_bi[2] * omega_bi[0] + uu[1]) / Ibody[1, 1]
ds[5] = (-(Ibody[1, 1] - Ibody[0, 0]) * omega_bi[0] * omega_bi[1] + uu[2]) / Ibody[2, 2]

pdb.set_trace()

# Output these values to the Runge-Kutta for integration

