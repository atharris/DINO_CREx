# This is a unit test to test the Equation of Motions (EOM) found within posVel.py.

import pdb
import numpy as np

def norm( input ):
  norm = np.sqrt(sum(np.square(input)))
  return norm

Position = np.array([70000., -35000., 15000.])
Velocity = np.array([4.5, 2., 0.4])
state = np.concatenate((Position, Velocity))

n_state = np.size(state)
phi = np.identity(n_state)
n_secondaries = 2

mu_primary = 132712428000.0

f_primary = -mu_primary * state[0:3] / norm(state[0:3]) ** 3

# gravitational force from secondary bodies
f_3rd_bodies = 0

# set the size of the r_spacecraftition_secondaries between
# secondary bodies and primary body
r_secondaries_primary = np.zeros((3, n_secondaries))
pdb.set_trace()

# loop through the secondary bodies
for ii in range(n_secondaries):
    # determine distance from secondary to primary body
    r_stateArray = np.zeros(3)
    stateSpice = pyswice.new_doubleArray(6)
    lt = pyswice.new_doubleArray(1)
    pyswice.spkezr_c(bodies[secondary_indices[ii]], et, ref_frame,
                     abcorr, bodies[primary_index], stateSpice, lt)
    for i in range(3):
        r_stateArray[i] = pyswice.doubleArray_getitem(stateSpice, i)
    r_secondaries_primary[:, ii] = r_stateArray

    # calculate the "third body" force
    f_3rd_bodies += -mu_secondaries[ii] * \
                    ((state[0:3] - r_secondaries_primary[:, ii]) / \
                     np.linalg.norm(state[0:3] - r_secondaries_primary[:, ii]) ** 3 + \
                     r_secondaries_primary[:, ii] / np.linalg.norm(r_secondaries_primary[:, ii]) ** 3)

# r_spacecraftition of sun with respect to primary body
r_sunArray = np.zeros(3)
stateSpice = pyswice.new_doubleArray(6)
lt = pyswice.new_doubleArray(1)
pyswice.spkezr_c(bodies[secondary_indices[ii]], et, ref_frame,
                 abcorr, bodies[primary_index], stateSpice, lt)
for i in range(3):
    r_sunArray[i] = pyswice.doubleArray_getitem(stateSpice, i)
r_sun = r_sunArray

# SRP force
f_SRP = cR * kSRP * (state[0:3] - r_sun) / np.linalg.norm(state[0:3] - r_sun) ** 3

# total force (acceleration) vector
f = f_primary + f_3rd_bodies + f_SRP

# args for the A matrix function
args = (state[0:n_state], n_secondaries, mu_primary, mu_secondaries, kSRP, cR, r_sun,
        r_secondaries_primary)

# A matrix calculation
A = matrixA(args)

# calculate the derivative of the STM
dPhi = np.dot(A, phi)
dPhi = np.reshape(dPhi, n_state * n_state)

# acceleration vector to be returned to the integrator
dState = [state[3], state[4], state[5], f[0], f[1], f[2]]
dState += list(dPhi)
pdb.set_trace()
return dState
