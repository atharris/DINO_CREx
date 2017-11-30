import sys, os, inspect
filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))
dinoName = 'DINO_CREx'
splitPath = path.split(dinoName)
batchPath = splitPath[0] + dinoName + '/DINObatch/RngRngRate/'
sys.path.append(batchPath)

import numpy as np
import pdb

from rngRngRtBatch import fncG
test_beacon_state = np.zeros([1,6])
test_state = np.zeros([1,6])
test_beacon_state[0,:] = np.array([2e6,3e7,10,0.01,0.02,0.03])
test_state[0,:] = np.array([1e3,2e4,2,7,2.5,0.3])
extras = dict(obs_beacons='1')
input = (test_state,test_beacon_state,extras)
tol = 1e-12
a = len(extras['obs_beacons'])
G = fncG(input)
test_G = np.zeros((1, 2))

# calculate the difference between the positions and velocites
test_r_diff = test_state[0,0:3] - test_beacon_state[0,0:3]
test_v_diff = test_state[0,3:6] - test_beacon_state[0,3:6]
# calculate the range
test_rng = np.linalg.norm(test_r_diff)
test_G[0,0] = test_rng
# calculate the range rate
test_rng_rate = np.divide(np.inner(test_r_diff, test_v_diff), test_rng)
test_G[0,1] = test_rng_rate
check= abs((G[0,:]-test_G))
if check.any() < tol:
    print('Test passed')
else:
    print('Test Failed')