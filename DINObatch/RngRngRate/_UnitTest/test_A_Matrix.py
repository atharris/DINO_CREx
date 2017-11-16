import sys, os, inspect
filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))
bskName = 'Basilisk'
dinoName = 'DINO_CREx'
splitPath = path.split(dinoName)
dinoSpicePath = splitPath[0] + dinoName + '/DINObatch/SPICE/'
bskSpicePath = splitPath[0] + bskName + '/External/EphemerisData/'
bskPath = splitPath[0] + bskName + '/'
sys.path.append(bskPath + 'modules')
sys.path.append(bskPath + 'PythonModules')
bskSpicePath = bskPath + '/External/EphemerisData/'
sys.path.append(dinoSpicePath)
batchPath = splitPath[0] + dinoName + '/DINObatch/RngRngRate/'
sys.path.append(batchPath)
sys.path.append(bskSpicePath)

import numpy as np
import pdb

from posVel import matrixA

test_state = np.array([1e3,2e4,2,7,2.5,0.3])

r_secondary = np.zeros((3,2))
r_secondary [:,0] =  np.array([  7.94442478e+07 , -1.18858180e+08 , -5.15249938e+07])
r_secondary [:,1] =  np.array([  7.94486075e+07 , -1.18858897e+08 , -5.15257383e+07])

r_sun = np.array ([  7.94486075e+07 , -1.18858897e+08 , -5.15257383e+07])

mu_primary = 1.32712428e+11

mu_secondaries = ([398600.4415 , 43050.0])

SRP = 0.3**2/14. * 149597870.**2 * 1358. / 299792458. / 1000
cR = 1
n_secondary = 2

input = (test_state,n_secondary,mu_primary,mu_secondaries,SRP,cR,r_sun,r_secondary)

n_state = len(input[0])
r_spacecraft = np.expand_dims(input[0][0:3], axis=1)
n_secondaries = input[1]
mu_primary = input[2]
mu_secondaries = input[3]
kSRP = input[4]
cR = input[5]
r_sun = np.expand_dims(input[6], axis=1)
r_secondaries_primary = input[7]


test_A = np.zeros((6,6))

test_dFdR_p = -mu_primary * (
    np.identity(3) / np.linalg.norm(r_spacecraft) ** 3 -
    3 * np.dot(r_spacecraft, r_spacecraft.T) / np.linalg.norm(r_spacecraft) ** 5
    )

for ii in xrange(n_secondaries):
    r_secondary = np.expand_dims(r_secondaries_primary[:, ii], axis=1)
    test_dFdR_s = -mu_secondaries[ii] * (
    np.identity(3) / np.linalg.norm(r_spacecraft - r_secondary) ** 3 -
    3 * np.dot(r_spacecraft - r_secondary,(r_spacecraft - r_secondary).T)/
    np.linalg.norm(r_spacecraft - r_secondary) ** 5
)

test_dFdR_SRP = cR * kSRP * (np.identity(3) / np.linalg.norm(r_spacecraft - r_sun) ** 3 -
                        3 * np.dot(r_spacecraft - r_sun, (r_spacecraft - r_sun).T) /
                        np.linalg.norm(r_spacecraft - r_sun) ** 5)

test_dFdR = test_dFdR_p + test_dFdR_s + test_dFdR_SRP

test_A[0, 3] = 1
test_A[1, 4] = 1
test_A[2, 5] = 1
test_A[3:6, 0:3] = test_dFdR

A = matrixA(input)

ErrorMat = A - test_A

testFailure = False
tol = 1e-10
for i in range(np.shape(ErrorMat)[0]):
    for j in range(np.shape(ErrorMat)[1]):
        if np.abs(ErrorMat[i,j])> tol:
            testFailure = True
            print 'TEST FAILED \n'
            print 'Error at row' , i , 'column' , j

if testFailure == False:
    print 'TEST PASSED'

