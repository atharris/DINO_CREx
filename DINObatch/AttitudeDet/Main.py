# Attitude Determination and Control



# This is the main script to determine the attitude of the cubesat in relation to the object IDs. The script also drives
# the cubesat to the desired attitude.

__author__ = 'Kyle Coker'
__version__ = '$Revision$'[11:-4]
__date__ = '$Date$'[11:28:2017]

#######################################################################################################
#                                        IMPORT

import numpy as np
import pdb

import sys, os, inspect

filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))

dinoName = 'DINO_CREx'
splitPath = path.split(dinoName)
rigidPath = splitPath[0] + '/ Basilisk/ PythonModules/'
sys.path.append(rigidPath)

from RigidBodyKinematics import C2EP, C2MRP, MRP2C, BmatMRP
from rk4 import rk4

# Reference Attitude and Angular Speed of satellite

sBN = np.array([0.3, -0.4, 0.5])
wBN = np.array([1, 1.75, -2.2])*np.pi/180

# Scenario: Orbit Parameters
REarth = 6378.14
hLeo = 400
hGeo = 35786
rLeo = REarth + hLeo
rGeo = REarth + hGeo

# Earth's Gravitational Parameter
mu = 398600
# True Anomaly of the spacecraft
thetaD1 = np.sqrt(mu/(rLeo**3))
# True Anomaly of the Beacon
thetaD2 = np.sqrt(mu/(rGeo**3))

# Initial States of the Spacecraft
RadSat = np.array([rLeo, 0, 0])
VelSat = np.array([0, rLeo*thetaD1, 0])
sigsat = np.array([0.1, 0, 0.2])

# Intial States of the Beacon
RadObj = np.array([rGeo, 0, 0])
VelObj = np.array([0, rGeo*thetaD2, 0])
sigobj = np.array([0, 0.1, 0.2])


def norm(input):
  norm = np.sqrt(sum(np.square(input)))
  return norm


# Defining inputs coming in from the Spice data of the beacons and the spacecraft
input1 = np.array([[RadSat, VelSat], [sigobj, sigsat]])
input2 = np.array([[RadObj, VelObj]])

# Gains for the Control Law
K = 10.
P = 10.

rsat = input1[0, 0]
vsat = input1[0, 1]

for jj in range(len(input2)):

    robj = input2[jj, 0]
    vobj = input2[jj, 1]

    # Reference Attitude and Angular speed of beacon

    sobj = input1[1, 0]
    ssat = input1[1, 1]

    # Body-Inertial Frame Rotation Matrix
    bn = MRP2C(sBN)

    # Hill-Inertial Frame Rotation Matrix
    hn1 = MRP2C(ssat)
    hn2 = MRP2C(sobj)

    # Inertial Position and Velocity of the Satellite
    rsatin = hn1.transpose().dot(rsat)
    vsatin = hn1.transpose().dot(vsat)

    # Inertial Position and Velocity of the Object
    robjin = hn2.transpose().dot(robj)
    vobjin = hn2.transpose().dot(vobj)

    # Inertial Position and Velocity from Satellite to Object
    rtargetin = robjin - rsatin
    vtargetin = vobjin - vsatin

    # Target unit vectors
    r1 = rtargetin/norm(rtargetin)
    r3 = np.cross(rtargetin, vtargetin)/norm(np.cross(rtargetin, vtargetin))
    r2 = np.cross(r3, r1)

    # Inertial - Target Rotation Matrix
    nr = np.matrix([r1.transpose(), r2.transpose(), r3.transpose()])
    rn = nr.transpose()

    # MRP of Inertial-Target
    srn = C2MRP(rn)

    # Check magnitude and switch to Shadow MRPs if necessary
    if norm(srn) > 1:
        srn = -srn/(norm(srn)**2)

    # Body - Target rotation Matrix
    br = bn*rn.transpose()

    # MRPs for Body-Target frame
    sbr = C2MRP(br)

    # Checking magnitude and switching to shadow set MRPs if necessary
    if norm(sbr) > 1:
        sbr = -sbr/(norm(sbr)**2)

    x = np.concatenate((sBN, wBN))

    # Time for the spacecraft to sweep to desired attitude (seconds)
    t = 130

    # Increment Time step (seconds)
    dt = 1

    time = range(0, t)
    wbn = wBN

    # Initializing Arrays
    control_torque = np.zeros([t, 3])

    for ii in xrange(t):

        # New Hill-Inertial MRPs
        shnnew1 = ssat
        shnnew2 = sobj

        # New Hill-Inertial Rotation Matrices
        hnnew1 = MRP2C(shnnew1)
        hnnew2 = MRP2C(shnnew2)

        # Updated position and velocity for Spacecraft
        rsinnew = hnnew1.transpose().dot(rsat)
        vsinnew = hnnew1.transpose().dot(vsat)

        # Updated position and velocity for Object
        roinnew = hnnew2.transpose().dot(robj)
        voinnew = hnnew2.transpose().dot(vobj)

        # Updated position and velocity for Spacecraft - Target
        rtinnew = roinnew - rsinnew
        vtinnew = voinnew - vsinnew

        # Updating the target unit vectors
        r1 = rtinnew / norm(rtinnew)
        r3 = np.cross(rtinnew, vtinnew) / norm(np.cross(rtinnew, vtinnew))
        r2 = np.cross(r3, r1)

        # Update for Inertial - Target Rotation Matrix
        nrnew = np.matrix([r1.transpose(), r2.transpose(), r3.transpose()])
        rnnew = nrnew.transpose()

        # Update for Inertial - Target MRPs
        srnnew = C2MRP(rnnew)

        # MRP for Inertial - Target rate of change
        sdeltarn = (srnnew - srn)/dt

        BMatrix = BmatMRP(srn)

        # Angular Rate of Inertial - Target
        wrn = (4/(1+np.dot(srn, srn))**2)*np.dot(BMatrix.transpose(), sdeltarn)

        # Angular Rate of Body - Target (i.e. The desired angular rate)
        wbr = wbn - np.dot(np.asarray(br), wrn)

        ## Runge-Kutta goes here

        # Time inputs for Runge-Kutta
        a = time[ii]
        bb = time[ii] + 1

        # Number of Steps for the Runge-Kutta to run through
        m = 1

        # Outputs from Runge-Kutta
        x = rk4(a, bb, x, m, sbr, wbr, K, P)

        # Update all previous values for next iteration

        sbn = x[0:3]
        bn = MRP2C(sbn)
        rn = rnnew
        srn = C2MRP(rn)
        nr = nrnew
        br = bn*nr
        wbn = x[3:6]

        sbr = C2MRP(br)
        s = norm(sbr)

        # Check the magnitude of the Body-Target MRPs and switch to shadow set if necessary
        if s > 1:
            sbr = - sbr/(s**2)

        # "Output" Values for the Executive Branch
        control_torque[ii, ...] = (-K*sbr - P*sbr)

    print(control_torque)













