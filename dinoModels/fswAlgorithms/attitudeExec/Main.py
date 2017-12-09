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
import matplotlib.pyplot as plt
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
sigsat = sBN

# Intial States of the Beacon
RadObj1 = np.array([rGeo, 0, 0])
VelObj1 = np.array([0, rGeo*thetaD2, 0])
sigobj1 = np.array([0, 0.1, 0.2])

RadObj2 = np.array([0, rGeo, rGeo])
VelObj2 = np.array([0, rGeo*thetaD2, 0])
sigobj2 = np.array([0.5, 0.5, 0.5])

def norm(input):
  norm = np.sqrt(sum(np.square(input)))
  return norm


# Defining inputs coming in from the Spice data of the beacons and the spacecraft
input1 = np.array([sigobj1, sigobj2])
input2 = np.array([[RadObj1, VelObj1], [RadObj2, VelObj2]])

# Gains for the Control Law
K = 2.
P = 2.
tol = 1E-12
option = 3

rsat = RadSat
vsat = VelSat

# Time for the spacecraft to sweep to desired attitude (seconds)
t = 130
# Initializing Arrays
control_torque = np.zeros([len(input2), t, 3])
sbrskp = np.zeros([len(input2), t, 3])
SigmaBN = np.zeros([len(input2), t, 3])
OmegaBN = np.zeros([len(input2), t, 3])
Observe = np.zeros([len(input2), t])
# Increment Time step (seconds)
dt = 1

tt=[]
sigsat = sBN
wbn = wBN
x = np.concatenate((sigsat, wbn))
for jj in range(len(input2)):

    robj = input2[jj, 0]
    vobj = input2[jj, 1]
    # Reference Attitude and Angular speed of beacon

    sobj = input1[jj]

    # Body-Inertial Frame Rotation Matrix
    bn = MRP2C(sigsat)

    # Hill-Inertial Frame Rotation Matrix
    hn1 = MRP2C(sigsat)
    hn2 = MRP2C(sobj)

    # Inertial Position and Velocity of the Satellite
    rsatin = hn1.transpose().dot(rsat)
    vsatin = hn1.transpose().dot(vsat)

    # Inertial Position and Velocity of the Object
    robjin = hn2.transpose().dot(robj)
    vobjin = hn2.transpose().dot(vobj)
    # Inertial Position and Velocity from Satellite to Object
    rtargetin = robj - rsat
    vtargetin = vobj - vsat

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



    time = range(0, t)




    for ii in xrange(t):

        # Updating the target unit vectors
        r1 = rtargetin / norm(rtargetin)
        r3 = np.cross(rtargetin, vtargetin) / norm(np.cross(rtargetin, vtargetin))
        r2 = np.cross(r3, r1)

        # Update for Inertial - Target Rotation Matrix
        nrnew = np.matrix([r1.transpose(), r2.transpose(), r3.transpose()])
        rnnew = nrnew.transpose()

        # Update for Inertial - Target MRPs
        srnnew = C2MRP(rnnew)

        # MRP for Inertial - Target rate of change
        sdeltarn = (srnnew - srn)/dt

        BMatrix = BmatMRP(srnnew)

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

        sigsat = x[0:3]
        bn = MRP2C(sigsat)
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

        if norm(sbr) < tol:
            observation = 1
        else:
            observation = 0

        if jj == 0:
            indextt = ii
        else:
            indextt = jj*t + jj*ii


        # "Output" Values for the Executive Branch
        Observe[jj, ii] = observation
        OmegaBN[jj, ii, ...] = wbn
        SigmaBN[jj, ii, ...] = sigsat
        control_torque[jj, ii, ...] = (-K*sbr - P*sbr)
        sbrskp[jj, ii, ...] = sbr
        tt.append(indextt)


control_torque.tolist()
OmegaBN.tolist()
SigmaBN.tolist()
sbrskp.tolist()
Observe.tolist()

option = 2

if option == 1:
    plt.plot(tt[0:t], control_torque[0, :, 0], tt[0:t], control_torque[0, :, 1], tt[0:t],
             control_torque[0, :, 2],
             tt[t:2*t], control_torque[1, :, 0], tt[t:2*t], control_torque[1, :, 1], tt[t:2*t],
             control_torque[1, :, 2])
    plt.legend(['$L_{1}$', '$L_{2}$', '$L_3$'])
    plt.ylabel('Control Torque (N $\dot{}$ m)')
    plt.title('For Gains: K =' + str(K) + ' , P = ' + str(P))
elif option == 2:
    plt.plot(tt[0:t], sbrskp[0, :, 0], tt[0:t], sbrskp[0, :, 1], tt[0:t], sbrskp[0, :, 2],
             tt[t:2*t], sbrskp[1, :, 0], tt[t:2*t], sbrskp[1, :, 1], tt[t:2*t], sbrskp[1, :, 2])
    plt.legend(['$\sigma_{1}$', '$\sigma_{2}$', '$\sigma_3$'])
    plt.ylabel('$\sigma_{BR}$')
    plt.title('For Gains: K =' + str(K) + ' , P = ' + str(P))
elif option == 3:
    plt.plot(tt[0:t], Observe[0, :], tt[t:2*t], Observe[1, :])
    plt.legend(['$Beacon_{1}$', '$Beacon_{2}$'])
    plt.ylabel('Observation')
    plt.title('For Gains: K =' + str(K) + ' , P = ' + str(P) + ' and Tol = ' + str(tol))
elif option == 4:
    plt.plot(tt[0:t], SigmaBN[0, :, 0], tt[0:t], SigmaBN[0, :, 1], tt[0:t], SigmaBN[0, :, 2],
             tt[t:2*t], SigmaBN[1, :, 0], tt[t:2*t], SigmaBN[1, :, 1], tt[t:2*t], SigmaBN[1, :, 2])
    plt.legend(['$\sigma_{1}$', '$\sigma_{2}$', '$\sigma_3$'])
    plt.ylabel('$\sigma_{BN}$')
    plt.title('For Gains: K =' + str(K) + ' , P = ' + str(P))
elif option == 5:
    plt.plot(tt[0:t], OmegaBN[0, :, 0], tt[0:t], OmegaBN[0, :, 1], tt[0:t], OmegaBN[0, :, 2],
             tt[t:2*t], OmegaBN[1, :, 0], tt[t:2*t], OmegaBN[1, :, 1], tt[t:2*t], OmegaBN[1, :, 2])
    plt.legend(['$\sigma_{1}$', '$\sigma_{2}$', '$\sigma_3$'])
    plt.ylabel('$\sigma_{BN}$')
    plt.title('For Gains: K =' + str(K) + ' , P = ' + str(P))


plt.grid(True)
plt.xlabel('Time (seconds)')

plt.show()





