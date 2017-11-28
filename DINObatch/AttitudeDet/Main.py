# Attitude Determination and Control



# This is the main script to determine the attitude of the cubesat in relation to the object IDs. The script also drives
# the cubesat to the desired attitude.

__author__ = 'Kyle Coker'
__version__ = '$Revision$'[11:-2]
__date__ = '$Date$'[11:5:2017]



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

#
REarth = 6378.14
hLeo = 400
hGeo = 35786
rLeo = REarth + hLeo
rGeo = REarth + hGeo

mu = 398600
thetaD1 = np.sqrt(mu/(rLeo**3))
thetaD2 = np.sqrt(mu/(rGeo**3))

RadSat = np.array([rLeo, 0, 0])
VelSat = np.array([0, rLeo*thetaD1, 0])

RadObj = np.array([rGeo, 0, 0])
VelObj = np.array([0, rGeo*thetaD2, 0])

sigobj = np.array([0, 0.1, 0.2])
sigsat = np.array([0.1, 0, 0.2])

def norm(input):
  norm = np.sqrt(sum(np.square(input)))
  BR = 4 
  return norm

input1 = np.array([[RadSat, VelSat], [sigobj, sigsat]])
input2 = np.array([[RadObj, VelObj]])

rsat = input1[0, 0]
vsat = input1[0, 1]

for jj in range(len(input2)):

    robj = input2[jj, 0]
    vobj = input2[jj, 1]

    # Reference Attitude and Angular speed of beacon

    sobj = input1[1, 0]
    ssat = input1[1, 1]

    bn = MRP2C(sBN)
    hn1 = MRP2C(ssat)
    hn2 = MRP2C(sobj)

    rsatin = hn1.transpose().dot(rsat)
    vsatin = hn1.transpose().dot(vsat)

    robjin = hn2.transpose().dot(robj)
    vobjin = hn2.transpose().dot(vobj)

    rtargetin = robjin - rsatin
    vtargetin = vobjin - vsatin

    r1 = rtargetin/norm(rtargetin)
    r3 = np.cross(rtargetin, vtargetin)/norm(np.cross(rtargetin, vtargetin))
    r2 = np.cross(r3, r1)

    nr = np.matrix([r1.transpose(), r2.transpose(), r3.transpose()])
    rn = nr.transpose()

    srn = C2MRP(rn)

    if norm(srn) > 1:
        srn = -srn/(norm(srn)**2)

    br = bn*rn.transpose()

    sbr = C2MRP(br)

    if norm(sbr) > 1:
        sbr = -sbr/(norm(sbr)**2)

    x = np.concatenate((sBN,wBN))

    K = 1.
    P = 1.

    t = 130
    dt = 1

    time = range(0, t)
    wbn = wBN
    sbrkp = np.zeros([3, t])
    sbnkp = np.zeros([3, t])
    wsckp = np.zeros([3, t])


    for ii in xrange(t):

        shnnew1 = ssat
        shnnew2 = sobj
        
        hnnew1 = MRP2C(shnnew1)
        hnnew2 = MRP2C(shnnew2)

        rsinnew = hnnew1.transpose().dot(rsat)
        roinnew = hnnew2.transpose().dot(robj)
        rtinnew = roinnew - rsinnew

        vsinnew = hnnew1.transpose().dot(vsat)
        voinnew = hnnew2.transpose().dot(vobj)
        vtinnew = voinnew - vsinnew

        r1 = rtinnew / norm(rtinnew)
        r3 = np.cross(rtinnew, vtinnew) / norm(np.cross(rtinnew, vtinnew))
        r2 = np.cross(r3, r1)


        nrnew = np.matrix([r1.transpose(), r2.transpose(), r3.transpose()])
        rnnew = nrnew.transpose()
        srnnew = C2MRP(rnnew)

        sdeltarn = (srnnew - srn)/dt

        BMatrix = BmatMRP(srn)

        wrn = (4/(1+np.dot(srn, srn))**2)*np.dot(BMatrix.transpose(), sdeltarn)
        wbr = wbn - np.dot(np.asarray(br), wrn)

        ## Runge-Kutta goes here
        a = time[ii]

        bb = time[ii] + 1
        m = 1
        x = rk4(a, bb, x, m, sbr, wbr)

        ##

        sbn = x[0:3]
        bn = MRP2C(sbn)
        rn = rnnew
        srn = C2MRP(rn)
        nr = nrnew
        br = bn*nr
        wbn = x[3:6]

        sbr = C2MRP(br)
        s = norm(sbr)

        if s > 1:
            sbr = - sbr/(s**2)

        sbrkp[..., ii] = sbr
        sbnkp[..., ii] = sbn
        wsckp[..., ii] = wbn

    print(sbrkp)













