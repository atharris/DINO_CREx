#! /usr/bin/env python3
# H+
#	Title   : camera_demo.py
#	Author  : Matt Muszynski
#	Date    : 03/23/17
#	Synopsis: 
#
#
#	$Date$
#	$Source$
#  @(#) $Revision$
#	$Locker$
#
#	Revisions:
#
# H-
# U+
#	Usage   :
#	Example	:
#	Output  :
# U-
# D+
#
# D-
###############################################################################

import sys
sys.path.insert(0, 'dependencies')
from datetime import datetime
import numpy as np
from numpy import deg2rad
import sqlite3
import pdb
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LogNorm
import matplotlib.colors as colors
import bodies as bod
import camera
from propagator import coe_2body
import numpy.linalg as la
from matplotlib.animation import FuncAnimation
from constants import au
import adcs
from adcs import Euler321_2DCM
from numpy.random import uniform

print('################## Initializing ##################')

makeAll = 1

###############################################################################
#
#       Create a spacecraft object to attach the camera to.
#       A lot of the details of the spacecraft don't matter because
#       they aren't used in the camera model. The really important one is
#       the state vector.
#
###############################################################################

sc = bod.sc(
        "SC", #name of this body
        "Earth", #central body
        2451545.0, #Epoch of coe presented here, in JD
        122164 , #a in km
        0, # e
        0, #inclination in deg
        0, #omega in deg
        0, #OMEGA in deg
        0, #Mean Anomaly at epoch in deg
        np.array([]), #state vector in HCI frame. [x,y,z,v_x,v_y,v_z]
        np.nan, #true anomaly measured at same time as state vector
        np.nan
        )

###############################################################################
#
#       Pull in canned QE and Transmission curves from DINO C-REx files
#
###############################################################################

#load tranmission curve for Canon 20D
_20D = np.load('tc/20D.npz')
tc = {}
tc['lambda'] = _20D['x']
tc['throughput'] = _20D['y']

#load QE curve for Hubble Space Telecope Advanced Camera for Surveys SITe CCD
ACS = np.load('qe/ACS.npz')
qe = {}
qe['lambda'] = ACS['x']
qe['throughput'] = ACS['y']




###############################################################################
#
#       Initialize camera
#
###############################################################################


msg = { 'bodies': [
	bod.earth,
	bod.luna,
	sc
	], 
	'addStars': 1,'rmOcc': 1, 'addBod': 1, 'psf': 1, 
	'raster': 1, 'photon': 0, 'dark': 0, 'read': 0, 'dt': 0.01}

cam = camera.camera(
	0.019968, 				#detector_height
	0.019968, 				#detector_width
	0.05, 			#focal_length
	512, 			#resolution_height
	512,			#resolution_width
	np.identity(3), #body2cameraDCM
	1000,		    #maximum magnitude
	-1000,			#minimum magnitude (for debugging)
	qe,
	tc,
	1,
	0.01**2, #effective area in m^2
	100, #dark current in electrons per second
	100, #std for read noise in electrons
	1000, #bin size
	2**32, #max bin depth
	1,
	sc,
	msg
	)

starCam = camera.camera(
	0.019968, 				#detector_height
	0.019968, 				#detector_width
	0.05, 			#focal_length
	512, 			#resolution_height
	512,			#resolution_width
	np.identity(3), #body2cameraDCM
	1000,		    #maximum magnitude
	-1000,			#minimum magnitude (for debugging)
	qe,
	tc,
	1,
	0.01**2, #effective area in m^2
	100, #dark current in electrons per second
	100, #std for read noise in electrons
	1000, #bin size
	2**32, #max bin depth
	1,
	sc,
	msg
	)

###############################################################################
#
#		Image 1 (Earth as Resolved Body, Moon as a Pt Source)
#
###############################################################################

if 0 or makeAll:

	sc.state = np.array([au/1000,-1e6, 0, 0, 0, 0])

	bod.earth.state = np.array([au/1000, 0, 0, 0, 0, 0])

	bod.luna.state = np.array([au/1000, 0, 0, 0, 0, 0]) + \
		np.array([34862, -601522, 0, 0, 0 ,0])

	sc.attitudeDCM = np.array([
		[ 0, 1, 0],
		[-1, 0, 0],
		[ 0, 0, 1]
		])
	msg['takeImage'] = 1
	cam.updateState()
	msg['takeImage'] = 0
	cam.updateState()

	plt.figure()
	plt.imshow(cam.images[len(cam.images)-1].detectorArray.reshape(512,512))


###############################################################################
#
#		Image 2 (Earth as Resolved Body, Moon as a Pt Source)
#
###############################################################################

if 0 or makeAll:

	sc.state = np.array([au/1000,-4e6, 0, 0, 0, 0])

	bod.earth.state = np.array([au/1000, 0, 0, 0, 0, 0])

	bod.luna.state = np.array([au/1000, 0, 0, 0, 0, 0]) + \
		np.array([34862, -3.601522E6, 0, 0, 0 ,0])

	sc.attitudeDCM = np.array([
		[ 0, 1, 0],
		[-1, 0, 0],
		[ 0, 0, 1]
		])
	msg['takeImage'] = 1
	cam.updateState()
	msg['takeImage'] = 0
	cam.updateState()

	plt.figure()
	plt.imshow(cam.images[len(cam.images)-1].detectorArray.reshape(512,512))


# ###############################################################################
# #
# #		Image 3 (Beacons as large resolved bodies)
# #
# ###############################################################################

if 0 or makeAll:
	sc.state = np.array([au/1000,-1e6, 0, 0, 0, 0])

	bod.earth.state = np.array([au/1000, 0, 0, 0, 0, 0])

	bod.luna.state = np.array([au/1000, 0, 0, 0, 0, 0]) + \
		np.array([34862, -9E5, 0, 0, 0 ,0])

	sc.attitudeDCM = np.array([
		[ 0, 1, 0],
		[-1, 0, 0],
		[ 0, 0, 1]
		])
	msg['takeImage'] = 1
	cam.updateState()
	msg['takeImage'] = 0
	cam.updateState()

	plt.figure()
	plt.imshow(cam.images[len(cam.images)-1].detectorArray.reshape(512,512))

###############################################################################
#
#		Image 4 (Blended Source)
#
###############################################################################

if 0 or makeAll:
	sc.state = np.array([au/1000,-5e5, 0, 0, 0, 0])

	bod.earth.state = np.array([au/1000, 0, 0, 0, 0, 0])

	bod.luna.state = np.array([au/1000, 0, 0, 0, 0, 0]) + \
		np.array([-2500, 0, 4000, 0, 0 ,0])

	tmpMoonR = bod.luna.r_eq
	bod.luna.r_eq = bod.earth.r_eq
	sc.attitudeDCM = np.array([
		[ 0, 1, 0],
		[-1, 0, 0],
		[ 0, 0, 1]
		])
	msg['takeImage'] = 1
	cam.updateState()
	msg['takeImage'] = 0
	cam.updateState()

	plt.figure()
	plt.imshow(cam.images[len(cam.images)-1].detectorArray.reshape(512,512))

###############################################################################
#
#		Image 5 (Corners)
#
###############################################################################

if 0 or makeAll:
	sc.state = np.array([au/1000,-5e5, 0, 0, 0, 0])

	bod.earth.state = np.array([au/1000, 0, 0, 0, 0, 0]) + \
		np.array([1.025e5, 0, -1e5, 0, 0 ,0])

	bod.luna.state = np.array([au/1000, 0, 0, 0, 0, 0]) + \
		np.array([-9.5e4, 0, 9.5e4, 0, 0 ,0])

	sc.attitudeDCM = np.array([
		[ 0, 1, 0],
		[-1, 0, 0],
		[ 0, 0, 1]
		])
	msg['takeImage'] = 1
	cam.updateState()
	msg['takeImage'] = 0
	cam.updateState()

	plt.figure()
	plt.imshow(cam.images[len(cam.images)-1].detectorArray.reshape(512,512))

###############################################################################
#
#		Image 6 (Not Blended, but ROIs overlap)
#
###############################################################################

if 0 or makeAll:
	sc.state = np.array([au/1000,-5e5, 0, 0, 0, 0])

	bod.earth.state = np.array([au/1000, 0, 0, 0, 0, 0])

	bod.luna.state = np.array([au/1000, 0, 0, 0, 0, 0]) + \
		np.array([-7000, 0, 4000, 0, 0 ,0])

	tmpMoonR = bod.luna.r_eq
	bod.luna.r_eq = bod.earth.r_eq
	sc.attitudeDCM = np.array([
		[ 0, 1, 0],
		[-1, 0, 0],
		[ 0, 0, 1]
		])
	msg['takeImage'] = 1
	cam.updateState()
	msg['takeImage'] = 0
	cam.updateState()

	plt.figure()
	plt.imshow(cam.images[len(cam.images)-1].detectorArray.reshape(512,512))

###############################################################################
#
#		Image 7 (Hot Pixels)
#
###############################################################################

if 0 or makeAll:
	sc.attitudeDCM = np.array([
		[ 0, 1, 0],
		[-1, 0, 0],
		[ 0, 0, 1]
		])

	msg['addBod'] = 0
	msg['rmOcc'] = 0

	msg['takeImage'] = 1
	starCam.updateState()
	msg['takeImage'] = 0
	starCam.updateState()

	location = np.array([  67238,  204676,    8502,    3127,   66954,   37175,
	        198832,  143308,  171596,  148843,   36628,  139528,
	        211838,  168193,   94000,  153682,   83715,  151847,
	        110000,   86087,  137886,  245761,  126063,   34410,
	         71959,   65053,  198715,   26613,   81186,  238110,
	        116165,   13751,  222347,   39720,   70698,   12145,
	        229304,  175423,  114801,  230443,  255514,   70676,
	        177123,   85416,   28670,   79794,   76463,   56513,
	        251937,  176083,  223210,   69272,  166395,  135999,
	        214500,   82304,  135686,   95969,  222475,  197569,
	         44762,   41205,   48567,   20308,   50656,  108693,
	        226143,  180566,   71751,  184761,  113617,  138538,
	         89182,  164795,   70876,  166511,   80919,  106320,
	        130322,   33213,  193733,   64612,   54932,  180376,
	        122599,   75537,   83099,  103414,  185201,   89593,
	        237452,  121161,   85407,  224852,  176546,    3958,
	        229694,   58432,   80575,  139825])

	val = np.array([  1.57159880e+09,   1.98555011e+09,   2.85284607e+09,
         1.80362116e+09,   8.79758079e+08,   5.44505254e+08,
         6.35932276e+08,   1.24185567e+09,   3.12418006e+09,
         5.66645652e+08,   3.58988980e+09,   4.46815068e+08,
         3.69032977e+09,   4.12124175e+09,   4.27757839e+09,
         2.39130333e+09,   3.90653374e+09,   1.94928756e+08,
         3.92622941e+09,   1.05477115e+09,   2.77483925e+09,
         3.20591270e+09,   2.52172439e+09,   2.19811753e+08,
         3.06862038e+09,   1.76635087e+08,   3.22916703e+09,
         2.02511379e+09,   3.53739037e+09,   9.74116192e+08,
         9.01350595e+07,   5.18526954e+08,   3.14735436e+09,
         1.15382492e+09,   8.03171441e+08,   3.20822869e+09,
         3.56730659e+09,   1.78021667e+09,   1.96956049e+09,
         3.07352237e+09,   1.38404406e+09,   1.32903194e+09,
         3.20982962e+09,   7.11457694e+08,   2.23783572e+09,
         1.52547654e+09,   6.43052018e+08,   1.81298153e+09,
         3.84288041e+08,   3.88439744e+09,   3.02720398e+09,
         5.00463803e+06,   3.62870716e+09,   3.62430177e+09,
         1.53281628e+09,   3.07906954e+09,   3.47965714e+09,
         2.25428634e+09,   3.21823533e+09,   2.42523942e+09,
         1.56843470e+09,   3.35810259e+09,   9.43395350e+08,
         7.16231196e+08,   2.03147121e+09,   1.47226573e+08,
         3.94249594e+09,   3.08576703e+09,   1.80063064e+09,
         2.53298980e+09,   3.22935951e+09,   1.18027328e+08,
         1.93788056e+09,   3.23944939e+09,   8.20264454e+08,
         3.38227856e+09,   2.53068987e+09,   2.77471913e+09,
         1.12539277e+09,   3.20333961e+09,   4.09704746e+09,
         2.61630188e+09,   3.71144967e+08,   1.69423439e+09,
         2.39992578e+09,   2.92584988e+09,   3.64039911e+09,
         1.18798527e+09,   3.31310671e+09,   3.83692086e+09,
         3.90837016e+09,   2.21377588e+09,   4.01242718e+09,
         4.07931836e+09,   2.25817344e+09,   1.65685211e+09,
         3.69398739e+09,   4.16022443e+08,   2.70393188e+09,
         4.08810387e+09])

	for i in range(0,len(location)):
		starCam.images[len(starCam.images)-1].detectorArray[location[i]] = \
		val[i]
	plt.figure()
	plt.imshow(starCam.images[len(starCam.images)-1].detectorArray.reshape(512,512))


###############################################################################
#
#		Image 7 (Hot Pixels)
#
###############################################################################

if 0 or makeAll:
	sc.attitudeDCM = np.array([
		[ 0, 1, 0],
		[-1, 0, 0],
		[ 0, 0, 1]
		])

	msg['addBod'] = 0
	msg['rmOcc'] = 0

	msg['takeImage'] = 1
	starCam.updateState()
	msg['takeImage'] = 0
	starCam.updateState()

	location = np.array([  67238,  204676,    8502,    3127,   66954,   37175,
	        198832,  143308,  171596,  148843,   36628,  139528,
	        211838,  168193,   94000,  153682,   83715,  151847,
	        110000,   86087,  137886,  245761,  126063,   34410,
	         71959,   65053,  198715,   26613,   81186,  238110,
	        116165,   13751,  222347,   39720,   70698,   12145,
	        229304,  175423,  114801,  230443,  255514,   70676,
	        177123,   85416,   28670,   79794,   76463,   56513,
	        251937,  176083,  223210,   69272,  166395,  135999,
	        214500,   82304,  135686,   95969,  222475,  197569,
	         44762,   41205,   48567,   20308,   50656,  108693,
	        226143,  180566,   71751,  184761,  113617,  138538,
	         89182,  164795,   70876,  166511,   80919,  106320,
	        130322,   33213,  193733,   64612,   54932,  180376,
	        122599,   75537,   83099,  103414,  185201,   89593,
	        237452,  121161,   85407,  224852,  176546,    3958,
	        229694,   58432,   80575,  139825])


	for i in range(0,len(location)):
		starCam.images[len(starCam.images)-1].detectorArray[location[i]] = \
		0
	plt.figure()
	plt.imshow(starCam.images[len(starCam.images)-1].detectorArray.reshape(512,512))


pdb.set_trace()











