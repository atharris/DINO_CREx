#! /usr/bin/python
# /usr/bin/env python
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
import sqlite3
import pdb
import matplotlib.pyplot as plt
import camera
from constants import au
from adcs import Euler321_2DCM

print('################## Initializing ##################')

makeAll = 1


###############################################################################
#
#       Pull in canned QE and Transmission curves from DINO C-REx files
#
###############################################################################

print('############### Reading TC/QE Data ###############')

#load tranmission curve for Canon 20D
tc = np.load('tc/20D.npz')

#load QE curve for Hubble Space Telecope Advanced Camera for Surveys SITe CCD
qe = np.load('qe/ACS.npz')





###############################################################################
#
#       Initializing Beacons
#
###############################################################################

print('############## Initializing Beacons ##############')

earth = camera.beacon()
earth.r_eq = 6378.137
earth.id = 'Earth'
earth.albedo = 0.434

moon = camera.beacon()
moon.r_eq = 1738.1
moon.id = 'Earth'
moon.albedo = 0.12

earth.state = np.array([au,0,0,0,0,0])
moon.state = earth.state + 250000*np.array([0,1,0,0,0,0])
scState = earth.state - 250000*np.array([1,0,0,0,0,0])
scDCM = np.identity(3)

bodies = [earth,moon]

###############################################################################
#
#       Initializing Cameras
#
###############################################################################

print('############## Initializing Cameras ##############')

msg = {
	'addStars': 0,'rmOcc': 1, 'addBod': 1, 'psf': 1, 
	'raster': 1, 'photon': 1, 'dark': 1, 'read': 1}
takeImage = 0

#create camera with no stars in it for demos that don't need them
#They will run significantly faster without them.
cam = camera.camera(
	0.019968, 			#detector_height
	0.019968, 			#detector_width
	0.05, 				#focal_length
	512, 				#resolution_height
	512,				#resolution_width
	np.identity(3), 	#body2cameraDCM
	1000,		    	#maximum magnitude
	-1000,				#minimum magnitude (for debugging)
	qe,					#quantum efficiency curve loaded above
	tc,					#transmission efficiency curve loaded above
	1,					#lambda bin size
	0.01**2, 			#effective area in m^2
	100, 				#dark current in electrons per second
	100, 				#std for read noise in electrons
	100, 				#bin size
	2**32, 				#max bin depth
	1,					#standard deviation 
	0.001, 				#simulation timestep
	scState,			#position state of s/c
	scDCM,				#intertal 2 body DCM for s/c
	bodies,				#bodies to track in images
	takeImage,			#takeImage message
	debug=msg,			#debug message
	db='db/tycho.db'	#stellar database
	)


msg = {
	'addStars': 1,'rmOcc': 1, 'addBod': 1, 'psf': 1, 
	'raster': 1, 'photon': 1, 'dark': 1, 'read': 1}

#create camera with all stars in it for demos that do need them
#They will run significantly faster without them.
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
	10, #bin size
	2**32, #max bin depth
	1,
	0.01, 				#simulation timestep
	scState,			#position state of s/c
	scDCM,				#intertal 2 body DCM for s/c
	bodies,				#bodies to track in images
	takeImage,			#takeImage message
	debug=msg,			#debug message
	db='db/tycho.db'	#stellar database
	)


###############################################################################
#
#		Image 1 (Earth as Resolved Body, Moon as a Pt Source)
#
###############################################################################

if 0 or makeAll:
	print('##################### Image 1 ####################')

	cam.scState = np.array([au/1000,-1e6, 0, 0, 0, 0])

	earth.state = np.array([au/1000, 0, 0, 0, 0, 0])

	moon.state = np.array([au/1000, 0, 0, 0, 0, 0]) + \
		np.array([34862, -601522, 0, 0, 0 ,0])

	cam.scDCM = np.array([
		[ 0, 1, 0],
		[-1, 0, 0],
		[ 0, 0, 1]
		])
	cam.takeImage = 1
	cam.updateState()
	cam.takeImage = 0
	cam.updateState()

	plt.figure()
	plt.imshow(cam.images[len(cam.images)-1].detectorArray)
	plt.title('Image 1 (Earth as Resolved Body, Moon as a Pt Source)')

###############################################################################
#
#		Image 2 (Earth as Resolved Body, Moon as a Pt Source)
#
###############################################################################

if 0 or makeAll:
	print('##################### Image 2 ####################')

	cam.scState = np.array([au/1000,-4e6, 0, 0, 0, 0])

	earth.state = np.array([au/1000, 0, 0, 0, 0, 0])

	moon.state = np.array([au/1000, 0, 0, 0, 0, 0]) + \
		np.array([34862, -3.601522E6, 0, 0, 0 ,0])

	cam.scDCM = np.array([
		[ 0, 1, 0],
		[-1, 0, 0],
		[ 0, 0, 1]
		])
	cam.takeImage = 1
	cam.updateState()
	cam.takeImage = 0
	cam.updateState()

	plt.figure()
	plt.imshow(cam.images[len(cam.images)-1].detectorArray)
	plt.title('Image 2 (Earth as Resolved Body, Moon as a Pt Source)')

###############################################################################
#
#		Image 3 (Beacons as large resolved bodies)
#
###############################################################################

if 0 or makeAll:
	print('##################### Image 3 ####################')

	cam.scState = np.array([au/1000,-1e6, 0, 0, 0, 0])

	earth.state = np.array([au/1000, 0, 0, 0, 0, 0])

	moon.state = np.array([au/1000, 0, 0, 0, 0, 0]) + \
		np.array([34862, -9E5, 0, 0, 0 ,0])

	cam.scDCM = np.array([
		[ 0, 1, 0],
		[-1, 0, 0],
		[ 0, 0, 1]
		])
	cam.takeImage = 1
	cam.updateState()
	cam.takeImage = 0
	cam.updateState()

	plt.figure()
	plt.imshow(cam.images[len(cam.images)-1].detectorArray)
	plt.title('Image 3 (Beacons as large resolved bodies)')

###############################################################################
#
#		Image 4 (Blended Source)
#
###############################################################################

if 0 or makeAll:
	print('##################### Image 4 ####################')
	cam.scState = np.array([au/1000,-5e5, 0, 0, 0, 0])

	earth.state = np.array([au/1000, 0, 0, 0, 0, 0])

	moon.state = np.array([au/1000, 0, 0, 0, 0, 0]) + \
		np.array([-2500, 0, 4000, 0, 0 ,0])

	tmpMoonR = moon.r_eq
	moon.r_eq = earth.r_eq
	cam.scDCM = np.array([
		[ 0, 1, 0],
		[-1, 0, 0],
		[ 0, 0, 1]
		])
	cam.takeImage = 1
	cam.updateState()
	cam.takeImage = 0
	cam.updateState()

	plt.figure()
	plt.imshow(cam.images[len(cam.images)-1].detectorArray)
	plt.title('Image 4 (Blended Source)')

###############################################################################
#
#		Image 5 (Corners)
#
###############################################################################

if 0 or makeAll:
	print('##################### Image 5 ####################')
	cam.scState = np.array([au/1000,-5e5, 0, 0, 0, 0])

	earth.state = np.array([au/1000, 0, 0, 0, 0, 0]) + \
		np.array([1.025e5, 0, -1e5, 0, 0 ,0])

	moon.state = np.array([au/1000, 0, 0, 0, 0, 0]) + \
		np.array([-9.5e4, 0, 9.5e4, 0, 0 ,0])

	cam.scDCM = np.array([
		[ 0, 1, 0],
		[-1, 0, 0],
		[ 0, 0, 1]
		])
	cam.takeImage = 1
	cam.updateState()
	cam.takeImage = 0
	cam.updateState()

	plt.figure()
	plt.imshow(cam.images[len(cam.images)-1].detectorArray)
	plt.title('Image 5 (Corners)')

###############################################################################
#
#		Image 6 (Not Blended, but ROIs overlap)
#
###############################################################################

if 0 or makeAll:
	print('##################### Image 6 ####################')
	cam.scState = np.array([au/1000,-5e5, 0, 0, 0, 0])

	earth.state = np.array([au/1000, 0, 0, 0, 0, 0])

	moon.state = np.array([au/1000, 0, 0, 0, 0, 0]) + \
		np.array([-7000, 0, 4000, 0, 0 ,0])

	tmpMoonR = moon.r_eq
	moon.r_eq = earth.r_eq
	cam.scDCM = np.array([
		[ 0, 1, 0],
		[-1, 0, 0],
		[ 0, 0, 1]
		])
	cam.takeImage = 1
	cam.updateState()
	cam.takeImage = 0
	cam.updateState()

	plt.figure()
	plt.imshow(cam.images[len(cam.images)-1].detectorArray)
	plt.title('Image 6 (Not Blended, but ROIs overlap)')

###############################################################################
#
#		Image 7 (simple slew)
#
###############################################################################

if 0 or makeAll:
	print('##################### Image 7 ####################')
	msg['addBod'] = 0
	msg['rmOcc'] = 0

	starCam.takeImage = 1

	for i in range(0,10):
		alpha = i*0.1

		starCam.scDCM = Euler321_2DCM(
			np.deg2rad(alpha),
			np.deg2rad(0),
			np.deg2rad(0)
			)
		starCam.updateState()

	starCam.takeImage = 0
	starCam.updateState()

	plt.figure()
	plt.imshow(starCam.images[len(starCam.images)-1].detectorArray)
	plt.title('Image 9 (Simple Slew)')

pdb.set_trace()











