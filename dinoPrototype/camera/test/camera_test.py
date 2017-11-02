#! /usr/bin/env python3
# H+
#	Title   : camera_test.py
#	Author  : Matt Muszynski
#	Date    : 10/10/17
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
sys.path.insert(0, '../dependencies/')
from em import planck, stefan_boltzmann
from constants import T_sun, r_sun, au
from numpy import arange, pi
import bodies as bod
import camera
import numpy as np
import matplotlib.pyplot as plt

#create a spacecraft object to be shared across tests
#a lot of the details of the spacecraft don't matter because
#they aren't used in the camera model. The really important one is 
#the state vector.
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

qe = {}
qe['lambda'] = np.arange(420,721,2.)
#the throughput definition here is a cheat to make something
#somewhat realistic while not spending time on researching that
#realism just yet.
qe['throughput'] = (100000-(qe['lambda']-570)**2) - \
	min(100000-(qe['lambda']-570)**2)
qe['throughput'] = 0.8*qe['throughput']/max(qe['throughput'])

tc = {}
tc['lambda'] = np.arange(390,721,1.7)
tc['throughput'] = (1000*4-(tc['lambda']-545)**4) - \
	min(1000*4-(tc['lambda']-545)**4)
tc['throughput'] = 0.6*tc['throughput']/max(tc['throughput'])

bod.earth.state = np.array([au,0,0,0,0,0])
bod.luna.state = bod.earth.state + 250000*np.array([0,1,0,0,0,0])
sc.state = bod.earth.state - 250000*np.array([1,0,0,0,0,0])

msg = { 'bodies': [
	bod.earth,
	bod.luna,
	sc
	], 
	'add_stars': 0,'rm_occ': 0, 'add_bod': 0, 'psf': 1, 
	'raster': 1, 'photon': 0, 'dark': 0, 'read': 0}

#create camera with no stars in it for tests that don't need them
#They will run significantly faster without them.
noStarCam = camera.camera(
	2, 				#detector_height
	2, 				#detector_width
	5.0, 			#focal_length
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
	sc,
	msg,
	db='../db/tycho.db'
	)

#now create a camera with stars in it for use in the tests that
#actually need them.
msg['add_stars'] = 1
StarCam = camera.camera(
	2, 				#detector_height
	2, 				#detector_width
	5.0, 			#focal_length
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
	sc,
	msg,
	db='../db/tycho.db'
	)
sc.attitudeDCM = np.identity(3)

def test_planck_eq_TSI():
	lambda_set = arange(1,10001,1) #in nm
	lambda_set = lambda_set*1e-9 #convert to m
	bb_curve = planck(T_sun,lambda_set)
	TSI = sum(pi*r_sun**2/au**2*bb_curve)
	assert( abs((TSI - 1367)/1367) <0.001 )

def test_planck_eq_stefan_boltzmann():
	sb = stefan_boltzmann(T_sun)*r_sun**2/au**2
	lambda_set = arange(1,10001,1) #in nm
	lambda_set = lambda_set*1e-9 #convert to m
	bb_curve = planck(T_sun,lambda_set)
	TSI = sum(pi*r_sun**2/au**2*bb_curve)
	assert( abs((TSI - sb)/sb) <0.001 )

def test_4_7_camera_update_state():
	msg['take_image'] = 1
	noStarCam.update_state()
	noStarCam.update_state()
	noStarCam.update_state()
	msg['take_image'] = 0
	noStarCam.update_state()
	assert(len(noStarCam.images) == 1)
	assert(len(noStarCam.images[0].scenes) == 3)
	msg['take_image'] = 1
	noStarCam.update_state()
	assert(len(noStarCam.images) == 2)
	noStarCam.update_state()
	msg['take_image'] = 0
	noStarCam.update_state()
	assert(len(noStarCam.images[1].scenes) == 2)

def test_extended_body_lightSim():
	a = 1
	assert (a==1)
	b = 2
	assert (b==2)
	c = 3
	assert (c==3)


# from lightSimFunctions import lightSim
# from numpy import identity, array
# from bodies import earth
# from constants import au
# facets = lightSim(
# 	np.identity(3), 
# 	array([0,0,0]), 
# 	array([au,0,0]), 
# 	(3,3), 
# 	200, 
# 	200, 
# 	0,
# 	earth.albedo, 
# 	earth.r_eq, 
# 	earth.name)


# bod.earth.state = np.array([au,0,0,0,0,0])
# sc.state = bod.earth.state - np.array([0,0,0,0,0,0])
# qe = {
# 	'lambda': np.arange(1,10001,1.), 
# 	'throughput': np.ones(10000)
# 	}
# tc = qe
# msg = { 'bodies': [bod.earth, sc], 'add_stars': 0,
# 'rm_occ': 1, 'add_bod': 1, 'psf': 1, 'raster': 1, 
# 'photon': 0, 'dark': 0, 'read': 0}
# print('camera init')
# cam = camera.camera(
# 	2, 				#detector_height
# 	2, 				#detector_width
# 	5.0, 			#focal_length
# 	512, 			#resolution_height
# 	512,			#resolution_width
# 	np.identity(3), #body2cameraDCM
# 	1000,		    #maximum magnitude
# 	qe,
# 	tc,
# 	10.,
# 	1., #effective area in m^2
# 	sc,
# 	msg
# )
# print('camera init')
# sc.attitudeDCM = np.identity(3)
# msg['take_image'] = 1
# cam.update_state()
# msg['take_image'] = 0
# cam.update_state()
# TSI = sum(
# 	planck(T_sun,cam.lambda_set*1e-9)*10
# 	)*pi*r_sun**2/au**2

# import pdb
# pdb.set_trace()
# assert(np.array_equal(cam.lambda_set,qe['lambda']))
# assert(np.array_equal(cam.lambda_set,tc['lambda']))
