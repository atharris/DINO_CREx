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

def test_extended_body_lightSim():

	assert (1==1)
	assert (2==2)
	assert (3==3)


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
