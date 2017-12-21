#! /usr/bin/env python3
# H+
#	Title   : hot_dark_test.py
#	Authors : Matt Muszynski
#				Ishaan Patel
#	Date    : 11/01/17
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
import pdb
import random

random.seed(777)

#create a spacecraft object to be shared across tests
#a lot of the details of the spacecraft don't matter because
#they aren't used in the camera model. The really important one is
#the state vector.
sc = bod.sc(
	"SC", 			# name of this body
	"Earth", 		# central body
	2451545.0, 		# Epoch of coe presented here, in JD
	122164 , 		# a in km
	0, 				# e
	0, 				# inclination in deg
	0, 				# omega in deg
	0, 				# OMEGA in deg
	0, 				# Mean Anomaly at epoch in deg
	np.array([]), 	# state vector in HCI frame. [x,y,z,v_x,v_y,v_z]
	np.nan, 		# true anomaly measured at same time as state vector
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
	'raster': 1, 'photon': 0, 'dark': 0, 'read': 0, 'hot_dark': 0}

# Initialize Hot Dark Array
init_hot_dark = np.ones( (1, 262144) )
sc.attitudeDCM = np.identity(3)

def test_calculate_hot_dark():
	########################################
	## 	Initialize Cameras
	########################################
	cam1 = camera.camera(
        2.0,            # detector_height
        2.0,            # detector_width
        5.0,            # focal_length
        512,            # resolution_height
        512,            # resolution_width
        np.identity(3), # body2cameraDCM
        1000,           # maximum magnitude
        -1000,          # minimum magnitude (for debugging)
        qe,
        tc,
        1,
        0.01**2, 		# effective area in m^2
        100, 			# dark current in electrons per second
        100, 			# std for read noise in electrons
		init_hot_dark,	# hot_dark_cam1
        sc,
        msg,
        db='../db/tycho.db'
        )

    	cam2= camera.camera(
 		3.0,
        4.0,
        7.0,
        512,
        512,
        np.identity(3),
        1000,
        -1000,
        qe,
        tc,
        1,
        0.01**2,
		100, 			# dark current in electrons per second
		100, 			# std for read noise in electrons
		init_hot_dark,  # hot_dark_cam2,
		sc,
        msg,
        db='../db/tycho.db'
        )

        cam3 = camera.camera(
        6.0,
        4.0,
        7.0,
        512,
        512,
        np.identity(3),
        1000,
        -1000,
        qe,
        tc,
        1,
        0.01**2,
        100, 			# dark current in electrons per second
        100, 			# std for read noise in electrons
		init_hot_dark,  # hot_dark_cam3,
        sc,
        msg,
        db='../db/tycho.db'
        )

########################################
## 	Take Two Images without Hot_dark
########################################
	msg['hot_dark'] = 0
# 	First Image
	msg['take_image'] = 1
	cam1.update_state()
	cam2.update_state()
	cam3.update_state()
	msg['take_image'] = 0
	cam1.update_state()
	cam2.update_state()
	cam3.update_state()
# 	Second Image
	msg['take_image'] = 1
	cam1.update_state()
	cam2.update_state()
	cam3.update_state()
	msg['take_image'] = 0
	cam1.update_state()
	cam2.update_state()
	cam3.update_state()
########################################
## 	Take Two Images with Hot_Dark
########################################
	msg['hot_dark'] = 1
# 	Image Three (+ Hot_Dark)
	msg['take_image'] = 1
	cam1.update_state()
	cam2.update_state()
	cam3.update_state()
	msg['take_image'] = 0
	cam1.update_state()
	cam2.update_state()
	cam3.update_state()
# 	Image Four (+ Hot_Dark)
	msg['take_image'] = 1
	cam1.update_state()
	cam2.update_state()
	cam3.update_state()
	msg['take_image'] = 0
	cam1.update_state()
	cam2.update_state()
	cam3.update_state()

########################################
## 	Camera 1 Properties
########################################
	# Image 1
	cam1_img1 = cam1.images[0].detector_array
	cam11_mean = np.mean(cam1_img1)
	cam11_var  = np.var(cam1_img1)
	# Image 2
	cam1_img2 = cam1.images[1].detector_array
	cam12_mean = np.mean(cam1_img2)
	cam12_var  = np.var(cam1_img2)
	# Image 3
	cam1_img3 = cam1.images[2].detector_array
	cam13_mean = np.mean(cam1_img3)
	cam13_var = np.var(cam1_img3)
	# Image 4
	cam1_img4 = cam1.images[3].detector_array
	cam14_mean = np.mean(cam1_img4)
	cam14_var = np.var(cam1_img4)

	pdb.set_trace()
	# Images without Hot_Dark should be different from those with Hot_dark
	assert ( cam1.images[0].detector_array.any != cam1.images[2].detector_array.any )
	assert ( cam1.images[1].detector_array.any != cam1.images[3].detector_array.any )

	assert ( cam2.images[0].detector_array.any != cam2.images[2].detector_array.any )
	assert ( cam2.images[1].detector_array.any != cam2.images[3].detector_array.any )

	assert ( cam3.images[0].detector_array.any != cam3.images[2].detector_array.any )
	assert ( cam3.images[1].detector_array.any != cam3.images[3].detector_array.any )


	# First Two Images from each camera should have the same MEAN & VARIANCE since NO Hot-Dark Pixels added
	assert ( abs(np.mean(cam1.images[0].detector_array) - np.mean(cam1.images[1].detector_array) ) <= 1e-12 )
	assert ( abs(np.mean(cam2.images[0].detector_array) - np.mean(cam2.images[1].detector_array) ) <= 1e-12 )
	assert ( abs(np.mean(cam3.images[0].detector_array) - np.mean(cam3.images[1].detector_array) ) <= 1e-12 )

	assert ( abs(np.std(cam1.images[0].detector_array) - np.std(cam1.images[1].detector_array) ) <= 1e-12 )
	assert ( abs(np.std(cam2.images[0].detector_array) - np.std(cam2.images[1].detector_array) ) <= 1e-12 )
	assert ( abs(np.std(cam3.images[0].detector_array) - np.std(cam3.images[1].detector_array) ) <= 1e-12 )


	# Third Image MEAN & VARIANCE should be different from First Image MEAN & VARIANCE
	assert ( abs(np.mean(cam1.images[2].detector_array) - np.mean(cam1.images[0].detector_array) ) >= 0.1 )
	assert ( abs(np.mean(cam2.images[2].detector_array) - np.mean(cam2.images[0].detector_array) ) >= 0.1 )
	assert ( abs(np.mean(cam3.images[2].detector_array) - np.mean(cam3.images[0].detector_array) ) >= 0.1 )

	assert ( abs(np.std(cam1.images[2].detector_array) - np.std(cam1.images[0].detector_array) ) >= 0.1 )
	assert ( abs(np.std(cam2.images[2].detector_array) - np.std(cam2.images[0].detector_array) ) >= 0.1 )
	assert ( abs(np.std(cam3.images[2].detector_array) - np.std(cam3.images[0].detector_array) ) >= 0.1 )


	# Fourth Image MEAN & VARIANCE should be different from SECOND Image MEAN & VARIANCE
	assert ( abs(np.mean(cam1.images[3].detector_array) - np.mean(cam1.images[1].detector_array) ) >= 0.1 )
	assert ( abs(np.mean(cam2.images[3].detector_array) - np.mean(cam2.images[1].detector_array) ) >= 0.1 )
	assert ( abs(np.mean(cam3.images[3].detector_array) - np.mean(cam3.images[1].detector_array) ) >= 0.1 )

	assert ( abs(np.std(cam1.images[3].detector_array) - np.std(cam1.images[1].detector_array) ) >= 0.1 )
	assert ( abs(np.std(cam2.images[3].detector_array) - np.std(cam2.images[1].detector_array) ) >= 0.1 )
	assert ( abs(np.std(cam3.images[3].detector_array) - np.std(cam3.images[1].detector_array) ) >= 0.1 )
