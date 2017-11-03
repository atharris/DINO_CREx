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
import pdb
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
starCam = camera.camera(
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

def test_4_1_loadAllStars():
	#load support dict that was calculated offline
	test_4_1_support_dict = np.load('camera_test_support_files/4.1.test_support.npy')[0]
	
	#check that all sums of arrays loaded into the camera at run time
	#match what they were when the support dict was made.
	assert(sum(starCam.T) == test_4_1_support_dict['Tsum'])
	assert(sum(starCam.n1) == test_4_1_support_dict['n1sum'])
	assert(sum(starCam.n2) == test_4_1_support_dict['n2sum'])
	assert(sum(starCam.n3) == test_4_1_support_dict['n3sum'])
	assert(sum(starCam.RA) == test_4_1_support_dict['RAsum'])
	assert(sum(starCam.DE) == test_4_1_support_dict['DEsum'])
	assert(sum(starCam.VT) == test_4_1_support_dict['VTsum'])
	assert(sum(starCam.BVT) == test_4_1_support_dict['BVTsum'])

# def test_4_2_calculate_FOV():
# 	cam1 = camera.camera(
# 		1,2,2)
# 	cam2= camera.camera(
# 		3,4,5)
# 	cam3 = camera.camera(
# 		6,4,2)

# 	assert(cam1.angular_height == 15)
# 	assert(cam1.angular_width == 16)
# 	assert(cam1.angular_diagonal == 18)

# 	assert(cam2.angular_height == 54)
# 	assert(cam2.angular_width == 45)
# 	assert(cam2.angular_diagonal == 25)

# 	assert(cam2.angular_height == 13)
# 	assert(cam2.angular_width == 17)
# 	assert(cam2.angular_diagonal == 18)

def test_4_7_cameraUpdateState():
	msg['take_image'] = 1
	noStarCam.update_state()
	noStarCam.update_state()
	noStarCam.update_state()
	assert(len(noStarCam.images) == 1)
	assert(len(noStarCam.images[0].scenes) == 0)
	msg['take_image'] = 0
	noStarCam.update_state()
	assert(len(noStarCam.images) == 1)
	assert(len(noStarCam.images[0].scenes) == 3)
	msg['take_image'] = 1
	noStarCam.update_state()
	assert(len(noStarCam.images) == 2)
	assert(len(noStarCam.images[0].scenes) == 3)
	assert(len(noStarCam.images[1].scenes) == 0)
	noStarCam.update_state()
	msg['take_image'] = 0
	noStarCam.update_state()
	assert(len(noStarCam.images[1].scenes) == 2)


def test_4_9_imageRemoveOccultations():
	#enforce position of earth and location of sc.
	#this way, earth is in the exact center of the FOV
	bod.earth.state = np.array([au,0,0,0,0,0])
	sc.state = bod.earth.state - 250000*np.array([1,0,0,0,0,0])
	sc.attitudeDCM = np.identity(3)

	#take an image pointed at the earth
	#remove the stars occulted by the earth
	#but don't add the earth back in
	msg['rm_occ'] = 1
	msg['add_bod'] = 0
	msg['take_image'] = 1
	starCam.update_state()
	msg['take_image'] = 0
	starCam.update_state()
	
	#take another image, pointed in the same place, but don't
	#remove occulted stars
	msg['rm_occ'] = 0
	msg['take_image'] = 1
	starCam.update_state()
	msg['take_image'] = 0
	starCam.update_state()

	#find distance between center of FOV and each star in p/l coords.
	pixPos = starCam.images[0].scenes[0].pixel - starCam.resolution_width/2
	linePos = starCam.images[0].scenes[0].line - starCam.resolution_height/2

	#physical size of each pixel in the p/l directions
	pixSize = float(starCam.detector_width)/\
		float(starCam.resolution_width)
	lineSize = float(starCam.detector_height)/\
		float(starCam.resolution_height)

	#physical distance between center of FOV and each star in p/l coords
	pixDist = pixPos*pixSize
	lineDist = linePos*lineSize

	#physical distance between center of FOV and each star
	diagDist = np.sqrt(pixDist**2 + lineDist**2)

	#angular distance between each star and FOV center
	angDist = np.arctan2(diagDist,starCam.focal_length)

	#angular distance between center of Earth and limp
	sc2earthDist = np.linalg.norm((sc.state - bod.earth.state)[0:3])
	center2limbAng = np.arctan2(bod.earth.r_eq,sc2earthDist)

	#now find use the second image to find all the stars that were
	#removed from the first
	removedStars = starCam.images[1].scenes[0].star_ids
	ind = removedStars > -1
	for star_id in starCam.images[0].scenes[0].star_ids:
		ind = np.logical_and(ind,removedStars != star_id)

	rmPix = starCam.images[1].scenes[0].pixel[ind] - starCam.resolution_width/2
	rmLine = starCam.images[1].scenes[0].line[ind] - starCam.resolution_height/2
	rmPixDist = rmPix*pixSize
	rmLineDist = rmLine*lineSize
	rmDiagDist = np.sqrt(rmPixDist**2 + rmLineDist**2)
	rmAngDist = np.arctan2(rmDiagDist,starCam.focal_length)

	#assert that the number of stars in the image with stars removed
	#plus the number of stars removed from it equal the number of
	#stars in the image where none were removed
	assert(len((starCam.images[0].scenes[0].star_ids)) + sum(ind) \
		== len((starCam.images[1].scenes[0].star_ids)))
	#assert that star closest to the center of the FOV in the image
	#with stars removed is farther from the center than the limb of
	#the earth
	assert( min(angDist) > center2limbAng )
	#assert that the removed star farthest from the center of the FOV
	#is closer than the limb of the earth.
	assert( max(rmAngDist) < center2limbAng )


def test_4_16_pixelLineConversion():
	#find distance between center of FOV and each star in p/l coords.
	pixPos = starCam.images[1].scenes[0].pixel - starCam.resolution_width/2
	linePos = starCam.images[1].scenes[0].line - starCam.resolution_height/2

	#physical size of each pixel in the p/l directions
	pixSize = float(starCam.detector_width)/\
		float(starCam.resolution_width)
	lineSize = float(starCam.detector_height)/\
		float(starCam.resolution_height)

	pixDist = pixPos*pixSize
	lineDist = linePos*lineSize

	diagDist = np.sqrt(pixDist**2 + lineDist**2)

	#angualr distance computed by projecting physical distance
	#on focal plane through lens
	angDist = np.arctan2(diagDist,starCam.focal_length)

	#angular distance computed by dotting the camera boresight
	#vector with the unit vector to the star
	trueAngDist = np.arccos(starCam.images[1].scenes[0].c1)

	#assert that the worst error between the two angular distances
	#is still within machine precision
	assert( max(abs(angDist - trueAngDist)) < 1e-13 )


def test_4_17_planckEqStefanBoltzmann():
	sb = stefan_boltzmann(T_sun)*r_sun**2/au**2
	lambda_set = arange(1,10001,1) #in nm
	lambda_set = lambda_set*1e-9 #convert to m
	bb_curve = planck(T_sun,lambda_set)
	TSI = sum(pi*r_sun**2/au**2*bb_curve)
	assert( abs((TSI - sb)/sb) <0.001 )


def test_4_18_PlanckEqTSI():
	lambda_set = arange(1,10001,1) #in nm
	lambda_set = lambda_set*1e-9 #convert to m
	bb_curve = planck(T_sun,lambda_set)
	TSI = sum(pi*r_sun**2/au**2*bb_curve)
	assert( abs((TSI - 1367)/1367) <0.001 )

def test_4_21_mapSphere():
	import lightSimFunctions
	earthSurfaceArea = 4*np.pi*bod.earth.r_eq**2

	map100x100 = lightSimFunctions.mapSphere(100,100,bod.earth.r_eq)
	map250x250 = lightSimFunctions.mapSphere(250,250,bod.earth.r_eq)
	map500x500 = lightSimFunctions.mapSphere(500,500,bod.earth.r_eq)
	
	map100x100surfaceArea = len(map100x100[0])*map100x100[1]
	map250x250surfaceArea = len(map250x250[0])*map250x250[1]
	map500x500surfaceArea = len(map500x500[0])*map500x500[1]
	pdb.set_trace()
	assert ( (map100x100surfaceArea - earthSurfaceArea/2)/(earthSurfaceArea/2) < 1e-15 )
	assert ( (map250x250surfaceArea - earthSurfaceArea/2)/(earthSurfaceArea/2) < 1e-15 )
	assert ( (map500x500surfaceArea - earthSurfaceArea/2)/(earthSurfaceArea/2) < 1e-15 )

def test_4_23_lumos():
	import lightSimFunctions
	pdb.set_trace()
	assert ( 1 == 1 )
