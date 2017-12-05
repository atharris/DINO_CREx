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
from em import planck, stefanBoltzmann
from constants import T_sun, r_sun, au
from numpy import arange, pi
import bodies as bod
import camera
import numpy as np
import matplotlib.pyplot as plt
import pdb
from adcs import Euler321_2DCM
import lightSimFunctions
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
#somewhat realistic it is suffient for the test that need to be done here.
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
scState = bod.earth.state - 250000*np.array([1,0,0,0,0,0])
scDCM = np.identity(3)

bodies = [bod.earth,bod.luna]
msg = {
	'addStars': 0,'rmOcc': 0, 'addBod': 0, 'psf': 1, 
	'raster': 1, 'photon': 0, 'dark': 0, 'read': 0, 'dt': 0.01}
takeImage = 0

#create camera with no stars in it for tests that don't need them
#They will run significantly faster without them.
noStarCam = camera.camera(
	2, 					#detectorHeight
	2, 					#detectorWidth
	5.0, 				#focalLength
	512, 				#resolutionHeight
	512,				#resolutionWidth
	np.identity(3), 	#body2cameraDCM
	1000,		   	 	#maximum magnitude
	-1000,				#minimum magnitude (for debugging)
	qe,					#quantum efficiency dictionary
	tc,					#transmission curve dictionary
	1,					#wavelength bin size in nm
	0.01**2, 			#effective area in m^2
	100, 				#dark current in electrons per second
	100, 				#std for read noise in electrons
	100, 				#bin size
	2**32, 				#max bin depth
	1,					#sigma for gaussian psf
	0.01, 				#simulation timestep
	scState,			#position state of s/c
	scDCM,				#intertal 2 body DCM for s/c
	bodies,				#bodies to track in images
	takeImage,			#takeImage message
	debug=msg,			#debug message
	db='../db/tycho.db'	#stellar database
	)

# #now create a camera with stars in it for use in the tests that
# #actually need them.
msg['addStars'] = 1
starCam = camera.camera(
	2, 					#detectorHeight
	2, 					#detectorWidth
	5.0, 				#focalLength
	512, 				#resolutionHeight
	512,				#resolutionWidth
	np.identity(3), 	#body2cameraDCM
	1000,		   	 	#maximum magnitude
	-1000,				#minimum magnitude (for debugging)
	qe,					#quantum efficiency dictionary
	tc,					#transmission curve dictionary
	1,					#wavelength bin size in nm
	0.01**2, 			#effective area in m^2
	100, 				#dark current in electrons per second
	100, 				#std for read noise in electrons
	100, 				#bin size
	2**32, 				#max bin depth
	1,					#sigma for gaussian psf
	0.01, 				#simulation timestep
	scDCM,				#intertal 2 body DCM for s/c
	scState,			#position state of s/c
	bodies,				#bodies to track in images
	takeImage,			#takeImage message
	debug=msg,			#debug message
	db='../db/tycho.db'	#stellar database
	)
#create a camera with a tiny detector so we can find just a single
#star for doing PSF tests.

tinyCam = camera.camera(
	0.02,				#detectorHeight
	0.02,				#detectorWidth
	5.0, 				#focalLength
	512, 				#resolutionHeight
	512,				#resolutionWidth
	np.identity(3), 	#body2cameraDCM
	1000,		   	 	#maximum magnitude
	-1000,				#minimum magnitude (for debugging)
	qe,					#quantum efficiency dictionary
	tc,					#transmission curve dictionary
	1,					#wavelength bin size in nm
	0.01**2, 			#effective area in m^2
	100, 				#dark current in electrons per second
	100, 				#std for read noise in electrons
	100, 				#bin size
	2**32, 				#max bin depth
	1,					#sigma for gaussian psf
	0.01, 				#simulation timestep
	scState,			#position state of s/c
	scDCM,				#intertal 2 body DCM for s/c
	bodies,				#bodies to track in images
	takeImage,			#takeImage message
	debug=msg,			#debug message
	db='../db/tycho.db'	#stellar database
	)

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

# # # def test_4_2_calculate_FOV():
# # # 	cam1 = camera.camera(
# # # 		1,2,2)
# # # 	cam2= camera.camera(
# # # 		3,4,5)
# # # 	cam3 = camera.camera(
# # # 		6,4,2)

# # # 	assert(cam1.angularHeight == 15)
# # # 	assert(cam1.angularWidth == 16)
# # # 	assert(cam1.angularDiagonal == 18)

# # # 	assert(cam2.angularHeight == 54)
# # # 	assert(cam2.angularWidth == 45)
# # # 	assert(cam2.angularDiagonal == 25)

# # # 	assert(cam2.angularHeight == 13)
# # # 	assert(cam2.angularWidth == 17)
# # # 	assert(cam2.angularDiagonal == 18)

def test_4_7_cameraUpdateState():
	assert(len(noStarCam.images) == 0)
	noStarCam.takeImage = 1
	noStarCam.updateState()
	noStarCam.updateState()
	noStarCam.updateState()
	assert(len(noStarCam.images) == 1)
	assert(len(noStarCam.images[0].scenes) == 0)
	noStarCam.takeImage = 0
	noStarCam.updateState()
	assert(len(noStarCam.images) == 1)
	assert(len(noStarCam.images[0].scenes) == 3)
	noStarCam.takeImage = 1
	noStarCam.updateState()
	assert(len(noStarCam.images) == 2)
	assert(len(noStarCam.images[0].scenes) == 3)
	assert(len(noStarCam.images[1].scenes) == 0)
	noStarCam.updateState()
	noStarCam.takeImage = 0
	noStarCam.updateState()
	assert(len(noStarCam.images) == 2)
	assert(len(noStarCam.images[0].scenes) == 3)
	assert(len(noStarCam.images[1].scenes) == 2)

def test_4_8_findStarsInFOV():
	msg = { 'bodies': [
		bod.earth,
		bod.luna
		], 
		'addStars': 1,'rmOcc': 0, 'addBod': 0, 'psf': 1, 
		'raster': 1, 'photon': 0, 'dark': 0, 'read': 0}

	OriCam = camera.camera(
		1.5,				#detectorHeight
		1.5, 				#detectorWidth
		3.0, 				#focalLength
		512, 				#resolutionHeight
		512,				#resolutionWidth
		np.identity(3), 	#body2cameraDCM
		3.57,		   	 	#maximum magnitude
		-1000,				#minimum magnitude (for debugging)
		qe,					#quantum efficiency dictionary
		tc,					#transmission curve dictionary
		1,					#wavelength bin size in nm
		0.01**2, 			#effective area in m^2
		100, 				#dark current in electrons per second
		100, 				#std for read noise in electrons
		100, 				#bin size
		2**32, 				#max bin depth
		1,					#sigma for gaussian psf
		0.01,				#integration timestep
		scState,			#position state of s/c
		scDCM,				#intertal 2 body DCM for s/c
		bodies,				#bodies to track in images
		takeImage,			#takeImage message
		debug=msg,			#debug message
		db='../db/tycho.db'	#stellar database
		)
	OriCam.scDCM = Euler321_2DCM(
		np.deg2rad(85),
		np.deg2rad(0),
		np.deg2rad(0)
		)
	
	OriCam.takeImage = 1
	OriCam.updateState()
	OriCam.takeImage = 0
	OriCam.updateState()


	UMiCam = camera.camera(
		2.3,				#detectorHeight
		2.3,				#detectorWidth
		4.0, 				#focalLength
		512, 				#resolutionHeight
		512,				#resolutionWidth
		np.identity(3), 	#body2cameraDCM
		4,			   	 	#maximum magnitude
		-1000,				#minimum magnitude (for debugging)
		qe,					#quantum efficiency dictionary
		tc,					#transmission curve dictionary
		1,					#wavelength bin size in nm
		0.01**2, 			#effective area in m^2
		100, 				#dark current in electrons per second
		100, 				#std for read noise in electrons
		100, 				#bin size
		2**32, 				#max bin depth
		1,					#sigma for gaussian psf
		0.01,				#integration timestep
		scState,			#position state of s/c
		scDCM,				#intertal 2 body DCM for s/c
		bodies,				#bodies to track in images
		takeImage,			#takeImage message
		debug=msg,			#debug message
		db='../db/tycho.db'	#stellar database
		)
	UMiCam.scDCM = Euler321_2DCM(
		np.deg2rad(187),
		np.deg2rad(59),
		np.deg2rad(0)
		)

	UMiCam.takeImage = 1
	UMiCam.updateState()
	UMiCam.takeImage = 0
	UMiCam.updateState()


	plot = 1

	if plot:
		plt.figure()
		plt.imshow(OriCam.images[0].detectorArray.reshape(
			OriCam.resolutionHeight,
			OriCam.resolutionWidth
			))
		plt.figure()
		plt.imshow(UMiCam.images[0].detectorArray.reshape(
			UMiCam.resolutionHeight,
			UMiCam.resolutionWidth
			))

		plt.figure()
		plt.plot(OriCam.images[0].scenes[0].pixel,OriCam.images[0].scenes[0].line,'.')
		plt.xlim(0,OriCam.resolutionWidth)
		plt.ylim(0,OriCam.resolutionHeight)
		plt.gca().invert_yaxis()
		plt.axis('equal')
		plt.figure()
		plt.plot(UMiCam.images[0].scenes[0].pixel,UMiCam.images[0].scenes[0].line,'.')
		plt.xlim(0,UMiCam.resolutionWidth)
		plt.ylim(0,UMiCam.resolutionHeight)
		plt.gca().invert_yaxis()
		plt.axis('equal')

	savedData = np.load('camera_test_support_files/4.8.test_support.npy')
	savedUMiRA = savedData[0]
	savedUMiDE = savedData[1]
	savedOriRA = savedData[2]
	savedOriDE = savedData[3]

	assert(np.array_equal(UMiCam.images[0].RA,savedUMiRA))
	assert(np.array_equal(UMiCam.images[0].DE,savedUMiDE))
	assert(np.array_equal(OriCam.images[0].RA,savedOriRA))
	assert(np.array_equal(OriCam.images[0].DE,savedOriDE))



def test_4_9_imageRemoveOccultations():
	#enforce position of earth and location of sc.
	#this way, earth is in the exact center of the FOV
	bod.earth.state = np.array([au,0,0,0,0,0])
	starCam.scState = bod.earth.state - 250000*np.array([1,0,0,0,0,0])
	starCam.scDCM = np.identity(3)

	#take an image pointed at the earth
	#remove the stars occulted by the earth
	#but don't add the earth back in
	msg['rmOcc'] = 1
	msg['addBod'] = 0
	starCam.takeImage = 1
	starCam.updateState()
	starCam.takeImage = 0
	starCam.updateState()
	
	#take another image, pointed in the same place, but don't
	#remove occulted stars
	msg['rmOcc'] = 0
	starCam.takeImage = 1
	starCam.updateState()
	starCam.takeImage = 0
	starCam.updateState()

	#find distance between center of FOV and each star in p/l coords.
	pixPos = starCam.images[0].scenes[0].pixel - starCam.resolutionWidth/2
	linePos = starCam.images[0].scenes[0].line - starCam.resolutionHeight/2

	#physical size of each pixel in the p/l directions
	pixSize = float(starCam.detectorWidth)/\
		float(starCam.resolutionWidth)
	lineSize = float(starCam.detectorHeight)/\
		float(starCam.resolutionHeight)

	#physical distance between center of FOV and each star in p/l coords
	pixDist = pixPos*pixSize
	lineDist = linePos*lineSize

	#physical distance between center of FOV and each star
	diagDist = np.sqrt(pixDist**2 + lineDist**2)

	#angular distance between each star and FOV center
	angDist = np.arctan2(diagDist,starCam.focalLength)

	#angular distance between center of Earth and limb
	sc2earthDist = np.linalg.norm((starCam.scState - bod.earth.state)[0:3])
	center2limbAng = np.arctan2(bod.earth.r_eq,sc2earthDist)

	#assert that star closest to the center of the FOV in the image
	#with stars removed is farther from the center than the limb of
	#the earth
	assert( min(angDist) > center2limbAng )


	#now find use the second image to find all the stars that were
	#removed from the first
	removedStars = starCam.images[1].scenes[0].starIDs
	ind = removedStars > -1
	for starID in starCam.images[0].scenes[0].starIDs:
		ind = np.logical_and(ind,removedStars != starID)

	rmPix = starCam.images[1].scenes[0].pixel[ind] - starCam.resolutionWidth/2
	rmLine = starCam.images[1].scenes[0].line[ind] - starCam.resolutionHeight/2
	rmPixDist = rmPix*pixSize
	rmLineDist = rmLine*lineSize
	rmDiagDist = np.sqrt(rmPixDist**2 + rmLineDist**2)
	rmAngDist = np.arctan2(rmDiagDist,starCam.focalLength)

	#assert that the number of stars in the image with stars removed
	#plus the number of stars removed from it equal the number of
	#stars in the image where none were removed
	assert(len((starCam.images[0].scenes[0].starIDs)) + sum(ind) \
		== len((starCam.images[1].scenes[0].starIDs)))
	#assert that the removed star farthest from the center of the FOV
	#is closer than the limb of the earth.
	assert( max(rmAngDist) < center2limbAng )
	plot = 0
	if plot:

		plt.figure()
		plt.plot(starCam.images[0].scenes[0].pixel, starCam.images[0].scenes[0].line, '.', markersize='2' )
		plt.plot(256,256,marker='x',color='red',markersize='5')
		plt.xlim([0,512])
		plt.ylim([0,512])
		plt.xlabel('Pixel')
		plt.ylabel('Line')
		plt.title('Test 4.9 Full Star Field with Occulted Stars Removed')
		plt.axes().set_aspect('equal')

		plt.figure()
		plt.plot(starCam.images[0].scenes[0].pixel, starCam.images[0].scenes[0].line, '.', markersize='5' )
		plt.plot(256,256,marker='x',color='red',markersize='5')
		plt.xlim([156,356])
		plt.ylim([156,356])
		plt.xlabel('Pixel')
		plt.ylabel('Line')
		plt.title('Test 4.9 Partial Star Field with Occulted Stars Removed')
		plt.axes().set_aspect('equal')

		plt.figure()
		plt.plot(rmPix+256, rmLine+256, '.', markersize='5' )
		plt.plot(256,256,marker='x',color='red',markersize='5')
		plt.xlim([156,356])
		plt.ylim([156,356])
		plt.xlabel('Pixel')
		plt.ylabel('Line')
		plt.title('Test 4.9 Partial Star Field with Non-Occulted Stars Removed')
		plt.axes().set_aspect('equal')

		plt.figure()
		plt.plot(starCam.images[1].scenes[0].pixel, starCam.images[1].scenes[0].line, '.', markersize='2' )
		plt.plot(256,256,marker='x',color='red',markersize='5')
		plt.xlim([0,512])
		plt.ylim([0,512])
		plt.xlabel('Pixel')
		plt.ylabel('Line')
		plt.title('Test 4.9 Full Star Field')
		plt.axes().set_aspect('equal')

		plt.figure()
		plt.plot(rmPix+256, rmLine+256, '.', markersize=2)
		plt.plot(256,256,marker='x',color='red',markersize='5')
		plt.xlim([0,512])
		plt.ylim([0,512])
		plt.xlabel('Pixel')
		plt.ylabel('Line')
		plt.title('Test 4.9 Full Star Field with Non-Occulted Stars Removed')
		plt.axes().set_aspect('equal')
		pdb.set_trace()

# set up image for tests 10-12
tinyCam.scDCM = Euler321_2DCM(
	np.deg2rad(1.12551889),
	np.deg2rad(2.26739556),
	np.deg2rad(0)
	)
tinyCam.takeImage = 1
tinyCam.updateState()
tinyCam.takeImage = 0
tinyCam.updateState()

def test_4_10_PSF_sum_to_I():
	assert ( sum(tinyCam.images[0].scenes[0].psfI) - tinyCam.images[0].scenes[0].I < 1e-10)

def test_4_11_PSF_mean_at_256_256():
	assert ( abs(256 - np.mean(tinyCam.images[0].scenes[0].psfPixel)) < 1e-10 )
	assert ( abs(256 - np.mean(tinyCam.images[0].scenes[0].psfLine)) < 1e-10 )

def test_4_12_PSF_sigma_is_correct():
	import matplotlib.mlab as mlab

	#flatten into 2d PDF
	r = np.sqrt((256 - tinyCam.images[0].scenes[0].psfPixel)**2 + \
		(256 - tinyCam.images[0].scenes[0].psfLine)**2)

	pixel = tinyCam.images[0].scenes[0].psfPixel
	line = tinyCam.images[0].scenes[0].psfLine
	psfI = tinyCam.images[0].scenes[0].psfI

	#344.12312454 term accounts for the fact that the two gaussian curves
	#are not normalized in the same way. This test really only cares
	#that sigma is the same on both, so the fact that we have to spoof the
	#amplitude doesn't really bother us.
	assert(
		max(abs(
			344.12312454*mlab.bivariate_normal(pixel-256, line-256, sigmax=1.0, sigmay=1.0, mux=0.0, muy=0.0, sigmaxy=0.0)-psfI
			)) < 1e-9
			)

def test_4_13_rasterize():
	#need to take an image so we can access image.photonEnergy()
	tinyCam.takeImage = 1
	tinyCam.updateState()
	tinyCam.takeImage = 0
	tinyCam.updateState()
	pixPos = np.array([ 1.4,  4.9,  4.9,  4.4,  1.4,  4.2,  4.1,  1.9,  1.3,  3.2])
	linePos = np.array([ 4.2,  3.1,  2.4,  3.5,  2.8,  3.2,  2.4,  2.5,  4.1,  2.2])
	I = np.array([ 2.,  1.,  2.,  3.,  1.,  4.,  0.,  3.,  2.,  3.])
	rasterize = tinyCam.images[0].rasterize(5,5,pixPos,linePos,I)
	offlineRasterize = np.array(
		[ 
		0,  0,  0,  0,  0,  
		0,  0,  0,  0,  0,  
		0,  4,  0,  3,  2, 
		0,  0,  0,  0,  8,  
		0,  4,  0,  0,  0
		])
	# 1,4 2+2=4
	# 4,3 1+3+4=8
	# 4,2 2+0=2
	# 1,2 1+3=4
	# 3,2 3=3
	assert( np.array_equal(offlineRasterize,rasterize) )



def test_4_14_photon_energy():
	#need to take an image so we can access image.photonEnergy()
	tinyCam.takeImage = 1
	tinyCam.updateState()
	tinyCam.takeImage = 0
	tinyCam.updateState()

	from constants import h, c
	testValues = np.array(
		[ 5960.,  6032.,  8337.,  8759.,  1998.,  
		9567.,   633.,  3347., 5504.,  5353.])

	runTimeValues = tinyCam.images[0].photonEnergy(testValues)

	energiesComputedOffline = np.array(
		[  3.33296279e-29,   3.29317942e-29,   2.38268661e-29, 2.26789111e-29, 9.94217129e-29,   
		2.07635186e-29, 3.13814506e-28,   5.93500396e-29,   3.60909488e-29, 3.71090197e-29])
	assert (sum(abs(energiesComputedOffline - runTimeValues)) < 1e12)

def test_4_16_pixelLineConversion():
	#find distance between center of FOV and each star in p/l coords.
	pixPos = starCam.images[1].scenes[0].pixel - starCam.resolutionWidth/2
	linePos = starCam.images[1].scenes[0].line - starCam.resolutionHeight/2

	#physical size of each pixel in the p/l directions
	pixSize = float(starCam.detectorWidth)/\
		float(starCam.resolutionWidth)
	lineSize = float(starCam.detectorHeight)/\
		float(starCam.resolutionHeight)

	pixDist = pixPos*pixSize
	lineDist = linePos*lineSize

	diagDist = np.sqrt(pixDist**2 + lineDist**2)

	#angualr distance computed by projecting physical distance
	#on focal plane through lens
	angDist = np.arctan2(diagDist,starCam.focalLength)

	#angular distance computed by dotting the camera boresight
	#vector with the unit vector to the star
	trueAngDist = np.arccos(starCam.images[1].scenes[0].c1)

	plot = 0

	if plot:
		plt.figure()
		plt.plot(starCam.images[1].scenes[0].pixel, starCam.images[1].scenes[0].line, '.', markersize='2' )
		plt.plot(256,256,marker='x',color='red',markersize='5')
		plt.xlim([0,512])
		plt.ylim([0,512])
		plt.xlabel('Pixel')
		plt.ylabel('Line')
		plt.title('Test 4.10 Star Field')
		plt.axes().set_aspect('equal')

	#assert that the worst error between the two angular distances
	#is still within machine precision	
	assert( max(abs(angDist - trueAngDist)) < 1e-13 )


def test_4_17_planckEqStefanBoltzmann():
	sb = stefanBoltzmann(T_sun)*r_sun**2/au**2
	lambdaSet = arange(1,10001,1) #in nm
	lambdaSet = lambdaSet*1e-9 #convert to m
	bbCurve = planck(T_sun,lambdaSet)
	TSI = sum(pi*r_sun**2/au**2*bbCurve)
	assert( abs((TSI - sb)/sb) <0.001 )


def test_4_18_PlanckEqTSI():
	lambdaSet = arange(1,10001,1) #in nm
	lambdaSet = lambdaSet*1e-9 #convert to m
	bbCurve = planck(T_sun,lambdaSet)
	TSI = sum(pi*r_sun**2/au**2*bbCurve)
	assert( abs((TSI - 1367)/1367) <0.001 )

def test_4_20_checkFOV():
	#remove moon from bodies message
	noStarCam.bodies = [bod.earth]
	msg['addBod'] = 1

	#position earth and sc so earth is at the center of the FOV
	bod.earth.state = np.array([au/1000,0,0,0,0,0])
	noStarCam.scState = bod.earth.state - np.array([100000,0,0,0,0,0])
	noStarCam.takeImage = 1
	noStarCam.updateState()
	noStarCam.takeImage = 0
	noStarCam.updateState()

	#position earth and sc so earth is completelt out of the FOV
	bod.earth.state = np.array([au/1000,0,0,0,0,0])
	noStarCam.scState = bod.earth.state - np.array([0,100000,0,0,0,0])
	noStarCam.takeImage = 1
	noStarCam.updateState()
	noStarCam.takeImage = 0
	noStarCam.updateState()

	#position earth and sc so earth is completelt out of the FOV
	bod.earth.state = np.array([au/1000,0,0,0,0,0])
	noStarCam.scState = bod.earth.state - np.array([100000,23000,0,0,0,0])
	noStarCam.takeImage = 1
	noStarCam.updateState()
	noStarCam.takeImage = 0
	noStarCam.updateState()

	#position earth and sc so earth is completelt out of the FOV
	bod.earth.state = np.array([au/1000,0,0,0,0,0])
	noStarCam.scState = bod.earth.state - np.array([100000,0,23000,0,0,0])
	noStarCam.takeImage = 1
	noStarCam.updateState()
	noStarCam.takeImage = 0
	noStarCam.updateState()

	plot = 1
	if plot:
		plt.figure()
		plt.imshow(noStarCam.images[2].detectorArray.reshape(
			noStarCam.resolutionWidth,
			noStarCam.resolutionHeight), cmap='Greys_r')
		plt.xlabel('Pixel')
		plt.ylabel('Line')
		plt.title('Body Fully In View')

		plt.figure()
		plt.imshow(noStarCam.images[3].detectorArray.reshape(
			noStarCam.resolutionWidth,
			noStarCam.resolutionHeight), cmap='Greys_r')
		plt.xlabel('Pixel')
		plt.ylabel('Line')
		plt.title('Body Fully Out of View')

		plt.figure()
		plt.imshow(noStarCam.images[4].detectorArray.reshape(
			noStarCam.resolutionWidth,
			noStarCam.resolutionHeight), cmap='Greys_r')
		plt.xlabel('Pixel')
		plt.ylabel('Line')
		plt.title('Body Center Out of View (Pixel Dimension)')

		plt.figure()
		plt.imshow(noStarCam.images[5].detectorArray.reshape(
			noStarCam.resolutionWidth,
			noStarCam.resolutionHeight), cmap='Greys_r')
		plt.xlabel('Pixel')
		plt.ylabel('Line')
		plt.title('Body Center Out of View (Line Dimension)')

	assert(len(noStarCam.images[2].scenes[0].pixel) != 0)
	assert(len(noStarCam.images[3].scenes[0].pixel) == 0)
	assert(len(noStarCam.images[4].scenes[0].pixel) != 0)
	assert(len(noStarCam.images[5].scenes[0].pixel) != 0)

def test_4_21_mapSphere():
	import lightSimFunctions
	earthSurfaceArea = 4*np.pi*bod.earth.r_eq**2

	map100x100 = lightSimFunctions.mapSphere(100,100,bod.earth.r_eq)
	map250x250 = lightSimFunctions.mapSphere(250,250,bod.earth.r_eq)
	map500x500 = lightSimFunctions.mapSphere(500,500,bod.earth.r_eq)
	
	map100x100surfaceArea = len(map100x100[0])*map100x100[1]
	map250x250surfaceArea = len(map250x250[0])*map250x250[1]
	map500x500surfaceArea = len(map500x500[0])*map500x500[1]
	assert ( (map100x100surfaceArea - earthSurfaceArea/2)/(earthSurfaceArea/2) < 1e-15 )
	assert ( (map250x250surfaceArea - earthSurfaceArea/2)/(earthSurfaceArea/2) < 1e-15 )
	assert ( (map500x500surfaceArea - earthSurfaceArea/2)/(earthSurfaceArea/2) < 1e-15 )

def test_4_22_mapSphere():
	import lightSimFunctions
	earthSurfaceArea = 4*np.pi*bod.earth.r_eq**2

	map100x100 = lightSimFunctions.mapSphere(10,10,bod.earth.r_eq)
	map250x250 = lightSimFunctions.mapSphere(250,250,bod.earth.r_eq)
	map500x500 = lightSimFunctions.mapSphere(500,500,bod.earth.r_eq)
	
	lon100x100 = map100x100[0][:,1]
	lat100x100 = map100x100[0][:,0]
	lon250x250 = map100x100[0][:,1]
	lat250x250 = map100x100[0][:,0]
	lon500x500 = map100x100[0][:,1]
	lat500x500 = map100x100[0][:,0]

	map100x100surfaceArea = len(map100x100[0])*map100x100[1]
	map250x250surfaceArea = len(map250x250[0])*map250x250[1]
	map500x500surfaceArea = len(map500x500[0])*map500x500[1]

	assert ( sum(lat100x100 + np.flip(lat100x100,0)) < 1e-8 )
	assert ( sum(lat250x250 + np.flip(lat250x250,0)) < 1e-8 )
	assert ( sum(lat500x500 + np.flip(lat500x500,0)) < 1e-8 )
	assert ( sum(np.flip(lon100x100,0)+lon100x100 -180) < 1e-8 )
	assert ( sum(np.flip(lon250x250,0)+lon250x250 -180) < 1e-8 )
	assert ( sum(np.flip(lon500x500,0)+lon500x500 -180) < 1e-8 )

def test_4_23_lumos():
	import lightSimFunctions
	earthLumos = lightSimFunctions.lumos(
		np.array([au,0,0]),
		np.array([0,0,0]),
		bod.earth.albedo,
		bod.earth.r_eq,
		100,
		100)
	assert ( abs(earthLumos[2]*len(earthLumos[1]) - 2*np.pi*bod.earth.r_eq**2) < 10e-7 )


# def test_4_23_lumos():
# 	earthLumos = lightSimFunctions.lumos(
# 		np.array([au,0,0]),
# 		np.array([0,0,0]),
# 		bod.earth.albedo,
# 		bod.earth.r_eq,
# 		100,
# 		100)
# 	assert ( abs(earthLumos[2]*len(earthLumos[1]) - 2*np.pi*bod.earth.r_eq**2) < 10e-7 )

def test_4_24_facetAreaCamview_eq_pi_r_sq():
	earthLightSim = lightSimFunctions.lightSim(
		np.identity(3),
		np.array([0,0,0]),
		np.array([au,0,0]),
		(30,30),
		500,500,
		1,
		bod.earth.albedo,
		bod.earth.r_eq,
		'earth')
	assert (
		abs(
			sum(earthLightSim['facetArea'])[0] - np.pi*bod.earth.r_eq**2
			)/(np.pi*bod.earth.r_eq**2) < 0.02 
		)


# def test_4_x_lightSim():
# 	from lightSimFunctionsOld import lightSim
# 	earthLightSim = lightSim(
# 		np.identity(3),
# 		np.array([0,0,0]),
# 		np.array([au,0,0]),
# 		(20,20),
# 		100,
# 		100,
# 		True,
# 		bod.earth.albedo,
# 		bod.earth.r_eq,
# 		'Earth'
# 		)
# 	pdb.set_trace()	
# 	assert( 1 == 1 )
