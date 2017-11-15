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


###############################################################################
#
#	Initialize Object Model
#
###############################################################################


start_time = datetime.now()
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

offset = np.deg2rad(90)

bod.earth.state = np.array([au/1000,0,0,0,0,0])
bod.luna.state = bod.earth.state - np.array(
	[bod.luna.a*np.cos(np.pi/36 + offset),bod.luna.a*np.sin(np.pi/36 + offset),0,0,0,0]
	)
sc.state = bod.earth.state - np.array(
		[1000000*np.cos(offset)
		,1000000*np.sin(offset)
		,0
		,0,0,0]
		)

print("Initialize SC: " + str(datetime.now() - start_time))

start_time = datetime.now()

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

msg = { 'bodies': [
	bod.earth,
	bod.luna,
	sc
	], 
	'addStars': 1,'rmOcc': 1, 'addBod': 0, 'psf': 1, 
	'raster': 1, 'photon': 0, 'dark': 0, 'read': 0, 'dt': 0.1}

cam = camera.camera(
	2, 				#detector_height
	2, 				#detector_width
	5.0, 			#focal_length
	512, 			#resolutionHeight
	512,			#resolutionWidth
	np.identity(3), #body2cameraDCM
	1000,		    #maximum magnitude
	-1000,			#minimum magnitude (for debugging)
	qe,
	tc,
	1,
	0.01**2, #effective area in m^2
	100, #dark current in electrons per second
	100, #std for read noise in electrons
	100, #bin size
	2**16, #max bin depth
	sc,
	msg
	)

###############################################################################
#
#	Requirement 4.1
#
###############################################################################


###############################################################################
#
#	Requirement 4.2: Lambda Function Plots
#
###############################################################################


# fig, ax = plt.subplots()
# plt.plot(qe['lambda'],qe['throughput'],label='Quantum Effeciency')
# plt.plot(tc['lambda'],tc['throughput'],label='Transmission')
# plt.title('Uninterpolated Throughput Curves')
# plt.xlabel('Wavelength (nm)')
# plt.ylabel('Throughput (%)')
# plt.xlim([390,720])
# plt.ylim([0,1])
# plt.legend()
# ax.grid(True)


# fig, ax = plt.subplots()
# plt.plot(qe['lambda'],qe['throughput'],'.',label='Quantum Effeciency')
# plt.plot(tc['lambda'],tc['throughput'],'.',label='Transmission')
# plt.legend()
# plt.title('Uninterpolated Throughput Curves')
# plt.xlabel('Wavelength (nm)')
# plt.ylabel('Throughput (%)')
# plt.xlim([390,720])
# plt.ylim([0,1])

# ax.grid(True)

# fig, ax = plt.subplots()
# plt.plot(qe['lambda'],qe['throughput'],'.',label='Quantum Effeciency')
# plt.plot(tc['lambda'],tc['throughput'],'.',label='Transmission')
# plt.legend()
# plt.title('Uninterpolated Throughput Curves')
# plt.xlabel('Wavelength (nm)')
# plt.ylabel('Throughput (%)')
# plt.xlim([495-165/2,495+165/2])
# plt.ylim([0.35,0.85])
# ax.grid(True)

# fig, ax = plt.subplots()
# plt.plot(qe['lambda'],qe['throughput'],'.',label='Quantum Effeciency')
# plt.plot(tc['lambda'],tc['throughput'],'.',label='Transmission')
# plt.legend()
# plt.title('Uninterpolated Throughput Curves')
# plt.xlabel('Wavelength (nm)')
# plt.ylabel('Throughput (%)')
# plt.xlim([460,520])
# plt.ylim([0.525,0.675])
# ax.grid(True)

# fig, ax = plt.subplots()
# plt.plot(qe['lambda'],qe['throughput'],'.',label='Quantum Effeciency')
# plt.plot(tc['lambda'],tc['throughput'],'.',label='Transmission')
# plt.legend()
# plt.title('Uninterpolated Throughput Curves')
# plt.xlabel('Wavelength (nm)')
# plt.ylabel('Throughput (%)')
# plt.xlim([490,500])
# plt.ylim([0.57,0.62])
# ax.grid(True)

# fig, ax = plt.subplots()
# plt.plot(cam.qe['lambda'],cam.qe['throughput'],'.',label='Quantum Effeciency')
# plt.plot(cam.tc['lambda'],cam.tc['throughput'],'.',label='Transmission')
# plt.title('Interpolated Throughput Curves')
# plt.xlabel('Wavelength (nm)')
# plt.ylabel('Throughput (%)')
# plt.legend()
# plt.xlim([490,500])
# plt.ylim([0.57,0.62])
# ax.grid(True)

# fig, ax = plt.subplots()
# plt.plot(cam.qe['lambda'],cam.qe['throughput'],'.',label='Quantum Effeciency')
# plt.plot(cam.tc['lambda'],cam.tc['throughput'],'.',label='Transmission')
# plt.title('Interpolated Throughput Curves')
# plt.xlabel('Wavelength (nm)')
# plt.ylabel('Throughput (%)')
# plt.legend()
# ax.grid(True)
# plt.show()

# pdb.set_trace()

###############################################################################
#
#	Requirement 4.4.8: Pointing Error
#
###############################################################################

msg = { 'bodies': [
	bod.earth,
	bod.luna,
	sc
	], 
	'addStars': 1,'rmOcc': 1, 'addBod': 1, 'psf': 1, 
	'raster': 1, 'photon': 1, 'dark': 1, 'read': 1, 'dt': 0.1}

cam = camera.camera(
	1, 				#detector_height
	1, 				#detector_width
	5.0, 			#focal_length
	512, 			#resolutionHeight
	512,			#resolutionWidth
	np.identity(3), #body2cameraDCM
	1000,		    #maximum magnitude
	-1000,			#minimum magnitude (for debugging)
	qe,
	tc,
	1,
	0.01**2, #effective area in m^2
	100, #dark current in electrons per second
	100, #std for read noise in electrons
	100, #bin size
	2**16, #max bin depth
	sc,
	msg
	)

takeImageArr = np.concatenate([np.zeros(1),np.ones(20),np.zeros(1)])

psi = 0
theta = 0
phi = 0

sc.attitudeDCM = np.identity(3)
msg['takeImage'] = 1
cam.updateState()
msg['takeImage'] = 0
cam.updateState()

for takeImage in takeImageArr:
	sc.attitudeDCM = adcs.Euler321_2DCM(psi,theta,phi)
	psi += np.deg2rad(0.05)
	msg['takeImage'] = takeImage
	cam.updateState()

for takeImage in takeImageArr:
	sc.attitudeDCM = adcs.Euler321_2DCM(psi,theta,phi)
	theta += np.deg2rad(0.05)
	msg['takeImage'] = takeImage
	cam.updateState()

for takeImage in takeImageArr:
	sc.attitudeDCM = adcs.Euler321_2DCM(psi,theta,phi)
	phi += np.deg2rad(0.15)
	msg['takeImage'] = takeImage
	cam.updateState()

takeImageArr = np.concatenate([np.zeros(1),np.ones(40),np.zeros(1)])

for takeImage in takeImageArr:
	sc.attitudeDCM = adcs.Euler321_2DCM(psi,theta,phi)
	psi += np.deg2rad(0.05)/2
	theta += np.deg2rad(0.05)/2
	phi += np.deg2rad(0.15)/2

	msg['takeImage'] = takeImage
	cam.updateState()


psi = np.deg2rad(90)
theta = 0
phi = 0
sc.attitudeDCM = adcs.Euler321_2DCM(psi,theta,phi)
msg['takeImage'] = 1
cam.updateState()
msg['takeImage'] = 0
cam.updateState()

# for takeImage in takeImageArr:
# 	sc.attitudeDCM = adcs.Euler321_2DCM(psi,theta,phi)
# 	psi += np.deg2rad(0.05)/2
# 	theta += np.deg2rad(0.05)/2
# 	phi += np.deg2rad(0.15)/2
# 	msg['takeImage'] = takeImage
# 	cam.updateState()

plt.figure()
plt.imshow(cam.images[0].detectorArray.reshape(
	cam.resolutionHeight,cam.resolutionWidth),cmap='Greys_r')
plt.title('Still Image')
plt.xlabel('Pixel')
plt.ylabel('Line')

plt.figure()
plt.imshow(cam.images[1].detectorArray.reshape(
	cam.resolutionHeight,cam.resolutionWidth),cmap='Greys_r')
plt.title('0.1 Degree Yaw')
plt.xlabel('Pixel')
plt.ylabel('Line')

plt.figure()
plt.imshow(cam.images[2].detectorArray.reshape(
	cam.resolutionHeight,cam.resolutionWidth),cmap='Greys_r')
plt.title('0.1 Degree Pitch')
plt.xlabel('Pixel')
plt.ylabel('Line')

plt.figure()
plt.imshow(cam.images[3].detectorArray.reshape(
	cam.resolutionHeight,cam.resolutionWidth),cmap='Greys_r')
plt.title('0.1 Degree Roll')
plt.xlabel('Pixel')
plt.ylabel('Line')

plt.figure()
plt.imshow(cam.images[4].detectorArray.reshape(
	cam.resolutionHeight,cam.resolutionWidth),cmap='Greys_r')
plt.title('0.1 Degree Yaw, Pitch, and Roll')
plt.xlabel('Pixel')
plt.ylabel('Line')

plt.figure()
plt.imshow(cam.images[5].detectorArray.reshape(
	cam.resolutionHeight,cam.resolutionWidth),cmap='Greys_r')
plt.title('Still Image With Visible Beacons')
plt.xlabel('Pixel')
plt.ylabel('Line')

# pdb.set_trace()
###############################################################################
#
#	Breann's Angular Size Question
#
###############################################################################

msg = { 'bodies': [
	bod.earth,
	bod.luna,
	sc
	], 
	'addStars': 1,'rmOcc': 0, 'addBod': 0, 'psf': 1, 
	'raster': 1, 'photon': 0, 'dark': 0, 'read': 0, 'dt': 0.1}

cam = camera.camera(
	1, 				#detector_height
	1, 				#detector_width
	5.0, 			#focal_length
	512, 			#resolutionHeight
	512,			#resolutionWidth
	np.identity(3), #body2cameraDCM
	1000,		    #maximum magnitude
	-1000,			#minimum magnitude (for debugging)
	qe,
	tc,
	1,
	0.01**2, #effective area in m^2
	100, #dark current in electrons per second
	100, #std for read noise in electrons
	100, #bin size
	2**16, #max bin depth
	sc,
	msg
	)


cam.msg['addStars'] = 0
cam.msg['rmOcc'] = 0
cam.msg['addBod'] = 0
cam.msg['psf'] = 1
cam.msg['raster'] = 1
cam.msg['photon'] = 1
cam.msg['dark'] = 1
cam.msg['read'] = 1
cam.msg['takeImage'] = 0
cam.msg['bodies'] = [bod.earth,sc]
offset = 0
bod.earth.state = np.array([au/1000,0,0,0,0,0])
sc.state = bod.earth.state - np.array(
		[1000000*np.cos(offset)
		,1000000*np.sin(offset)
		,0
		,0,0,0]
		)

psi = offset
theta = np.deg2rad(0)
phi = np.deg2rad(0)
attitudeDCM = adcs.Euler321_2DCM(psi,theta,phi)
sc.attitudeDCM = attitudeDCM
cam.msg['takeImage'] = 1
cam.updateState()
cam.msg['takeImage'] = 0
cam.updateState()

sc.state = bod.earth.state - np.array(
		[500000*np.cos(offset)
		,500000*np.sin(offset)
		,0
		,0,0,0]
		)

psi = offset
theta = np.deg2rad(0)
phi = np.deg2rad(0)
attitudeDCM = adcs.Euler321_2DCM(psi,theta,phi)
sc.attitudeDCM = attitudeDCM
cam.msg['takeImage'] = 1
cam.updateState()
cam.msg['takeImage'] = 0
cam.updateState()
sc.state = bod.earth.state - np.array(
		[100000*np.cos(offset)
		,100000*np.sin(offset)
		,0
		,0,0,0]
		)



psi = offset
theta = np.deg2rad(0)
phi = np.deg2rad(0)
attitudeDCM = adcs.Euler321_2DCM(psi,theta,phi)
sc.attitudeDCM = attitudeDCM
cam.msg['takeImage'] = 1
cam.updateState()
cam.msg['takeImage'] = 0
cam.updateState()

plt.figure()
plt.imshow(cam.images[0].detectorArray.reshape(cam.resolutionHeight,cam.resolutionWidth))
plt.xlabel('Pixel')
plt.ylabel('Line')
plt.title('Earth at 1e6 km. Angular Size ≈ 0.73 deg')
plt.figure()
plt.imshow(cam.images[1].detectorArray.reshape(cam.resolutionHeight,cam.resolutionWidth))
plt.xlabel('Pixel')
plt.ylabel('Line')
plt.title('Earth at 5e5 km. Angular Size ≈ 1.46 deg')
plt.figure()
plt.imshow(cam.images[2].detectorArray.reshape(cam.resolutionHeight,cam.resolutionWidth))
plt.xlabel('Pixel')
plt.ylabel('Line')
plt.title('Earth at 1e5 km. Angular Size ≈ 7.27 deg')

minPix = np.array([min(cam.images[0].scenes[0].pixel), min(cam.images[1].scenes[0].pixel), min(cam.images[2].scenes[0].pixel)])
maxPix = np.array([max(cam.images[0].scenes[0].pixel), max(cam.images[1].scenes[0].pixel), max(cam.images[2].scenes[0].pixel)])
pixDiff = maxPix - minPix
percentOfDetector = pixDiff/cam.resolutionWidth
angularSizeOnDetector = percentOfDetector*cam.angularWidth

dist2Earth = np.array([1000000, 500000, 100000])
angularDistance = np.rad2deg(np.arctan2(2*bod.earth.r_eq,dist2Earth))

pdb.set_trace()

###############################################################################
#
#	Support for Image Processing
#
###############################################################################

# msg = { 'bodies': [
# 	bod.earth,
# 	bod.luna,
# 	sc
# 	], 
# 	'addStars': 0,'rmOcc': 1, 'addBod': 1, 'psf': 1, 
# 	'raster': 1, 'photon': 0, 'dark': 0, 'read': 0, 'dt': 0.01}

# cam = camera.camera(
# 	1, 				#detector_height
# 	1, 				#detector_width
# 	5.0, 			#focal_length
# 	512, 			#resolutionHeight
# 	512,			#resolutionWidth
# 	np.identity(3), #body2cameraDCM
# 	1000,		    #maximum magnitude
# 	-1000,			#minimum magnitude (for debugging)
# 	qe,
# 	tc,
# 	1,
# 	0.01**2, #effective area in m^2
# 	100, #dark current in electrons per second
# 	100, #std for read noise in electrons
	# 100, #bin size
	# 2**16, #max bin depth
# 	sc,
# 	msg
# 	)


# #Force camera FOV tobe that of the camera in the DRM
# cam.angularHeight=10
# cam.angularWidth=12

# offsets = [
# 	np.deg2rad(0),
# 	np.deg2rad(45),
# 	np.deg2rad(90),
# 	np.deg2rad(135)
# 	]

# filename = [
# 	'0_deg.npz',
# 	'45_deg.npz',
# 	'90_deg.npz',
# 	'135_deg.npz'
# 	]

# cam.msg['addStars'] = 0
# cam.msg['rmOcc'] = 0
# cam.msg['addBod'] = 1
# cam.msg['psf'] = 1
# cam.msg['raster'] = 1
# cam.msg['photon'] = 1
# cam.msg['dark'] = 1
# cam.msg['read'] = 1
# cam.msg['takeImage'] = 0

# for i in range(0,len(offsets)):
# 	print(offset)
# 	offset = offsets[i]
# 	bod.earth.state = np.array([au/1000,0,0,0,0,0])
# 	bod.luna.state = bod.earth.state - np.array(
# 		[bod.luna.a*np.cos(np.pi/36 + offset),bod.luna.a*np.sin(np.pi/36 + offset),0,0,0,0]
# 		)
# 	sc.state = bod.earth.state - np.array(
# 			[1000000*np.cos(offset)
# 			,1000000*np.sin(offset)
# 			,0
# 			,0,0,0]
# 			)

# 	psi = offset
# 	theta = np.deg2rad(0)
# 	phi = np.deg2rad(0)
# 	attitudeDCM = adcs.Euler321_2DCM(psi,theta,phi)
# 	sc.attitudeDCM = attitudeDCM
# 	cam.msg['takeImage'] = 1
# 	cam.updateState()
# 	cam.msg['takeImage'] = 0
# 	cam.updateState()
# 	plt.figure()
# 	plt.imshow(cam.images[i].detectorArray.reshape(
# 		cam.resolutionHeight,cam.resolutionWidth))

# 	np.savez(filename[i],
# 		sun_pos = np.array([0,0,0]),
# 		earth_pos = bod.earth.state[0:3],
# 		moon_pos = bod.luna.state[0:3],
# 		sc_pos = sc.state[0:3],
# 		sc_dcm = attitudeDCM,
# 		detectorArray = cam.images[i].detectorArray
# 		)
# print('##########################################')
# pdb.set_trace()
# cam.msg['addStars'] = 1
# cam.msg['rmOcc'] = 0
# cam.msg['addBod'] = 0
# cam.msg['psf'] = 1
# cam.msg['raster'] = 1
# cam.msg['photon'] = 1
# cam.msg['dark'] = 1
# cam.msg['read'] = 1
# cam.msg['takeImage'] = 0
# #Need to reinitiate the camera so we can have stars

# cam = camera.camera(
# 	1, 				#detector_height
# 	1, 				#detector_width
# 	5.0, 			#focal_length
# 	512, 			#resolutionHeight
# 	512,			#resolutionWidth
# 	np.identity(3), #body2cameraDCM
# 	1000,		    #maximum magnitude
# 	-1000,			#minimum magnitude (for debugging)
# 	qe,
# 	tc,
# 	1,
# 	0.01**2, #effective area in m^2
# 	100, #dark current in electrons per second
# 	100, #std for read noise in electrons
	# 100, #bin size
	# 2**16, #max bin depth
# 	sc,
# 	msg
# 	)
# #Force camera FOV tobe that of the camera in the DRM
# cam.angularHeight=10
# cam.angularWidth=12


# cam.msg['takeImage'] = 1
# cam.updateState()
# msg['takeImage'] = 0
# cam.updateState()
# pdb.set_trace()
# np.savez('stars_only.npz',
# 	sun_pos = np.array([0,0,0]),
# 	sc_pos = sc.state[0:3],
# 	sc_dcm = attitudeDCM,
# 	detectorArray = cam.images[0].detectorArray
# 	)
# plt.imshow(cam.images[0].detectorArray.reshape(
# 	cam.resolutionHeight,cam.resolutionWidth))
# plt.figure()
# plt.plot(cam.images[0].scenes[0].pixel,cam.images[0].scenes[0].line,'.')
# pdb.set_trace()

###############################################################################
#
#	Verification that IP, Nav, and Cam are all doing pixel/line conversion
#	the same way.
#
###############################################################################

# sc.state = np.array([1000, 0, 0])
# bod.earth.state = np.array([1000, 1000, 0])
# #need to spoof earth radius so checkFOV doesn't fail due to
# #sc being inside beacon
# bod.earth.r_eq = 10
# bod.luna.state = np.array([1200, 1000, 400])
# bod.luna.r_eq = 10
# msg['addStars'] = 0
# msg['addBod']	= 1
# msg['bodies'] = [bod.earth,bod.luna,sc]
# sc.attitudeDCM = np.array([
# 	[0 , 1, 0],
# 	[-1, 0, 0],
# 	[ 0, 0, 1]
# 	])

# cam = camera.camera(
# 	5120, 				#detector_height
# 	5120, 				#detector_width
# 	100, 			#focal_length
# 	1024, 			#resolutionHeight
# 	1024,			#resolutionWidth
# 	np.identity(3), #body2cameraDCM
# 	1000,		    #maximum magnitude
# 	-1000,			#minimum magnitude (for debugging)
# 	qe,
# 	tc,
# 	1,
# 	0.01**2, #effective area in m^2
# 	100, #dark current in electrons per second
# 	100, #std for read noise in electrons
	# 100, #bin size
	# 2**16, #max bin depth
# 	sc,
# 	msg
# 	)

# cam.msg['takeImage'] = 1
# cam.updateState()
# cam.msg['takeImage'] = 0
# cam.updateState()

# plt.imshow(cam.images[0].detectorArray.reshape(cam.resolutionHeight,cam.resolutionWidth))
# print('Pixel: ' + str(cam.images[0].scenes[0].pixel))
# print('Line: ' + str(cam.images[0].scenes[0].line))

# pix_diff = cam.images[0].scenes[0].pixel[1] - cam.images[0].scenes[0].pixel[0]
# line_diff = cam.images[0].scenes[0].line[1] - cam.images[0].scenes[0].line[0]

# pix_diff_percent = pix_diff/cam.resolutionWidth
# line_diff_percent = line_diff/cam.resolutionHeight

# pix_ang_separation = pix_diff_percent*cam.angularWidth
# line_ang_separation = line_diff_percent*cam.angularHeight

# total_diff = np.sqrt(pix_ang_separation**2 + line_ang_separation**2)

# n_arr = np.vstack([cam.images[0].n1,cam.images[0].n2,cam.images[0].n3])
# c_arr = np.vstack([cam.images[0].c1,cam.images[0].c2,cam.images[0].c3])

# n_arr = np.transpose(n_arr)
# c_arr = np.transpose(c_arr)

# cam2earth = bod.earth.state - sc.state
# cam2moon = bod.luna.state - sc.state

# earth_moon_angle = np.rad2deg(np.arccos(cam2earth.dot(cam2moon)\
# 	/np.linalg.norm(cam2moon)\
# 	/np.linalg.norm(cam2earth)))


# bod.luna.state = np.array([1200, 1000, 450])

# cam = camera.camera(
# 	5120, 				#detector_height
# 	5120, 				#detector_width
# 	100, 			#focal_length
# 	512, 			#resolutionHeight
# 	2048,			#resolutionWidth
# 	np.identity(3), #body2cameraDCM
# 	1000,		    #maximum magnitude
# 	-1000,			#minimum magnitude (for debugging)
# 	qe,
# 	tc,
# 	1,
# 	0.01**2, #effective area in m^2
# 	100, #dark current in electrons per second
# 	100, #std for read noise in electrons
	# 100, #bin size
	# 2**16, #max bin depth
# 	sc,
# 	msg
# 	)

# cam.msg['takeImage'] = 1
# cam.updateState()
# cam.msg['takeImage'] = 0
# cam.updateState()

# plt.imshow(cam.images[0].detectorArray.reshape(cam.resolutionHeight,cam.resolutionWidth))
# print('Pixel: ' + str(cam.images[0].scenes[0].pixel))
# print('Line: ' + str(cam.images[0].scenes[0].line))

# bod.luna.state = np.array([1200, 1000, 400])

# cam = camera.camera(
# 	5120, 				#detector_height
# 	5120, 				#detector_width
# 	100, 			#focal_length
# 	512, 			#resolutionHeight
# 	2048,			#resolutionWidth
# 	np.identity(3), #body2cameraDCM
# 	1000,		    #maximum magnitude
# 	-1000,			#minimum magnitude (for debugging)
# 	qe,
# 	tc,
# 	1,
# 	0.01**2, #effective area in m^2
# 	100, #dark current in electrons per second
# 	100, #std for read noise in electrons
	# 100, #bin size
	# 2**16, #max bin depth
# 	sc,
# 	msg
# 	)

# cam.msg['takeImage'] = 1
# cam.updateState()
# cam.msg['takeImage'] = 0
# cam.updateState()

# plt.imshow(cam.images[0].detectorArray.reshape(cam.resolutionHeight,cam.resolutionWidth))
# print('Pixel: ' + str(cam.images[0].scenes[0].pixel))
# print('Line: ' + str(cam.images[0].scenes[0].line))

pdb.set_trace()
#spoof an array of take image messages.
#one scene per ms, 500 scenes = 0.5s exposure
takeImage = np.arange(100)
takeImage[takeImage < 1] = 0
takeImage[takeImage > 12] = 0
takeImage = takeImage/takeImage
takeImage[np.isnan(takeImage)] = 0

psi = offset - 0.1
theta = np.deg2rad(0)
phi = np.deg2rad(0)


for i in range(0,len(takeImage)):
	psi += np.deg2rad(0.1)
	msg['takeImage'] = takeImage[i]
	attitudeDCM = adcs.Euler321_2DCM(psi,theta,phi)
	sc.attitudeDCM = attitudeDCM
	cam.updateState()
	# theta = theta - np.deg2rad(360/10000)
	# phi = phi - np.deg2rad(360/1000)
	# psi = psi - np.deg2rad(360/1000)
	print(i)
	print(phi)

pdb.set_trace()
# plt.figure()
# plt.gca().invert_yaxis()
# for eachScene in cam.images[0].scenes:
# 	plt.plot(eachScene.pixel,eachScene.line,'.')

# plt.figure()
# plt.gca().invert_yaxis()
# for eachScene in cam.images[0].scenes:
# 	plt.plot(eachScene.psfPixel,eachScene.psfLine,'.')

plt.figure()
nonzero = cam.images[0].detectorArray[cam.images[0].detectorArray !=0]


plt.imshow(cam.images[3].detectorArray.reshape(cam.resolutionHeight,cam.resolutionWidth),cmap='Greys_r')
plt.xlabel('Pixel')
plt.ylabel('Line')
plt.title('1 Degree Yaw')
# 	colors.PowerNorm(gamma=1./2.),cmap='PuBu_r')#,
#	norm=LogNorm(vmin=min(nonzero), vmax=max(nonzero)))

pdb.set_trace()

x = np.linspace(0, cam.resolutionHeight -1, cam.resolutionHeight)
y = np.linspace(0, cam.resolutionWidth - 1,cam.resolutionWidth)
x,y = np.meshgrid(x,y)

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_surface(x, y, cam.images[0].detectorArray.reshape(
	cam.resolutionHeight,cam.resolutionWidth), 
	cmap=cm.coolwarm, linewidth=0.2, antialiased=True)
pdb.set_trace()

if msg['raster']:
	plt.figure(figsize=(7, 7))
	im = plt.imshow(img.detectorArray.reshape(50,50))
	plt.title('Rasterized Stars with PSF and Read Noise')
	plt.xlabel('Pixel')
	plt.ylabel('Line')
else:
	fig, ax = plt.subplots(figsize=(7, 7))
	ax.set_facecolor('black')
	ax.plot(img.pixel,img.line,'w.',markersize=1)
	ax.set_xlim(0, cam.resolutionWidth)
	ax.set_ylim(0, cam.resolutionHeight)
	ax.set_title('Stars with PSF')
	ax.set_xlabel('Pixel')
	ax.set_ylabel('Line')
	ax.axis("equal")

plt.figure()
plt.plot(cam.qe['lambda'],cam.qe['throughput'],label='Quantum Effeciency')
plt.plot(cam.tc['lambda'],cam.tc['throughput'],label='Transmission')

plt.title('QE and Transmission Curves')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Percent Throughput')
plt.ylim([0,1])
plt.legend()

plt.figure()
plt.plot(cam.qe['lambda'],cam.qe['throughput'],label='Quantum Effeciency')
plt.plot(cam.tc['lambda'],cam.tc['throughput'],label='Transmission')
plt.plot(cam.lambda_set,cam.sensitivityCurve,'--',color='black',label='Sensitivity')

# plt.plot(cam.lambda_set,cam.sensitivityCurve*bb,label='Incident Light')
plt.title('QE, Transmission, and Sensitivity Curves')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Percent Throughput')
plt.ylim([0,1])
plt.legend()

plt.figure()
lambda_set = np.arange(1e-3,3,1e-3)*1e-6
bb3000 = img.planck(3000,lambda_set)
bb4000 = img.planck(4000,lambda_set)
bb5000 = img.planck(5000,lambda_set)
plt.plot(lambda_set*1e9,bb3000/max(bb5000),'r',label='3000 K')
plt.plot(lambda_set*1e9,bb4000/max(bb5000),'g',label='4000 K')
plt.plot(lambda_set*1e9,bb5000/max(bb5000),'b',label='5000 K')
plt.ylim([0,1.05])
plt.title('Normalized Fluxes for Three Stars')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Dimensionless Flux')
plt.legend()

plt.figure()
plt.plot(lambda_set*1e9,bb3000/max(bb5000),'r',label='3000 K')
plt.plot(lambda_set*1e9,bb4000/max(bb5000),'g',label='4000 K')
plt.plot(lambda_set*1e9,bb5000/max(bb5000),'b',label='5000 K')
plt.plot(cam.lambda_set,cam.sensitivityCurve,'--',color='black',label='Sensitivity Curve')
plt.ylim([0,1.05])
plt.title('Normalized Fluxes for Three Stars with Sensitivity Curve')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Dimensionless Flux/Percent Throughput')
plt.legend()

plt.figure()
bb3000 = img.planck(3000,cam.lambda_set*1e-9)
bb4000 = img.planck(4000,cam.lambda_set*1e-9)
bb5000 = img.planck(5000,cam.lambda_set*1e-9)
plt.plot(cam.lambda_set,bb3000/max(bb5000),'r',label='3000 K')
plt.plot(cam.lambda_set,bb4000/max(bb5000),'g',label='4000 K')
plt.plot(cam.lambda_set,bb5000/max(bb5000),'b',label='5000 K')
plt.plot(cam.lambda_set,cam.sensitivityCurve,'--',color='black',label='Sensitivity Curve')
plt.title('Normalized Fluxes for Three Stars with Sensitivity Curve (Zoomed)')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Dimensionless Flux/Percent Throughput')
plt.ylim([0,1.1])
plt.legend(loc=2)


plt.figure()
plt.plot(cam.lambda_set,bb3000/max(bb5000)*cam.sensitivityCurve,'r',label='3000 K')
plt.plot(cam.lambda_set,bb4000/max(bb5000)*cam.sensitivityCurve,'g',label='4000 K')
plt.plot(cam.lambda_set,bb5000/max(bb5000)*cam.sensitivityCurve,'b',label='5000 K')
plt.title('Camera Response to Three Stars')
plt.ylabel('Dimensionless Response Function')
plt.xlabel('Wavelength (nm)')
plt.ylim([0,0.5])
plt.legend()
plt.show()
pdb.set_trace()
