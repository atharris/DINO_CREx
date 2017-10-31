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
import bodies as bod
import camera
from propagator import coe_2body
import numpy.linalg as la
from matplotlib.animation import FuncAnimation


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
print("Initialize SC: " + str(datetime.now() - start_time))

start_time = datetime.now()

qe = {}
qe['lambda'] = np.arange(420,721,5.)
#the throughput definition here is a cheat to make something
#somewhat realistic while not spending time on researching that
#realism just yet.
qe['throughput'] = (100000-(qe['lambda']-570)**2) - \
	min(100000-(qe['lambda']-570)**2)
qe['throughput'] = 0.8*qe['throughput']/max(qe['throughput'])

tc = {}
tc['lambda'] = np.arange(390,702,2.7)
tc['throughput'] = (1000*4-(tc['lambda']-545)**4) - \
	min(1000*4-(tc['lambda']-545)**4)
tc['throughput'] = 0.6*tc['throughput']/max(tc['throughput'])


cam = camera.camera(
	0.5, 				#detector_height
	0.5, 				#detector_width
	3.0, 				#focal_length
	50, 				#resolution_height
	50,					#resolution_width
	np.identity(3),  	#body2cameraDCM
	1000,			    #maximum magnitude
	qe,
	tc,
	1.7,
	sc
	)

print("Initialize Camera: " + str(datetime.now() - start_time))

msg = { 'bodies': [bod.sun,bod.earth,bod.luna, sc], 
	'rm_occ': 1, 'add_bod': 0, 'psf': 0, 'raster': 0, 
	'photon': 0, 'dark': 0, 'read': 0, 'hot_dark': 0}

states = coe_2body(msg['bodies'],0)
bod.earth.M_at_epoch += 90



start_time = datetime.now()
img = camera.image(
	np.deg2rad(180),	#alpha
	np.deg2rad(0), 	    #beta
	np.deg2rad(0), 		#gamma
	cam, 	 			#camera
	msg					#msg placeholder
	)
print("Initialize Image: " + str(datetime.now() - start_time))

if msg['raster']:
	plt.figure(figsize=(7, 7))
	im = plt.imshow(img.detector_array.reshape(50,50))
	plt.title('Rasterized Stars with PSF and Read Noise')
	plt.xlabel('Pixel')
	plt.ylabel('Line')
else:
	fig, ax = plt.subplots(figsize=(7, 7))
	ax.set_facecolor('black')
	ax.plot(img.pixel,img.line,'w.',markersize=1)
	ax.set_xlim(0, cam.resolution_width)
	ax.set_ylim(0, cam.resolution_height)
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
plt.plot(cam.lambda_set,cam.sensitivity_curve,'--',color='black',label='Sensitivity')

# plt.plot(cam.lambda_set,cam.sensitivity_curve*bb,label='Incident Light')
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
plt.plot(cam.lambda_set,cam.sensitivity_curve,'--',color='black',label='Sensitivity Curve')
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
plt.plot(cam.lambda_set,cam.sensitivity_curve,'--',color='black',label='Sensitivity Curve')
plt.title('Normalized Fluxes for Three Stars with Sensitivity Curve (Zoomed)')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Dimensionless Flux/Percent Throughput')
plt.ylim([0,1.1])
plt.legend(loc=2)


plt.figure()
plt.plot(cam.lambda_set,bb3000/max(bb5000)*cam.sensitivity_curve,'r',label='3000 K')
plt.plot(cam.lambda_set,bb4000/max(bb5000)*cam.sensitivity_curve,'g',label='4000 K')
plt.plot(cam.lambda_set,bb5000/max(bb5000)*cam.sensitivity_curve,'b',label='5000 K')
plt.title('Camera Response to Three Stars')
plt.ylabel('Dimensionless Response Function')
plt.xlabel('Wavelength (nm)')
plt.ylim([0,0.5])
plt.legend()
plt.show()
pdb.set_trace()
