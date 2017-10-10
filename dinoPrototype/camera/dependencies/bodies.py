#! /usr/bin/env python3
# H+
#	Title   : view_tester.py
#	Author  : Matt Muszynski
#	Date    : 06/30/17
#	Synopsis: Functions for orbit stuff
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

import numpy as np
print_data = 0

class sc:
	def __init__(
		self, 
		name,
		central_body, 
		coe_epoch,
		a,
		e,
		i,
		omega,
		OMEGA,
		M_at_epoch,
		state,
		nu,
		I
		):
		self.name = name
		self.central_body = central_body
		self.coe_epoch = coe_epoch
		self.a = a
		self.e = e
		self.i = i
		self.omega = omega
		self.OMEGA = OMEGA
		self.M_at_epoch = M_at_epoch
		self.state = state
		self.nu = nu
		self.I = I


class body:
	def __init__(
		self, 
		name,
		central_body, 
		mu, 
		coe_epoch,
		a,
		e,
		i,
		omega,
		OMEGA,
		M_at_epoch,
		r_eq,
		r_pole,
		h,
		p_rot,
		state,
		nu
		):
		self.name = name
		self.central_body = central_body
		self.mu = mu
		self.coe_epoch = coe_epoch
		self.a = a
		self.e = e
		self.i = i
		self.omega = omega
		self.OMEGA = OMEGA
		self.M_at_epoch = M_at_epoch
		self.r_eq = r_eq
		self.r_pole = r_pole
		self.h = h
		self.p_rot = p_rot
		self.state = state
		self.nu = nu
		self.surface = self.create_surface(r_eq,r_pole)

	def create_surface(self, r_eq, r_pole):
		from numpy import mgrid, pi, sin, cos
		grid = mgrid[0.1:2*pi:100j, 0.1:pi:100j]
		theta = grid[0].reshape(1,10000)[0]
		phi = grid[1].reshape(1,10000)[0]

		return {
		'n1': r_eq*sin(phi)*cos(theta),
		'n2': r_eq*sin(phi)*sin(theta),
		'n3': r_pole*cos(phi)
		}



###############################################################################
#
#	Sun
#
###############################################################################

sun = body(
	"Sun", #name of this body
	"None", #central body
	1.32712428e11, #mu in km^3/s^2
	np.nan, #Epoch of coe presented here, in JD
	np.nan, #a in km
	np.nan, # e
	np.nan, #inclination in deg
	np.nan, #omega in deg
	np.nan, #OMEGA in deg
	np.nan, #Mean Anomaly at epoch in deg
	695700, #equitorial radius in km
	695700, #polar radius in km
	np.nan, #scale height in km
	np.nan, #sidereal rotation period in s
	np.array([]), #state vector in HCI frame. [x,y,z,v_x,v_y,v_z] 
	np.nan #true anomaly measured at same time as state vector
	)
###############################################################################
#
#	Mercury
#
###############################################################################

mercury = body(
	"Mercury", #name of this body
	"Sun", #central body
	2.2032e4, #mu in km^3/s^2
	2451545.0, #Epoch of coe presented here, in JD
	57909083, #a in km
	0.205631752, # e
	7.00498625, #inclination in deg
	77.45611904, #omega in deg
	48.33089304, #OMEGA in deg
	252.25090551, #Mean Anomaly at epoch in deg
	2439.0, #equitorial radius in km
	np.nan, #polar radius in km
	np.nan, #scale height in km
	np.nan, #sidereal rotation period in s
	np.array([]), #state vector in HCI frame. [x,y,z,v_x,v_y,v_z] 
	np.array([]) #true anomaly measured at same time as state vector
	)
###############################################################################
#
#	Venus
#
###############################################################################

venus = body(
	"Venus", #name of this body
	"Sun", #central body
	3.257e5, #mu in km^3/s^2
	2451545.0, #Epoch of coe presented here, in JD
	108208601, #a in km
	0.006771882, # e
	3.39446619, #inclination in deg
	131.56370724, #omega in deg
	76.67992019, #OMEGA in deg
	181.97980084, #Mean Anomaly at epoch in deg
	6052.0, #equitorial radius in km
	np.nan, #polar radius in km
	np.nan, #scale height in km
	np.nan, #sidereal rotation period in s
	np.array([]), #state vector in HCI frame. [x,y,z,v_x,v_y,v_z] 
	np.array([]) #true anomaly measured at same time as state vector
	)

###############################################################################
#
#	Earth
#
###############################################################################

earth = body(
	"Earth", #name of this body
	"Sun", #central body
	3.986004415e5, #mu in km^3/s^2
	2451545.0, #Epoch of coe presented here, in JD
	149598023, #a in km
	0.01671022, # e
	0, #inclination in deg
	102.93734808, #omega in deg
	0, #OMEGA in deg
	100.46644851, #Mean Anomaly at epoch in deg
	6378.137, #equitorial radius in km
	6356.752, #polar radius in km
	8.5, #scale height in km
	23.9345*3600,  #sidereal rotation period in s
	np.array([]), #state vector in HCI frame. [x,y,z,v_x,v_y,v_z] 
	np.array([]) #true anomaly measured at same time as state vector
	)

###############################################################################
#
#	Luna
#
###############################################################################

luna = body(
	"Luna", #name of this body
	"Earth", #central body
	4902.799, #mu in km^3/s^2
	2451545.0, #Epoch of coe presented here, in JD
	384400, #a in km
	0.0554, # e
	5.145396, #inclination in deg
	318.15, #omega in deg
	125.08, #OMEGA in deg
	135.37, #Mean Anomaly at epoch in deg
	1738.1, #equitorial radius in km
	1736.0, #polar radius in km
	np.nan, #scale height in km
	np.nan, #sidereal rotation period in s
	np.array([]), #state vector in HCI frame. [x,y,z,v_x,v_y,v_z] 
	np.array([]) #true anomaly measured at same time as state vector
	)

###############################################################################
#
#	Mars
#
###############################################################################

mars = body(
	"Mars", #name of this body
	"Sun", #central body
	4.305e4, #mu in km^3/s^2
	2451545.0, #Epoch of coe presented here, in JD
	227939186, #a in km
	0.09341233, # e
	1.85061, #inclination in deg
	336.04084, #omega in deg
	49.57854, #OMEGA in deg
	355.45332, #Mean Anomaly at epoch in deg
	3397.2, #equitorial radius in km
	np.nan, #polar radius in km
	np.nan, #scale height in km
	np.nan, #sidereal rotation period in s
	np.array([]), #state vector in HCI frame. [x,y,z,v_x,v_y,v_z] 
	np.array([]) #true anomaly measured at same time as state vector
	)

###############################################################################
#
#	Jupiter
#
###############################################################################

jupiter = body(
	"Jupiter", #name of this body
	"Sun", #central body
	1.268e8, #mu in km^3/s^2
	2451545.0, #Epoch of coe presented here, in JD
	778298361, #a in km
	0.04849851, # e
	1.30326966, #inclination in deg
	14.331309, #omega in deg
	100.46444064, #OMEGA in deg
	34.35148392, #Mean Anomaly at epoch in deg
	71492.0, #equitorial radius in km
	np.nan, #polar radius in km
	np.nan, #scale height in km
	np.nan, #sidereal rotation period in s
	np.array([]), #state vector in HCI frame. [x,y,z,v_x,v_y,v_z] 
	np.array([]) #true anomaly measured at same time as state vector
	)

###############################################################################
#
#	Saturn
#
###############################################################################

saturn = body(
	"Saturn", #name of this body
	"Sun", #central body
	3.794e7, #mu in km^3/s^2
	2451545.0, #Epoch of coe presented here, in JD
	1429394133, #a in km
	0.055508622, # e
	2.48887810, #inclination in deg
	93.05678728, #omega in deg
	113.66442370, #OMEGA in deg
	50.07747138, #Mean Anomaly at epoch in deg
	60268.0, #equitorial radius in km
	np.nan, #polar radius in km
	np.nan, #scale height in km
	np.nan, #sidereal rotation period in s
	np.array([]), #state vector in HCI frame. [x,y,z,v_x,v_y,v_z] 
	np.array([]) #true anomaly measured at same time as state vector
	)

# ###############################################################################
# #
# #	Uranus
# #
# ###############################################################################

# uranus = body(
# 	"Uranus", #name of this body
# 	"Sun", #central body
# 	0.00490e6, #mu in km^3/s^2
# 	2451545.0, #Epoch of coe presented here, in JD
# 	384400, #a in km
# 	0.0554, # e
# 	5.16, #inclination in deg
# 	318.15, #omega in deg
# 	125.08, #OMEGA in deg
# 	135.37, #Mean Anomaly at epoch in deg
# 	1738.1, #equitorial radius in km
# 	1736.0, #polar radius in km
# 	np.nan, #scale height in km
# 	np.nan, #sidereal rotation period in s
# 	np.array([]), #state vector in HCI frame. [x,y,z,v_x,v_y,v_z] 
# 	np.array([]) #true anomaly measured at same time as state vector
# 	)

# ###############################################################################
# #
# #	Neptune
# #
# ###############################################################################

# neptune = body(
# 	"Neptune", #name of this body
# 	"Sun", #central body
# 	0.00490e6, #mu in km^3/s^2
# 	2451545.0, #Epoch of coe presented here, in JD
# 	384400, #a in km
# 	0.0554, # e
# 	5.16, #inclination in deg
# 	318.15, #omega in deg
# 	125.08, #OMEGA in deg
# 	135.37, #Mean Anomaly at epoch in deg
# 	1738.1, #equitorial radius in km
# 	1736.0, #polar radius in km
# 	np.nan, #scale height in km
# 	np.nan, #sidereal rotation period in s
# 	np.array([]), #state vector in HCI frame. [x,y,z,v_x,v_y,v_z] 
# 	np.array([]) #true anomaly measured at same time as state vector
# 	)

# ###############################################################################
# #
# #	Pluto
# #
# ###############################################################################

# pluto = body(
# 	"Pluto", #name of this body
# 	"Sun", #central body
# 	0.00490e6, #mu in km^3/s^2
# 	2451545.0, #Epoch of coe presented here, in JD
# 	384400, #a in km
# 	0.0554, # e
# 	5.16, #inclination in deg
# 	318.15, #omega in deg
# 	125.08, #OMEGA in deg
# 	135.37, #Mean Anomaly at epoch in deg
# 	1738.1, #equitorial radius in km
# 	1736.0, #polar radius in km
# 	np.nan, #scale height in km
# 	np.nan, #sidereal rotation period in s
# 	np.array([]), #state vector in HCI frame. [x,y,z,v_x,v_y,v_z] 
# 	np.array([]) #true anomaly measured at same time as state vector
# 	)


###############################################################################
#
#	Print celestial body data. I wrote this to test that I was using
#		classes right. Left it in as a test.
#
###############################################################################


if print_data == 1:
	for obj in gc.get_objects():
	    if isinstance(obj, body):
	        print(obj.name)
	        print("\tCentral Body: " + str(obj.central_body))
	        print("\tGravitational Parameter: " + str(obj.mu) + " km^3/s^2")
	        print("\tSemimajor Axis: " + str(obj.a) + " km")
	        print("\tEccentricity: " + str(obj.e))
	        print("\tInclination: " + str(obj.i) + " Degrees")
	        print("\tArgument of Perigee: " + str(obj.omega) + " Degrees")
	        print("\tRAAN: " + str(obj.OMEGA) + " Degrees")
	        print("\tMean Anomaly at Epoch: " + str(obj.M_at_epoch) + " Degrees")
	        print("\tRadius: " + str(obj.r) + " km")
	        print("\tScale Height: " + str(obj.h) + " km")
	        print("\tRotation Period: " + str(obj.p_rot) + " s")






