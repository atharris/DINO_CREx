#! /usr/bin/env python3
# H+
#	Title   : em.py
#	Author  : Matt Muszynski
#	Date    : 09/11/17
#	Synopsis: Functions for modeling E&M Phenomena
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
#	Usage   
#	Example	:
#	Output  :
# U-
# D+
#
# D-
###############################################################################
# import pdb
# import matplotlib.pyplot as plt

###############################################################################
#	planck() is a function that calculates a planck blackbody function
#
#	Inputs:
#		T: temperature of the star in Kelvin
#		lam: wavelength bins at which to calculate in METERS(!)
#
#	Outputs:
#		I: intensity of light at each wavelengh bin in W/m^2/nm/sr
#
#	Notes:
#
###############################################################################

def planck(T,lam):
	from constants import h, c, k_B
	from numpy import exp

	top_part = 2*h*c**2
	bottom_part = lam**5*(exp(h*c/(lam*k_B*T))-1)

	I = top_part/bottom_part*1e-9 #1e-9 to convert to per nm

	return I

###############################################################################
#	stefan_boltzmann() is a function that total flux from a star given 
#
#	Inputs:
#		T: temperature of the star in Kelvin
#
#	Outputs:
#		F: total flux at stellar surface in W/m^2
#
#	Notes:
#
###############################################################################

def stefan_boltzmann(T):
	from constants import sigma
	return sigma*T**4


#debugging stuff. Removed this once I had to import it into camera.py
#
# from numpy import arange, pi
# from constants import au, r_sun, sigma

# from datetime import datetime
# start_time = datetime.now()

# T = 5778
# lam = arange(0.1,10,10/250)*1e-6

# for i in range(0,300000):
# 	I = planck(T,lam)
# print(datetime.now() - start_time)
# print("Integrated Planck Function")
# #1/10 nm bin size, 1/pi accounts for sr^-1 in the planck function output.
# print(sum(I)*r_sun**2/au**2/pi/10) 
# print("Stephan-Boltzmann Flux")
# print(stefan_boltzmann(T)*r_sun**2/au**2)


