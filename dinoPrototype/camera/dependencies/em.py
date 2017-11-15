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
#	stefanBoltzmann() is a function that total flux from a star given 
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

def stefanBoltzmann(T):
	from constants import sigma
	return sigma*T**4

###############################################################################
#	photon_energy() calculates the energy for a photon of wavelength lambda.
#
#	Inputs:
#		lam: wavelength (in meters!) of photon in question. May be scalar or numpy 
#			array of wavelength values
#
#	Outputs:
#		E: energy of a single photon at wavelength lam in Joules. Will be 
#			same data type as lam (i.e. scalar if lam is scalar, numpy array if
#			lam is numpy array)
#
#	Notes:
#
###############################################################################

def photon_energy(lam):
	from constants import h, c
	return h*c/lam



