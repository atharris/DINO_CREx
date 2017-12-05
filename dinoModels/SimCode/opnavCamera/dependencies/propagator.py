#! /usr/bin/env python3
# H+
#	Title   : propagator.py
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
#	Usage   
#	Example	:
#	Output  :
# U-
# D+
#
# D-
###############################################################################

###############################################################################
#	coe_2body() is a quick and dirty two body propagator to demo the DINO
#		C-REx camera module without needing to connect to BSK.
#
#	Standard Dependencies:
#		Numpy
#	DINO C-REx Dependencies:
#		orbits.py: Module originally written for Steve Nerem's ASEN 5050 course
#			in Fall 2016.
#		bodies.py: Quick and dirty definition of celestial bodies to demo this
#			module without tying it to BSK.
#
#	Inputs:
#		bodies: a python list of bodies from the body.py module
#		
#	Outputs:
#		pos_vectors_in_hci:
#
###############################################################################

def coe_2body(bodies,t):
	from bodies import sun, earth, luna
	from numpy import rad2deg, concatenate, array, sqrt
	from orbits import coe2rv, anomalies

	pos_vectors_in_hci = array([])
	for body in bodies:
		if body.name == 'Sun':
			body.state = [0,0,0,0,0,0]
		else:
			mu = eval(body.central_body.lower()).mu
			central_body_state = eval(body.central_body.lower()).state
			n = sqrt(mu/body.a**3)
			M = body.M_at_epoch + rad2deg(t*n)
			nu = anomalies('M',M,body.a,body.e,silent=1)['nu']
			body.state = \
				coe2rv(body.a,body.e,body.i,body.OMEGA,body.omega,nu,mu=mu,silent=1) + \
				central_body_state
		pos_vectors_in_hci = concatenate([pos_vectors_in_hci,body.state])
	
	return pos_vectors_in_hci











