#! /usr/bin/env python3
# H+
#	Title   : orbits.py
#	Author  : Matt Muszynski
#	Date    : 09/04/16
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

def tilde(threevec):
	import numpy as np
	x1, x2, x3 = threevec
	return np.array([
		[  0,-x3, x2],
		[ x3,  0,-x1],
		[-x2, x1,  0]
		])

def symbolic_tilde(vector):
	import sympy as sym
	x1 = vector[0]
	x2 = vector[1]
	x3 = vector[2]
	return sym.Matrix([
		[  0,-x3, x2],
		[ x3,  0,-x1],
		[-x2, x1,  0]
		])

def cayley_transform(matrix):
	import numpy as np
	import numpy.linalg as la
	#define identity in a generalized way so function can be used
	#in n dimensions
	ident = np.identity(matrix.shape[0])
	return (ident-matrix)*la.inv(ident+matrix)

def euler_params2dcm(beta_vec):
	import numpy as np
	b0 = beta_vec.item(0)
	b1 = beta_vec.item(1)
	b2 = beta_vec.item(2)
	b3 = beta_vec.item(3)
	return np.matrix([
		[b0**2+b1**2-b2**2-b3**2,         2*(b1*b2+b0*b3),         2*(b1*b3-b0*b2)],
		[        2*(b1*b2-b0*b3), b0**2-b1**2+b2**2-b3**2,         2*(b2*b3+b0*b1)],
		[        2*(b1*b3+b0*b2),         2*(b2*b3-b0*b1), b0**2-b1**2-b2**2+b3**2]
		])


def crp2dcp(crp_vec):
	import numpy as np
	q1 = crp_vec.item(0)
	q2 = crp_vec.item(1)
	q3 = crp_vec.item(2)

	q_vec = np.matrix([q1,q2,q3]).T

	return 1/(1+q_vec.T*q_vec).item(0)*\
	(
		(1-q_vec.T*q_vec).item(0)*np.identity(3) + \
		2*q_vec*q_vec.T - 
		2*tilde(q_vec)
	)


def Euler321_2DCM(psi,theta,phi):
	from util import r1, r2, r3
	return r1(phi).dot(r2(theta).dot(r3(psi)))

# def MRP2DCM(sigma):
# 	from numpy import array, identity
# 	# sigma_tilde = tilde(sigma)

# 	# C = identity(3) + 8*sigma_tilde.dot(sigma_tilde) - \
# 	# 4*(1-sigma.dot(sigma))*sigma_tilde/(1 + sigma.dot(sigma))**2

# 	# return C
# 	sigsq = sigma.dot(sigma)
# 	sig1 = sigma[0]
# 	sig2 = sigma[1]
# 	sig3 = sigma[2]
# 	return array([
# 		[
# 		4*(sig1**2-sig2**2-sig3**2)+(1-sigsq)**2,
# 		8*sig1*sig2 + 4*sig3*(1-sigsq),
# 		8*sig1*sig3 - 4*sig2*(1-sigsq),
# 		],
# 		[
# 		8*sig2*sig1 - 4*sig3*(1-sigsq),
# 		4*(-sig1**2+sig2**2-sig3**2)+(1-sigsq)**2,
# 		8*sig3*sig1 + 4*sig2*(1-sigsq),
# 		],
# 		[
# 		8*sig3*sig1 + 4*sig2*(1-sigsq),
# 		8*sig3*sig2 - 4*sig1*(1-sigsq),
# 		4*(-sig1**2-sig2**2+sig3**2)+(1-sigsq)**2,
# 		]
# 		])/(1+sigsq)**2

# def DCM2MRP(DCM):
# 	from numpy import sqrt, array, trace

# 	beta_0_sq = 0.25*(1 + trace(DCM))
# 	beta_1_sq = 0.25*(1 + 2*DCM[1,1] - trace(DCM))
# 	beta_2_sq = 0.25*(1 + 2*DCM[0,0] - trace(DCM))
# 	beta_3_sq = 0.25*(1 + 2*DCM[2,2] - trace(DCM))

# 	beta_0 = 0.5*sqrt(DCM[0,0] + DCM[1,1] + DCM[2,2] + 1)
# 	beta_1 = (DCM[1,2]-DCM[2,1])/4*beta_0
# 	beta_2 = (DCM[2,0]-DCM[0,2])/4*beta_0
# 	beta_3 = (DCM[0,1]-DCM[1,0])/4*beta_0

# 	sig1 = beta_1/(1+beta_0)
# 	sig2 = beta_2/(1+beta_0)
# 	sig3 = beta_3/(1+beta_0)

# 	return array([sig1,sig2,sig3])







