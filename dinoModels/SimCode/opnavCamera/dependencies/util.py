#! /usr/bin/env python3
# H+
#	Title   : view_tester.py
#	Author  : Matt Muszynski
#	Date    : 06/30/17
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

###############################################################################
#	block_diag() is used to create a block diagonal matrix given a set of
#	smaller matrices. It takes a numpy array in the form:
#
#			[  0   1   2   3   4   5]
#			[  6   7   8   9  10  11]
#	 		[ 12  13  14  15  16  17]
#
#	and turns it into the form:
#
#			[  0   1   2   0   0   0]
#			[  6   7   8   0   0   0]
#	 		[ 12  13  14   0   0   0]
#	 		[  0   0   0   3   4   5]
#	 		[  0   0   0   9  10  11]
#	 		[  0   0   0  15  16  17]
#
#	Inputs:
#		in_mat: an nxm numpy array consisting of m/n nxn matrices.
#
#	Outputs:
#		out_mat: an mxm array with the same nxn matrices that were passed in
#		but now buffered with zeros to make it block diagonal.
#
#	Notes:
#		This is a helper function to manipulate very large sets of data. I
#		originally wrote it so I can do many coe2rv calculations all at once.
#
###############################################################################

def block_diag(in_mat):

	from numpy import arange, zeros
	#nxm array --> mxm array (n<m)
	n = len(in_mat)
	m = len(in_mat[0])

	if m%n:
		print("Error: input array must be nxm where n divides evenly into m.")
		return

	nm_helper = arange(n*m).reshape(n,m)
	mm_helper = arange(m**2)
	m_helper = arange(m)
	out_mat = zeros(m**2)
	
	for i in range(0,n):
		for j in range(n):
			out_mat[mm_helper%((m+1)*n) == m*i+j]  = \
				in_mat[i][m_helper%n ==j]
	out_mat = out_mat.reshape(m,m)

	return out_mat

################################################################################
#	Rotation matrices
#
# 	Reference: lecture 6 ASEN 5050 CU Boulder, Fall 2016, Slide 37
#
#	Angles all in radians!
#
###############################################################################

# rx =  np.matrix( \
# [ \
# [1.,  0.,         0.        ], \
# [0.,  np.cos(theta), np.sin(theta)], \
# [0., -np.sin(theta), np.cos(theta)]  \
# ] \
# )	

def rx (theta):
	from numpy import cos, sin, zeros, arange
	try:
		rx = zeros((3,3*len(theta)))
		helper = arange(3*3*len(theta)).reshape(3,3*len(theta))
	except:
		rx = zeros((3,3))
		helper = arange(3*3).reshape(3,3)

	rx[0][helper[0]%3 == 0] =  1
	rx[1][helper[0]%3 == 1] =  cos(theta)
	rx[1][helper[0]%3 == 2] =  sin(theta)
	rx[2][helper[0]%3 == 1] = -sin(theta)
	rx[2][helper[0]%3 == 2] =  cos(theta)

	try:
		rx = rx.reshape((3,3*len(theta)))
	except:
		rx = rx.reshape(3,3)

	return rx

# ry =  np.matrix( \
# [ \
# [np.cos(theta), 0., -np.sin(theta)], \
# [0.,            1.,  0.           ], \
# [np.sin(theta), 0., np.cos(theta)]  \
# ] \
# )	

def ry(theta):
	from numpy import cos, sin, zeros, arange
	try:
		ry = zeros((3,3*len(theta)))
		helper = arange(3*3*len(theta)).reshape(3,3*len(theta))
	except:
		ry = zeros((3,3))
		helper = arange(3*3).reshape(3,3)

	ry[0][helper[0]%3 == 0] =  cos(theta)
	ry[0][helper[0]%3 == 2] = -sin(theta)
	ry[1][helper[0]%3 == 1] =  1
	ry[2][helper[0]%3 == 0] =  sin(theta)
	ry[2][helper[0]%3 == 2] =  cos(theta)

	try:
		ry = ry.reshape((3,3*len(theta)))
	except:
		ry = ry.reshape(3,3)

	return ry

# rz = np.matrix( \
# [ \
# [ np.cos(theta),   np.sin(theta), 0.], \
# [-np.sin(theta),   np.cos(theta), 0.], \
# [ 0.,                0.,             1.]  \
# ] \
# )

def rz (theta):
	from numpy import cos, sin, zeros, arange
	try:
		rz = zeros((3,3*len(theta)))
		helper = arange(3*3*len(theta)).reshape(3,3*len(theta))
	except:
		rz = zeros((3,3))
		helper = arange(3*3).reshape(3,3)

	rz[0][helper[0]%3 == 0] =  cos(theta)
	rz[0][helper[0]%3 == 1] =  sin(theta)
	rz[1][helper[0]%3 == 0] = -sin(theta)
	rz[1][helper[0]%3 == 1] =  cos(theta)
	rz[2][helper[0]%3 == 2] =  1

	try:
		rz = rz.reshape((3,3*len(theta)))
	except:
		rz = rz.reshape(3,3)
	

	return rz

def r1 (theta):
	r1 = rx(theta)
	return r1

def r2 (theta):
	r2 = ry(theta)
	return r2

def r3 (theta):
	r3 = rz(theta)
	return r3


###############################################################################
#
# interpolateLambdaDependent() 
#
#	Inputs:
#		
#	Outputs:
#
#	Notes: Please forgive me for my crappy variable names in this method.
#		at least it's short and relatively simple...
#
###############################################################################

def interpolateLambdaDependent(ex,lambda_set):
	from numpy import array
	lam = ex['lambda']
	data = ex['throughput']

	int_ex = []
	lambda_set_ex = []
	for i in range(0,len(lambda_set)):
		#if this item in lambda_set is in the lambda array passed
		#by the user, just grab its data value and use it.
		if min(abs(lambda_set[i] - lam)) < 1e-8:
			for j in range(0,len(lam)):
				if lam[j] == lambda_set[i]:
					data_ex	= data[j]
		#if this item in lambda_set is less than the minimum of the
		#lambda array passed by the user then this curve has no
		#throughput at this wavelength. Set data to zero.					
		elif lambda_set[i] < min(lam):
			data_ex = 0
		#if this item in lambda_set is greater than the maximum of the
		#lambda array passed by the user then this curve has no
		#throughput at this wavelength. Set data to zero.					
		elif lambda_set[i] > max(lam):
			data_ex = 0
		else:
		#this is the meat of this method. If this item in lambda_set is
		#not already represented by a point in the 'lambda' array passed
		#by the user, then take the point just above and just below it
		#and do a linear interpolation between them to find a representation
		#of throughput at the given lambda.
			for j in range(0,len(lam)):
				if lam[j] < lambda_set[i]:
					lower_lam = lam[j]
					lower_data = data[j]
					upper_lam = lam[j+1]
					upper_data = data[j+1]	
					m = (upper_data-lower_data)/(upper_lam-lower_lam)		
			data_ex = lower_data + m*(lambda_set[i]-lower_lam)
		lambda_set_ex.insert(len(lambda_set_ex),lambda_set[i])
		int_ex.insert(len(int_ex),data_ex)
	return {
	'lambda': array(lambda_set_ex),
	'throughput': array(int_ex)
}
