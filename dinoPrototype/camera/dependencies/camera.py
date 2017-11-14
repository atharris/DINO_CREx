#! /usr/bin/env python3
# H+
#	Title   : camera.py
#	Author  : Matt Muszynski
#	Date    : 07/26/17
#	Synopsis: Functions for modeling a spacecraft camera
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
#	Camera is a class for the simulation of a spacecraft camera.
#
#	methods:
#		load_all_stars()
#		calculate_FOV()
#		interpolate_lambda()
#		find_lambda_set()
#		update_state()
#
#	User Variables:
#		detector_height: Physical height of detector
#			Can be any dimension so long as it matches detector_width and
#			focal_length
#		detector_width: Physical width of detector
#			Can be any dimension so long as it matches detector_height and
#			focal_length
#		focal_length: Physical distance between lens and focal point
#			Can be any dimension so long as it matches detector_width and
#			detector_width
#		resolution_height: # of pixels in same dimension as detector_height
#		resolution_width: # of pixels in same dimension as detector_width
#		angular_height: Angular height of the detector field of view. 
#			Calculated using calculate_fov()
#		angular_width: Angular width of the detector field of view. 
#			Calculated using calculate_fov()
#		body2cameraDCM: DCM describing orientation of camera relative to sc
#		max_mag: Maximum magnitude star visible to this camera.
#		qe: Quantum Effeciency Curve as a dict of two numpy arrays. 
#			One must be called 'lambda', and the other must be called 
#			must be called 'throughput'. qe['lambda'] must be an array of 
#			wavelength bins, and qe['throughput'] must be an array of 
#			percentages detected at the given 'lambda'. qe['lambda'] must
#			be in nm, and qe['throughput'] should be a number between 0 and
#			1. qe['lambda'] and tc['lambda'] do not necessarily need to
#			match as interpolate_lambda will interpolate qe['throughput']
#			so that they do.
#		tc: Transmission Curve as a dict of two numpy arrays. 
#			One must be called 'lambda', and the other must be called 
#			must be called 'throughput'. tc['lambda'] must be an array of 
#			wavelength bins, and tc['throughput'] must be an array of 
#			percentages transmitted at the given 'lambda'. tc['lambda'] must
#			be in nm, and tc['throughput'] should be a number between 0 and
#			1. qe['lambda'] and tc['lambda'] do not necessarily need to
#			match as interpolate_lambda will interpolate tc['throughput']
#			so that they do.
#		min_lambda_bin_size: the width of bins in wavelength space to be
#			returned from interpolate_lambda().
#		sc: Spacecraft object this camera is attached to
#		hot_dark: NOT YET IMPLEMENTED!
#
#	Computed Variables:
#		angular_height: Size of detector field of view on the sky.
#		angular_width: 
#		RA: Right ascension of all stars loaded from db
#		DE: Declination of all stars loaded from db
#		n1: First unit vector component in inertial space of all stars loaded 
#			from db via load_all_stars()
#		n2: Second unit vector component in inertial space of all stars loaded 
#			from db via load_all_stars()
#		n3: Third unit vector component in inertial space of all stars loaded 
#			from db via load_all_stars()
#		VT: Tycho Visual band magnitude of all stars loaded from db
#		BVT: Tycho (B-V) color index of all stars loaded from db
#
#	Keyword Arguments:
#		verbose: print more to standard out for debugging. Mostly timings.
#		db: File path to a stellar database file. If not provided, 
#			db/tycho_small.db is used.
#
###############################################################################
class camera:
	def __init__(
		self, 
		detector_height, 
		detector_width, 
		focal_length,
		resolution_height,
		resolution_width,
		body2cameraDCM,
		max_mag,
		min_mag,
		qe,
		tc,
		lambda_bin_size,
		effective_area,
		dark_current,
		read_sigma,
		dnBinSize,
		dnDepthMax,
		sc,
		msg,
		**kwargs
		):

		from em import planck
		from numpy import pi, array

		try:
			verbose = kwargs['verbose']
		except:
			verbose = 0

		try: 
			db = kwargs['db']
		except:
			db = 'db/tycho.db'
		
		FOV = self.calculate_FOV(
			focal_length,
			detector_height, 
			detector_width,
			verbose=verbose
			)

		if msg['add_stars']:
			allstars = self.load_all_stars(db, max_mag, min_mag, verbose=verbose)
		lambda_set = self.find_lambda_set(qe, tc, lambda_bin_size)
		qe = self.interpolate_lambda_dependent(qe,lambda_set)
		tc = self.interpolate_lambda_dependent(tc,lambda_set) 
		sensitivity_curve = qe['throughput']*tc['throughput']
		lambda_set = lambda_set[sensitivity_curve != 0]

		if msg['add_stars']:
			self.RA = allstars['RA']
			self.DE = allstars['DE']
			self.n1 = allstars['n1']
			self.n2 = allstars['n2']
			self.n3 = allstars['n3']
			self.VT = allstars['VT']
			self.BVT = allstars['BVT']
			self.star_id = allstars['star_id']
			self.T = allstars['T']
			self.reduction_term = allstars['reduction_term']
		else:
			self.RA = array([])
			self.DE = array([])
			self.n1 = array([])
			self.n2 = array([])
			self.n3 = array([])
			self.VT = array([])
			self.BVT = array([])
			self.star_id = array([])
			self.T = array([])
			self.reduction_term = array([])

		self.solar_bb = planck(5778,lambda_set*1e-9)
		self.lambda_set = lambda_set
		self.sensitivity_curve = sensitivity_curve[sensitivity_curve != 0]
		self.qe = qe
		self.tc = tc
		self.lambda_bin_size = lambda_bin_size
		self.resolution_height = resolution_height
		self.resolution_width = resolution_width		
		self.body2cameraDCM = body2cameraDCM
		self.max_mag = max_mag
		self.min_mag = min_mag
		self.angular_height = FOV['alpha']
		self.angular_width = FOV['beta']
		self.angular_diagonal = FOV['gamma']
		self.detector_height = detector_height
		self.detector_width = detector_width
		self.focal_length = focal_length
		self.read_sigma = read_sigma
		self.sc = sc
		self.effective_area = effective_area
		self.dark_current = dark_current
		self.images = {}
		self.msg = msg
		self.dnBinSize = dnBinSize
		self.dnDepthMax = dnDepthMax

	###########################################################################
	#	load_all_stars()
	#
	#		Inputs:
	#			none
	#		Outputs:
	#			allstars: A dict with RA, DE, inertial position vector, 
	#				Johnson V magnitude and Johnson (B-V) color index for
	#				all stars in the given database
	#			db: string. The database to pull star data from. Format is 
	#				in the DINO C-REx Image Generation documentation
	#
	#		Keyword Arguements:
	#			verbose: boolean. If set to 1, load_all_stars will present
	#				the user with timings for the database call
	#
	###########################################################################
	def load_all_stars(self,db, max_mag, min_mag, **kwargs):
		import sqlite3
		from numpy import array, sin, cos, deg2rad, logical_and

		try: 
			verbose = kwargs['verbose']
			start_time = datetime.now()
		except:
			verbose = 0

		conn = sqlite3.connect(db)
		c = conn.cursor()

		select_string = "SELECT RA, DE, VTmag, BTmag, id, reduction_term, " + \
			"computed_temperature from tycho_data"
		c.execute(select_string)

		if verbose:
			print("DB Query: " + str(datetime.now() - start_time))
			start_time = datetime.now()

		RA = []
		DE = []
		VT = []
		BT = []
		phi = []
		theta = []
		star_id = []
		reduction_term = []
		T = []

		for row in c:
			RA.append(float(row[0]))
			DE.append(float(row[1]))
			#this is NOT a cheat. This is the official way for
			#us to deal with the fact that Tycho-2 doesn't have
			#BT or VT for some stars ATM!
			try:
				VT.append(float(row[2]))
			except:
				VT.append(float(row[3]))
			try:
				BT.append(float(row[3]))
			except:
				BT.append(float(row[2]))
			theta.append(float(RA[-1]))
			phi.append(float(90-DE[-1]))
			star_id.append(float(row[4]))
			#this IS a cheat. force the flux from the star
			#to be zero if we haven't gotten a computed T in for
			#it yet. The temp doesn't matter if the reduction
			#term is zero
			try:
				reduction_term.append(float(row[5]))
			except:
				reduction_term.append(float(0))
			try:
				T.append(float(row[6]))
			except:
				T.append(float(5000))

		RA = array(RA)
		DE = array(DE)
		VT = array(VT)
		BT = array(BT)
		BVT = BT - VT
		phi = array(phi)
		theta = array(theta)
		star_id = array(star_id)
		reduction_term = array(reduction_term)
		T = array(T)

		ind = VT < max_mag
		ind = logical_and(ind, VT > min_mag)
		RA = RA[ind]
		DE = DE[ind]
		VT = VT[ind]
		BVT = BVT[ind]
		theta = theta[ind]
		phi = phi[ind]
		star_id = star_id[ind]
		T = T[ind]
		reduction_term = reduction_term[ind]
		#convert spherical coordinates to cartesian in inertial
		n1 = sin(deg2rad(phi))*cos(deg2rad(theta))
		n2 = sin(deg2rad(phi))*sin(deg2rad(theta))
		n3 = cos(deg2rad(phi))

		return {
			'RA': RA,
			'DE': DE,
			'VT': VT,
			'BVT': BVT,
			'n1': n1,
			'n2': n2,
			'n3': n3,
			'star_id': star_id,
			'reduction_term': reduction_term,
			'T': T	
			}

	###############################################################################
	#
	# calculate_FOV() finds the angular size of a camera FOV using freshman-level
	#	optics. 
	#
	#	Inputs:
	#		focal_length: Focal length of modeled lens 
	#		detector_height: y-direction dimension of the lens
	#		detector_width: x-direction dimension of the lens
	#		
	#	Outputs:
	#		A dict with the following entries:
	#		alhpa: angular size of FOV along y-direction (height) 
	#		beta: angular size of FOV along x-direction (width) 
	#
	#	Notes:
	#		Assumes lens is thin. Units for detector dimensions and focal length
	#		can be anything, but they must all match.
	#
	###############################################################################

	def calculate_FOV(self, focal_length, detector_height, detector_width, **kwargs):
		from numpy import sqrt, arctan2, rad2deg
		f = focal_length
		a = detector_height
		b = detector_width
		c = sqrt(a**2 + b**2)

		#angular distance of diagonal of FOV
		gamma = 2*arctan2(c/2,f)
		alpha = 2*arctan2(a/2,f)
		beta = 2*arctan2(b/2,f)

		return {
			'alpha':rad2deg(alpha),\
			'beta':rad2deg(beta),
			'gamma':rad2deg(gamma)
			}

	###############################################################################
	#
	# find_lambda_set() 
	#
	#	Inputs:
	#		
	#	Outputs:
	#
	#	Notes:
	#
	###############################################################################

	def find_lambda_set(self, qe, tc, lambda_bin_size):
		from numpy import concatenate, arange
		all_lambdas = concatenate([qe['lambda'],tc['lambda']])
		min_lambda = min(all_lambdas)
		max_lambda = max(all_lambdas)
		lambda_set = arange(min_lambda,max_lambda + lambda_bin_size,lambda_bin_size)
		return lambda_set

	###############################################################################
	#
	# interpolate_lambda_dependent() 
	#
	#	Inputs:
	#		
	#	Outputs:
	#
	#	Notes: I moved the actual funciton to util within the python framework
	#		because it was useful elsewhere, too.
	# 		I let the function defined within the camera object to keep this
	#		architecture in place for when we port to C++.
	#
	###############################################################################

	def interpolate_lambda_dependent(self, ex,lambda_set):
		from util import interpolate_lambda_dependent
		return interpolate_lambda_dependent(ex,lambda_set)

	###############################################################################
	#
	# hot_dark() 
	#
	#	Inputs:
	#		
	#	Outputs:
	#
	#	Notes: 
	#
	###############################################################################

	# def calculate_hot_dark(self):
	# 	resolution_height = self.resolution_height
	# 	resolution_width = self.resolution_width
	# 	current_hot_dark = self.hot_dark
	# 	#do your thing here
		
	# 	self.hot_dark_arr = hot_dark_array

	###############################################################################
	#
	# update_state() 
	#
	#	Inputs:
	#		
	#	Outputs:
	#
	#	Notes: 
	#
	###############################################################################

	def update_state(self):

		# self.hot_dark = self.calculate_hot_dark(dt)

		#initialize a counter so we can check if we have open images
		open_images = 0
		take_image = self.msg['take_image']
		#loop through all images counting the ones that are open.
		#if you find an open one, set its key to open_image_key.
		#if there's more than one open, the var will be overwritten, but
		#that's OK since the code will error out later.
		for key in self.images:
			if self.images[key].image_open == 1: open_images +=1
			open_image_key = key

		if open_images > 1.:
			#we should never have more than one image open at a time, so
			#if open_images > 1, we have a problem. Report an error and
			#return
			print('ERROR: Camera Object Has Multiple Open Images!!!')
			return
		elif open_images == 1.:
			#if there is exactly one open image and we are still imaging
			#we need to update the state.
			if take_image == 1:
				self.images[open_image_key].update_state()
			#if there is exactly one open image but the exec has told the
			#camera to stop imaging, just close it. We still neet do run
			#the image's update_state() method in order to get it to actually
			#create the image. It only collects attitude information until
			#this step.
			else:
				self.images[open_image_key].image_open = 0
				self.images[open_image_key].update_state()
		else:
			if take_image == 1.:
				#if we have no open image, count the number of open images
				#and use that as the key for the next one (python indexes at
				#zero, so if we want to index using autoincrementing numbers, 
				#this will give us the next number).
				new_key = len(self.images.keys())
				#initialize the image.
				self.images[new_key] = image(self,self.msg)
				#give the image its first update.
				self.images[new_key].update_state()

		return

###############################################################################
#	Image is a class for the image data product produced by a spacecraft
#		camera.
#
#	methods:
#
#	User Variables:
#		alpha, beta, gamme: angles for (3-1-3) Euler angle rotation from
#			inertial from to camera frame. SEE NOTES!
#		camera: The camera object that the image is tied to.
#		msg: A dictionary used to control which parts of the code are run.
#			the main purpose is to control with bodies are removed by 
#			remove_occultations, but the user can also control how images
#			are made in order to demo and debug the code in various ways
#
#	Computed Variables:
# 		RA: Similar to the RA varable that is held in the camera object,
#			but with only stars/facets that are in the FOC of this image. 
# 		DE: Similar to the DE varable that is held in the camera object,
#			but with only stars/facets that are in the FOC of this image.
# 		n1: Similar to the n1 varable that is held in the camera object,
#			but with only stars/facets that are in the FOC of this image.
# 		n2: Similar to the n2 varable that is held in the camera object,
#			but with only stars/facets that are in the FOC of this image.
# 		n3: Similar to the n3 varable that is held in the camera object,
#			but with only stars/facets that are in the FOC of this image.
# 		VT: Similar to the VT varable that is held in the camera object,
#			but with only stars/facets that are in the FOC of this image.
# 		BVT: Similar to the BVT varable that is held in the camera object,
#			but with only stars/facets that are in the FOC of this image.
# 		c1: Similar to the c1 varable that is held in the camera object,
#			but with only stars/facets that are in the FOC of this image.
# 		c2: Similar to the c2 varable that is held in the camera object,
#			but with only stars/facets that are in the FOC of this image.
# 		c3: Similar to the c3 varable that is held in the camera object,
#			but with only stars/facets that are in the FOC of this image.
# 		pixel: Similar to the pixel varable that is held in the camera object,
#			but with only stars/facets that are in the FOC of this image.
# 		line: Similar to the line varable that is held in the camera object,
#			but with only stars/facets that are in the FOC of this image.
# 		color:
#		detector_array:
#
#	Keyword Arguments:
#
#	Notes:
#		-alpha, beta, gamma really need to be removed here. Although the
#		Euler angle transformation is a lot easier to grasp intuitively, it
#		would make way a camera-to-body DCM that is stored in the camera
#		object along with an inertial to spacecraft DCM that is stored in the
#		spacecraft object. This way the image object will only need camera and
#		msg passed in.
#		-This class ultimately needs to be split into image and frame classes
#		where frames are totaled in order to make the final image. The frames
#		should be calculated with a margin added to the true FOV because once
#		we start to add saturation effects, it will matter a lot what is right
#		outside the FOV (i.e. if the sun is right next to the FOV, the 
#		background counts will be very high, even though the sun is not imaged
#		directly)
#
###############################################################################

class image:
	def __init__(
		self, 
		camera,
		msg,
		**kwargs
		):

		self.image_open = 1
		self.alpha = []
		self.beta = []
		self.gamma = []
		self.camera = camera
		self.msg = 0
		self.RA = 0
		self.DE = 0
		self.n1 = 0
		self.n2 = 0
		self.n3 = 0
		self.VT = 0
		self.BVT = 0
		self.c1 = 0
		self.c2 = 0
		self.c3 = 0
		self.DCM = []
		self.scenes = []
		self.frames = []
		self.detector_array = []
		self.star_id = []

	###############################################################################
	#
	# update_state() 
	#
	#	Inputs:
	#		
	#	Outputs:
	#
	#	Notes: 
	#
	###############################################################################

	def update_state(self):

		if self.camera.msg['take_image'] == 1:
			#if the image is still actively being taken, all we do is save
			#off a DCM for that scene. The image will only be created once
			#self.camera.msg['take_image'] is set to zero and
			#image.update_state() is run.
			self.DCM.append(self.camera.sc.attitudeDCM)
		else:
			from scipy.linalg import inv
			from numpy import arctan2, arcsin, array, zeros, pi, floor
			from em import planck
			import pdb

			#this finds a subset of stars that are in the FOV at some point
			#during the exposure. Some functionality of find_stars_in_FOV()
			#is not needed in this step (like converting to pixel/line), 
			#but it is done anyway so we can reuse the function.

		
			#use the first attitude of the exposure as the central attitude
			#for calculating subset of stars.
			DCM_0 = self.DCM[0]
			DCM_0_inv = inv(DCM_0)
			alpha_0 = arctan2(DCM_0[0,1],DCM_0[0,0])
			beta_0 = arcsin(DCM_0[0,2])
			gamma_0 = arctan2(DCM_0[1,2],DCM_0[2,2])

			#All of the other DCMs at this point are transformations from
			#the intertial reference frame to the camera frame at the time
			#this scene was recorded. This loop converts them to be
			#transformations from camera frame at the beginning of the exposure
			#(i.e. DCM_0) to the camera frame at the time of this exposure.
			for i in range(0,len(self.DCM)):
				self.DCM[i] = self.DCM[i].dot(DCM_0_inv)
				self.alpha.append(arctan2(self.DCM[i][0,1],self.DCM[i][0,0]))
				self.beta.append(arcsin(self.DCM[i][0,2]))
				self.gamma.append(arctan2(self.DCM[i][1,2],self.DCM[i][2,2]))

			#this is getting sloppy... When we call find_stars_in_FOV()
			#to find all the stars in the exposure at any time, we want
			#to find all physical objects. This way we only render bodies
			#once per exposure. All the rest is just filtering them out
			#per scene and adding psf/noise characteristics of the image.
			full_exposure_msg = dict(self.camera.msg)
			full_exposure_msg['psf'] = 0
			full_exposure_msg['raster'] = 0
			full_exposure_msg['photon'] = 0
			full_exposure_msg['dark'] = 0
			full_exposure_msg['read'] = 0

			FOV = self.find_stars_in_FOV(
				alpha_0, #center Euler angle of First FOV of exposure
				beta_0,  #center Euler angle of First FOV of exposure
				gamma_0, #center Euler angle of First FOV of exposure
				max(self.alpha), #extreme euler angle of exposure
				max(self.beta),  #extreme euler angle of exposure
				min(self.alpha), #extreme euler angle of exposure
				min(self.beta),  #extreme euler angle of exposure
				#Don't know what twist (gamma) angle will be for
				#any given exposure, so assume worst case scenario
				#for bounding stars
				self.camera.angular_diagonal/2,
				self.camera.angular_diagonal/2,
				self.camera.resolution_height,
				self.camera.resolution_width,
				self.camera.RA, self.camera.DE, 
				self.camera.n1, self.camera.n2, self.camera.n3,
				self.camera.VT, self.camera.BVT, self.camera.T,
				self.camera.T, #spoof so the fcn won't break :-(
				self.camera.reduction_term,
				self.camera.max_mag, 
				self.camera.star_id,
				full_exposure_msg
				)

			self.RA = FOV['RA']
			self.DE = FOV['DE']
			self.n1 = FOV['n1']
			self.n2 = FOV['n2']
			self.n3 = FOV['n3']
			self.c1 = FOV['c1']
			self.c2 = FOV['c2']
			self.c3 = FOV['c3']
			self.T = FOV['T']
			self.star_id = FOV['star_id']
			self.reduction_term = FOV['reduction_term']
			I = []
			import matplotlib.pyplot as plt


			for i in range(0,len(self.T)): 

				T = self.T[i]
				if T == 5778: 
					flux_per_m2_per_nm_per_sr_at_star = self.camera.solar_bb
				else:
					flux_per_m2_per_nm_per_sr_at_star = planck(T,self.camera.lambda_set*1e-9)
				reduction_term = self.reduction_term[i]
				
				flux_per_m2_per_nm_per_sr_at_obs = flux_per_m2_per_nm_per_sr_at_star*reduction_term
				flux_per_m2_per_nm = pi*flux_per_m2_per_nm_per_sr_at_obs
				flux_per_m2 = flux_per_m2_per_nm*self.camera.lambda_bin_size
				flux = flux_per_m2*self.camera.effective_area
				photons_per_sec = flux/self.photon_energy(self.camera.lambda_set*1e-9)
				electrons_per_sec = photons_per_sec*self.camera.sensitivity_curve
				I.append(sum(electrons_per_sec))
			
			self.I = array(I)
			self.star_id = FOV['star_id']


			#since the camera object removes occulted objects
			#and adds body renderings, we don't need to do that
			#all again here, so we force those parts of the msg
			# to be zero.
			scene_msg = dict(self.camera.msg) 
			scene_msg['rm_occ'] = 0
			scene_msg['add_bod'] = 0
			#create one scene per DCM we collected above
			for i in range(0,len(self.DCM)):
				FOV = self.find_stars_in_FOV(
					self.alpha[i] + alpha_0, 
					self.beta[i] + beta_0, 
					self.gamma[i] + gamma_0, 
					0,
					0,
					0,
					0,
					self.camera.angular_height/2,
					self.camera.angular_width/2,
					self.camera.resolution_height,
					self.camera.resolution_width,
					self.RA, self.DE, 
					self.n1, self.n2, self.n3,
					self.VT, self.BVT, self.T,
					self.I,
					self.reduction_term,
					self.camera.max_mag, 
					self.star_id,
					scene_msg
					)
				self.scenes.append(
					scene(
						FOV['RA'],
						FOV['DE'],
						FOV['n1'],
						FOV['n2'],
						FOV['n3'],
						FOV['c1'],
						FOV['c2'],
						FOV['c3'],
						FOV['I'],
						FOV['pixel'],
						FOV['line'],
						FOV['star_id']
						)
					)
			i = 0
			for each_scene in self.scenes:
				print(i)
				i+=1
				if self.camera.msg['psf']:
					psf = self.psf(1)

					pixel = psf['x'].reshape(len(psf['x']),1) + each_scene.pixel
					line = psf['y'].reshape(len(psf['y']),1) + each_scene.line
					I = psf['I'].reshape(len(psf['I']),1)*each_scene.I

					each_scene.psf_pixel = pixel.reshape(len(each_scene.pixel)*len(psf['x']))
					each_scene.psf_line = line.reshape(len(each_scene.line)*len(psf['y']))

					each_scene.psf_I = I.reshape(len(psf['I'])*len(each_scene.I))


				if self.camera.msg['photon']:
					each_scene.psf_I = self.add_poisson_noise(each_scene.psf_I)

				if self.camera.msg['raster']:
					each_scene.detector_array = \
						self.rasterize(
							self.camera.resolution_width,
							self.camera.resolution_height,
							each_scene.psf_pixel,
							each_scene.psf_line,
							each_scene.psf_I
							)*self.camera.msg['dt']

					if self.camera.msg['dark']:
						each_scene.detector_array = \
							each_scene.detector_array + \
							self.add_poisson_noise(
								zeros(len(each_scene.detector_array)) + \
								self.camera.dark_current
								)

			self.detector_array = 0
			for each_scene in self.scenes:
				self.detector_array += each_scene.detector_array

			if self.camera.msg['read']:
				self.detector_array = self.add_read_noise(
						self.detector_array,
						self.camera.read_sigma
						)

			self.detector_array[self.detector_array < 0] = 0
			
			#convert phton counts to DN
			self.detector_array = floor(self.detector_array/self.camera.dnBinSize)

			#cut off pixels that are saturated
			self.detector_array[self.detector_array > self.camera.dnDepthMax] = self.camera.dnDepthMax
			# if self.camera.msg['hot_dark']:
			# 	self.detector_array = self.detector_array*self.camera.hot_dark
				# if 	self.camera.msg['raster'] + self.camera.msg['photon'] + \
				# 	self.camera.msg['read']:
				# 	self.detector_array = detector_array*hot_dark



	###########################################################################
	#
	# find_stars_in_FOV() finds the pixel and line coordinates of a 
	#
	#	Inputs:
	#		
	#	Outputs:
	#		A dict with the following entries:
	# 			pixel:
	# 			line:
	# 			RA:
	# 			DE:
	# 			n1:
	# 			n2:
	# 			n3:
	# 			VT:
	# 			BVT:
	# 			c1:
	# 			c2:
	# 			c3:
	#
	#	Notes:
	#
	###########################################################################
	def find_stars_in_FOV(
		self, 
		alpha, 
		beta, 
		gamma, 
		alpha_max,
		beta_max,
		alpha_min,
		beta_min,
		half_alpha,
		half_beta,
		alpha_resolution, 
		beta_resolution,
		RA, DE, 
		n1, n2, n3,
		VT, BVT,T,I,
		reduction_term,
		max_mag, 
		star_ids,
		msg,
		**kwargs):

		from numpy import array, deg2rad, sin, cos, append, sqrt, zeros, ones, logical_and
		from numpy.linalg import norm

		import pdb
		pdb.set_trace()
		if len(msg['bodies']) != 0:

			n = 1	
			bodies = msg['bodies']

			#calculate how far each body is from the sc
			for body in bodies: body.distFromSc = norm(body.state[0:3] - self.camera.sc.state[0:3])
			#sort bodies by how far they are from the sc
			#this needs to be done 
			bodies.sort(key=lambda x:x.distFromSc, reverse=True)

			for body in bodies:
				n+=1
				if body.name == 'SC': continue
				if msg['rm_occ']:
					occ_check = self.remove_occultations(body,n1,n2,n3)
					n1 = n1[occ_check]
					n2 = n2[occ_check]
					n3 = n3[occ_check]
					RA = RA[occ_check]
					DE = DE[occ_check]
					star_ids = star_ids[occ_check]
					T = T[occ_check]
					reduction_term = reduction_term[occ_check]
					I = I[occ_check]


				if msg['add_bod']:
					print(body.name)
					from lightSimFunctions import lightSim
					DCM = self.camera.body2cameraDCM.dot(self.camera.sc.attitudeDCM)
					facets = lightSim(
						DCM,
						self.camera.sc.state[0:3],
						body.state[0:3],
						(
							self.camera.angular_height,
							self.camera.angular_width
							),
						200,
						200,
						False,
						body.albedo,
						body.r_eq,
						body.name
						)
					if facets == -1: continue #if true, then lightSim FOV check failed
					#position from center of body to facet

					surf_n1 = facets['bodypos'][:,0].astype(float)
					surf_n2 = facets['bodypos'][:,1].astype(float)
					surf_n3 = facets['bodypos'][:,2].astype(float)

					surf_n1_tmp = surf_n1	
					surf_n2_tmp = surf_n2	
					surf_n3_tmp = surf_n3	
					#position from [0,0,0] HCI to facet
					surf_n1 += body.state[0]
					surf_n2 += body.state[1]
					surf_n3 += body.state[2]

					#camera to facet vectors
					surf_n1 -= self.camera.sc.state[0]
					surf_n2 -= self.camera.sc.state[1]
					surf_n3 -= self.camera.sc.state[2]

					surf_r = surf_n1**2 + surf_n2**2 + surf_n3**2
					#turn to unit vectors
					surf_n1 = surf_n1/sqrt(surf_r)
					surf_n2 = surf_n2/sqrt(surf_r)
					surf_n3 = surf_n3/sqrt(surf_r)

					n1 = append(n1, surf_n1)
					n2 = append(n2, surf_n2)
					n3 = append(n3, surf_n3)
					RA = append(RA, zeros(len(surf_n1)))
					DE = append(DE, zeros(len(surf_n1)))
					VT = append(VT, zeros(len(surf_n1)))
					BVT = append(BVT, zeros(len(surf_n1)))
					T = append(T, zeros(len(surf_n1)) + 5778)
					star_ids = append(star_ids, zeros(len(surf_n1)))

					if len(surf_n1) > 1:
						print('>1')
						reduction_term = append(
							reduction_term,
							facets['netAlbedo']*\
							facets['facetArea'])
						I = append(I,
							facets['netAlbedo']*\
							facets['facetArea'])
					else:
						print('1')
						reduction_term = append(
							reduction_term,
							sum(facets['netAlbedo']*\
							facets['facetArea']))
						I = append(I,
							sum(facets['netAlbedo']*\
							facets['facetArea']))
						print('2')

		c2_max = alpha_max + sin(deg2rad(half_alpha))
		c2_min = alpha_min - sin(deg2rad(half_alpha))
		c3_max = beta_max + sin(deg2rad(half_beta))
		c3_min = beta_min - sin(deg2rad(half_beta))

		#rotate all stars into camera frame coordinates
		#this can be done far more efficiently with a DCM
		#but I'm not really ready to test it yet
		c1 = \
		 n1*cos(beta)*cos(alpha) \
		+n2*cos(beta)*sin(alpha) \
		-n3*sin(beta)

		c2 = \
		 n1*(
		 	sin(gamma)*sin(beta)*cos(alpha) - \
		 	cos(gamma)*sin(alpha)
		 	) \
		+n2*(
			sin(gamma)*sin(beta)*sin(alpha) + \
			cos(gamma)*cos(alpha)
			) \
		+n3*sin(gamma)*cos(beta)

		c3 = \
		 n1*(
		 	cos(gamma)*sin(beta)*cos(alpha) + \
		 	sin(gamma)*sin(alpha)
		 	) \
		+n2*(
			cos(gamma)*sin(beta)*sin(alpha) - \
			sin(gamma)*cos(alpha)
			) \
		+n3*cos(gamma)*cos(beta)

		#Remove stars outside the FOV in the c2 direction
		ind = abs(c2/c1*self.camera.focal_length) < self.camera.detector_width/2
		#Remove stars outside the FOV in the c3 direction
		ind = logical_and(ind,abs(c3/c1*self.camera.focal_length) < self.camera.detector_height/2)
		#remove stars in the anti-boresight direction
		ind = logical_and(ind,c1 > 0)

		#apply logic from above
		RA = RA[ind]
		DE = DE[ind]
		n1 = n1[ind]
		n2 = n2[ind]
		n3 = n3[ind]
		c1 = c1[ind]
		c2 = c2[ind]
		c3 = c3[ind]
		star_ids = star_ids[ind]
		T = T[ind]
		reduction_term = reduction_term[ind]
		I = I[ind]

		#using similar triangles
		pixel = -self.camera.focal_length*c2/c1*\
			beta_resolution/self.camera.detector_width + \
			beta_resolution/2

		line = -self.camera.focal_length*c3/c1*\
			alpha_resolution/self.camera.detector_height + \
			alpha_resolution/2

		return {
			'pixel': pixel,
			'line': line,
			'RA': RA,
			'DE': DE,
			'n1': n1,
			'n2': n2,
			'n3': n3,
			'c1': c1,
			'c2': c2,
			'c3': c3,
			'star_id' : star_ids,
			'I': I,
			'T': T,
			'reduction_term': reduction_term
		}
	###########################################################################
	#
	# remove_occultations() 
	#
	#	Inputs:
	#		body: Object. Instantiation of the body class from bodies.py.
	#			Variables used from that class include equatorial radius,
	#			polar radius, and state (for the position to the body).
	#
	#		n1, n2, n3: Numpy Float Arrays. All should be the same size. Each
	#			star/body fact that is still in the FOV at this point should
	#			have an entry in each, and each index MUST correspond to data
	#			from the same star/facet!
	#
	#
	#	Outputs:
	#		occ_check: Numpy Array. A list of booleans, each corresponding to 
	#
	###########################################################################
	def remove_occultations(self,body,n1,n2,n3):
		from numpy import array, stack, einsum, logical_or
		from numpy.linalg import norm
		#this needs to be multiplied by a transformation matrix in order
		#to account for planetary attitude, but I haven't thought through
		#that just yet.
		A = array([
			[body.r_eq**-2,0,0],
			[0,body.r_eq**-2,0],
			[0,0,body.r_pole**-2]
			])

		v = stack([n1,n2,n3])
		#agnostic to coordinate frame origin. Just needs same axis
		#directions as v.
		y_0 = self.camera.sc.state[0:3] - body.state[0:3]
		#1xn
		vTAy_0 = (v.T.dot(A).dot(y_0)).T
		#1xn
		vTAv = einsum('ij,ji->i',v.T,A.dot(v))
		#scalar
		y_0TAy_0 = y_0.T.dot(A).dot(y_0)


		discriminant = vTAy_0**2 - vTAv*(y_0TAy_0 - 1)
		occ_check = logical_or(discriminant < 0,v.T.dot(y_0) > 0)
		return occ_check

	###########################################################################
	#
	# diego() 
	#
	#	Inputs:
	#		r, r_0, p, p_r: Floats. 
	#
	#	Outputs:
	#		I
	#
	###########################################################################
	def diego(self,r,r_0,p,p_r):
		I = (1+(r/r_0)**(p*(1+(r/p_r))))**-1
		return I

	###########################################################################
	#
	# gaussian() 
	#
	#	Inputs:
	#
	#	Outputs:
	#
	###########################################################################
	def gaussian(self,r,sigma):
		from numpy import exp
		return exp(-r**2/(2*sigma**2))

	###########################################################################
	#
	# psf() 
	#
	#	Inputs:
	#
	#	Outputs:
	#
	###########################################################################
	def psf(self,sigma):
		from numpy import arange, pi, ones
		from numpy import array, sin, cos, append, linspace
		r_3sig = sigma*3
		rad = arange(0.125,r_3sig+0.25,0.25)

		x = array([])
		y = array([])
		I = array([])

		for r in rad:
			c = 2*pi*r
			theta = linspace(0,2*pi,c*4)[0:-1]
			x = append(x,r*cos(theta))
			y = append(y,r*sin(theta))
			I = append(I,ones(len(theta))*self.gaussian(r,1))
		I = I/sum(I)

		return { "x":x, "y":y, "I":I }

	###########################################################################
	#
	# rasterize() 
	#
	#	Inputs:
	#
	#	Outputs:
	#
	###########################################################################
	def rasterize(self,pixel_resolution,line_resolution,pixel_coord,line_coord,intensity):
		from numpy import floor, zeros, array, arange, append 
		from numpy import concatenate, logical_and
		from pandas import DataFrame

		#adding PSF introduces some values that are not on the detector. Remove them here

		positive_coords = logical_and(pixel_coord > 0, line_coord > 0)
		pixel_coord = pixel_coord[positive_coords]
		line_coord = line_coord[positive_coords]
		intensity = intensity[positive_coords]

		not_too_big = logical_and(pixel_coord < pixel_resolution, line_coord < line_resolution)
		pixel_coord = pixel_coord[not_too_big]
		line_coord = line_coord[not_too_big]
		intensity = intensity[not_too_big]

		int_pix_coord = floor(pixel_coord).astype(int)
		int_line_coord = floor(line_coord).astype(int)

		detector_position = (line_resolution*int_line_coord + int_pix_coord)
		detector_position = append(detector_position,arange(pixel_resolution*line_resolution))
		intensity = append(intensity,zeros(pixel_resolution*line_resolution))

		data = concatenate([detector_position,intensity])
		data = data.reshape(2,int(len(data)/2)).T
		df = DataFrame(data,columns = ["Position","Intensity"])
		detector_array = df.groupby("Position").sum().values.T[0]

		return detector_array

	###########################################################################
	#
	# add_read_noise() 
	#
	#	Inputs:
	#
	#	Outputs:
	#
	###########################################################################
	def add_read_noise(self,detector_array,sigma):
		from numpy.random import randn
		return detector_array + sigma*randn(len(detector_array))

	###########################################################################
	#
	# add_poisson_noise() 
	#
	#	Inputs:
	#
	#	Outputs:
	#
	###########################################################################
	def add_poisson_noise(self, I):
		from numpy.random import poisson
		return I + poisson(I)

	###########################################################################
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
	###########################################################################

	def planck(self, T,lam):
		from constants import h, c, k_B
		from numpy import exp

		top_part = 2*h*c**2
		bottom_part = lam**5*(exp(h*c/(lam*k_B*T))-1)

		I = top_part/bottom_part*1e-9 #1e-9 to convert to per nm

		return I

	###########################################################################
	#	stefan_boltzmann() is a function that total flux from a star given 
	#
	#	Inputs:
	#		T: temperature of the star in Kelvin
	#
	#	Outputs:
	#		F: total flux at stellar surface in W/m^2
	#
	#	Notes:
	#		This function isn't used in DINO C-REx except as part of the test
	#		for planck(). The Planck function integrated over all wavelength
	#		space should be identically equal to the Stephan-Boltzmann function
	#
	###########################################################################

	def stefan_boltzmann(self, T):
		from constants import sigma
		return sigma*T**4

	###########################################################################
	#	photon_energy() 
	#
	#	Inputs:
	#
	#	Outputs:
	#
	#	Notes:
	#
	###########################################################################

	def photon_energy(self, lam):
		from em import photon_energy
		return photon_energy(lam)
		from constants import h, c
		return h*c/lam

	###########################################################################
	#	addZodiacalLight() 
	#
	#	Inputs: 
	#		attitudeDCM: DCM for the inertial to camera attitude.
	#		
	#	Outputs:
	#		zodiacalFlux: 
	#
	#	Notes:
	#
	# ###########################################################################
	def addZodiacalLight(self):
		from numpy import array, identity, flip, vstack, hstack, transpose
		from numpy import arctan2, pi, rad2deg, arcsin, argmin
		sunAngleRad = arctan2(self.camera.sc.attitudeDCM[1],self.camera.sc.attitudeDCM[0]) + pi
		sunAngleDeg = rad2deg(sunAngleRad)

		raBins = array([
			0,5,10,15,20,25,30,35,40,45,60,75,90,105,120,135,150,165,180
			])
		decBins = array([
			0,5,10,15,20,25,30,45,60,75
			])
		C12 = self.camera.sc.attitudeDCM[0,1]
		C11 = self.camera.sc.attitudeDCM[0,0]
		C13 = self.camera.sc.attitudeDCM[0,2]

		camRArad = arctan2(C12,C11)
		camDErad = -arcsin(C13)

		camRAdeg = abs(rad2deg(camRArad))
		camDEdeg = abs(rad2deg(camDErad))
		camDEbin = argmin(abs(camDEdeg-decBins))
		camRAbin = argmin(abs(camRAdeg-raBins))

		zodiacalLight = array([
			[0, 0, 0, 2450, 1260, 770, 500, 215, 117, 78],
			[0, 0, 0, 2300, 1200, 740, 490, 212, 117, 78],
			[0, 0, 3700, 1930, 1070, 675, 460, 206, 116, 78],
			[9000, 5300, 2690, 1450, 870, 590, 410, 196, 114, 78],
			[5000, 3500, 1880, 1100, 710, 495, 355, 185, 110, 77],
			[3000, 2210, 1350, 860, 585, 425, 320, 174, 106, 76],
			[1940, 1460, 955, 660, 480, 365, 285, 162, 102, 74],
			[1290, 990, 710, 530, 400, 310, 250, 151, 98, 73],
			[925, 735, 545, 415, 325, 264, 220, 140, 94, 72],
			[710, 570, 435, 345, 278, 228, 195, 130, 91, 70],
			[395, 345, 275, 228, 190, 163, 143, 105, 81, 67],
			[264, 248, 210, 177, 153, 134, 118, 91, 73, 64],
			[202, 196, 176, 151, 130, 115, 103, 81, 67, 62],
			[166, 164, 154, 133, 117, 104, 93, 75, 64, 60],
			[147, 145, 138, 120, 108, 98, 88, 70, 60, 58],
			[140, 139, 130, 115, 105, 95, 86, 70, 60, 57],
			[140, 139, 129, 116, 107, 99, 91, 75, 62, 56],
			[153, 150, 140, 129, 118, 110, 102, 81, 64, 56],
			[180, 166, 152, 139, 127, 116, 105, 82, 65, 56]
		])

		#zodiacalLight = zodiacalLight*1.28e-8*.156

		return zodiacalLight[camRAbin,camDEbin]

	###########################################################################
	#	addBkgdStarLight() 
	#
	#	Inputs: 
	#		attitudeDCM: DCM for the inertial to camera attitude.
	#		
	#	Outputs:
	#		bkgdFlux: 
	#
	#	Notes: the data represented in pioneer10BkdgStarLight comes from
	#	Lienert's 1997 reference on diffuse night sky brightness on p78 and 79.
	# 
	#	Units for the array are S10_sun, given by leinart to be 1.28e-8*.156
	#	W/m^2/sr
	#
	#	Rows of the array represent bins of declination, and are spaced by 10
	#	degrees. Columns are binned by right ascension and are spaced every 10
	#	degrees.
	#
	# ###########################################################################
	def addBkgdStarLight(self):

		from numpy import array, arcsin, arctan2, rad2deg, deg2rad, flip

		decBins = array([-80, -70, -60, -50, -40, -30, -24, -20, -16, -12, -8, -4, 0,4, 8, 12, 16, 20, 24, 30, 40, 50, 60, 70, 80])
		raBins = array([0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350])

		C12 = self.camera.sc.attitudeDCM[0,1]
		C11 = self.camera.sc.attitudeDCM[0,0]
		C13 = self.camera.sc.attitudeDCM[0,2]

		camRArad = arctan2(C12,C11)
		camDErad = -arcsin(C13)

		camRAdeg = rad2deg(camRArad)
		camDEdeg = rad2deg(camDErad)
		if camDEdeg < 0: camDEdeg
		if camRAdeg < 0: camRAdeg += 360


		#probably a better way to do this...
		#the place where the difference between camera DE and the decliation bins
		#given by Leinert is the smallest will be the delination bin that we want.
		camDEbin = min(abs(camDEdeg-decBins)) == abs(camDEdeg-decBins)
		camRAbin = min(abs(camRAdeg-raBins)) == abs(camRAdeg-raBins)



		pioneer10BkdgStarLight = array([
			[ 68,  55,  55,  56,  49,   35,  39,  40,  38,  42,  49,  46,  50,  56,  43,  51,  59,  56,  61,  71, 116, 192, 284, 162, 124],                  
			[ 74,  74,  42,  38,  36,   35,  34,  43,  39,  34,  44,  45,  35,  48,  47,  51,  55,  52,  63,  64, 108, 181, 261, 185, 121],
			[ 78,  67,  51,  44,  40,   31,  36,  40,  39,  39,  43,  38,  32,  53,  45,  40,  56,  58,  57,  79,  99, 158, 250, 176, 106],
			[ 78,  53,  48,  44,  39,   38,  34,  46,  45,  41,  40,  51,  43,  47,  64,  56,  61,  56,  66,  82,  97, 187, 265, 157, 110],
			[ 73,  69,  49,  42,  36,   36,  45,  42,  40,  44,  45,  48,  53,  49,  59,  54,  61,  62,  60,  74, 116, 177, 240, 150, 109],
			[ 72,  63,  56,  41,  48,   50,  44,  41,  44,  51,  49,  56,  60,  59,  61,  70,  74,  79,  78,  85, 141, 168, 205, 140, 111],
			[ 77,  71,  53,  52,  48,   50,  51,  52,  52,  57,  60,  68,  68,  73,  88,  85,  82,  90, 103, 108, 133, 183, 175, 128, 117],
			[ 81,  90,  72,  62,  54,   58,  66,  64,  64,  74,  76,  75,  87,  82,  93,  95, 112, 122, 102, 104, 189, 164, 161, 137, 103],
			[ 81, 341,  64,  72,  64,   63,  80,  84,  90, 115, 134, 131, 136, 135, 135, 125, 146, 143, 130, 179, 214, 284, 125, 115,  85],
			[ 88, 226,  88,  93,  83,   90, 117, 123, 138, 147, 155, 158, 174, 174, 220, 238, 258, 218, 246, 227, 169, 137, 112,  91,  93],
			[ 98, 107, 114, 116, 135,  176, 193, 248, 240, 248, 225, 253, 283, 280, 284, 263, 229, 233, 207, 210, 145, 119, 108,  77,  84],
			[102, 111, 141, 182, 218,  315, 331, 321, 297, 312, 295, 267, 225, 206, 182, 184, 165, 135, 138, 105, 108,  93,  81,  73,  77],
			[110, 141, 176, 275, 308,  385, 327, 288, 249, 200, 166, 178, 173, 167, 155,  -1,  -1,  -1,  -1, 105,  75,  63,  75,  79,  73],
			[113, 142, 224, 337, 272,  274, 201, 156, 137,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  77,  52,  57,  74,  70],
			[119, 206, 372, 294, 398,  175, 114, 116, 126, 109, 101,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  53,  57,  80,  65],
			[116, 224, 490, 334, 193,  120, 104,  85,  72,  68,  65,  65,  66,  67,  68,  -1,  -1,  -1,  -1,  -1,  -1,  52,  55,  57,  55],
			[116, 265, 846, 250, 156,  100,  75,  74,  62,  65,  57,  68,  71,  71,  65,  62,  59,  60,  58,  55,  52,  46,  52,  46,  49],
			[141, 308, 744, 229, 116,   88,  76,  59,  63,  61,  56,  52,  55,  60,  68,  66,  58,  52,  52,  50,  60,  36,  51,  54,  52],
			[145, 340, 632, 230,  99,   78,  61,  61,  57,  58,  54,  56,  54,  54,  41,  47,  50,  31,  45,  38,  37,  43,  45,  38,  52],
			[137, 335, 450, 214, 120,   72,  67,  59,  63,  54,  51,  42,  43,  55,  45,  41,  41,  36,  41,  33,  37,  35,  42,  50,  51],
			[145, 410, 565, 262, 115,   91,  77,  59,  53,  58,  51,  49,  44,  48,  40,  34,  39,  42,  38,  38,  34,  37,  45,  50,  58],
			[154, 351, 516, 256, 144,   95,  73,  68,  66,  71,  58,  52,  47,  50,  49,  44,  57,  44,  40,  40,  40,  43,  43,  48,  59],
			[153, 337, 559, 284, 173,  102,  87,  81,  80,  77,  54,  53,  58,  48,  53,  47,  42,  47,  54,  38,  35,  44,  40,  54,  61],
			[156, 310, 538, 481, 249,  134, 110,  88,  97,  79,  77,  62,  71,  60,  54,  50,  52,  43,  46,  41,  55,  44,  49,  61,  69],
			[146, 262, 537, 357, 242,  200, 152, 137, 106, 101,  82,  83,  77,  77,  64,  64,  54,  54,  55,  53,  47,  51,  53,  60,  68],
			[128, 213, 373, 571, 485,  366, 217, 145, 141, 145, 108, 110, 101,  89,  87,  79,  75,  68,  67,  59,  57,  59,  56,  65,  71],
			[127, 190, 284, 494, 550,  571, 344, 356, 301, 172, 144, 119, 136, 120, 123, 119, 116,  90,  85,  89,  81,  66,  64,  78,  75],
			[115, 155, 205, 284, 626, 1393, 649, 498, 474, 261, 199, 162, 204, 225, 201, 166, 157, 148, 132, 120,  94,  96,  73,  82,  76],
			[101, 106, 150, 135, 262,  461, 481, 534, 542, 515, 509, 343, 210, 254, 382, 341, 237, 218, 206, 176, 145, 111,  96,  91,  81],
			[103, 113, 129, 191, 141,  183, 211, 246, 272, 274, 260, 257, 331, 366, 374, 271, 283, 250, 327, 387, 228, 168, 117,  90,  89],
			[102, 107,  93,  99, 112,  112, 117, 134, 144, 143, 162, 149, 178, 215, 237, 323, 366, 407, 316, 421, 403, 245, 163, 102,  93],
			[ 93,  88,  78,  80,  77,   81,  89,  85,  89,  91,  91, 109, 128, 138, 145, 174, 168, 227, 240, 288, 280, 292, 187, 133,  94],
			[ 87,  81,  72,  65,  64,   73,  70,  73,  72,  77, 100,  75,  85,  93, 109, 107, 131, 119, 153, 192, 279, 286, 225, 136, 105],
			[ 75,  76,  56,  54,  60,   42,  66,  60,  56,  67,  67,  63,  81,  78,  78, 103,  85,  93, 105, 142, 188, 384, 259, 167, 119],
			[ 66,  69,  57,  52,  51,   53,  43,  48,  58,  64,  57,  47,  49,  67,  65,  58,  63,  76,  76,  97, 171, 268, 241, 178, 122],
			[ 66,  64,  59,  43,  46,   56,  48,  49,  43,  50,  44,  47,  62,  48,  55,  59,  61,  69,  72,  84, 104, 203, 256, 166, 110] 
			])



		return pioneer10BkdgStarLight[camRAbin,camDEbin][0]

###############################################################################
#	scene is a class for collecting each integration step during an exposure.
#	These will then be summed to make frames, which are the physically relevant
#	unit of an image. frames are then summed to create the full images.
#
#
###############################################################################

class scene:
	def __init__(
		self,
		RA,
		DE,
		n1, 
		n2,
		n3,
		c1, 
		c2,
		c3,
		I,
		pixel,
		line,
		star_ids,
		**kwargs
		):

		self.RA = RA
		self.DE = DE
		self.n1 = n1
		self.n2 = n2
		self.n3 = n3
		self.c1 = c1
		self.c2 = c2
		self.c3 = c3
		self.I = I
		self.pixel = pixel
		self.line = line
		self.star_ids = star_ids






