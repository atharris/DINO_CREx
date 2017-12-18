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
#	Camera is a class for the simulation of a spacecraft camera. It holds all
#	physical characteristics of the camera and all stars. The Camera object
#	controls when image objects are opened and closed based on the takeImage
#	message read in from the navigation executive.
#
#	Methods:
#		loadAllStars()
#		calculateFOV()
#		findLambdaSet()
#		interpolateLambda()
#		calculateHotDark()
#		updateState()
#
#	User Variables:
#		detectorHeight: Physical height of detector
#			Can be any dimension so long as it matches detectorWidth and
#			focalLength
#		detectorWidth: Physical width of detector
#			Can be any dimension so long as it matches detectorHeight and
#			focalLength
#		focalLength: Physical distance between lens and focal point
#			Can be any dimension so long as it matches detectorWidth and
#			detectorWidth
#		resolutionHeight: # of pixels in same dimension as detectorHeight
#		resolutionWidth: # of pixels in same dimension as detectorWidth
#		body2cameraDCM: DCM describing orientation of camera relative to sc
#		maxMag: Maximum magnitude star visible to this camera.
#		maxMag: Minimum magnitude star visible to this camera (not physically
#			realistic, but useful for debugging).
#		qe: Quantum Effeciency Curve as a dict of two numpy arrays. 
#			One must be called 'lam', and the other must be called 
#			must be called 'throughput'. qe['lam'] must be an array of 
#			wavelength bins, and qe['throughput'] must be an array of 
#			percentages detected at the given 'lam'. qe['lam'] must
#			be in nm, and qe['throughput'] should be a number between 0 and
#			1. qe['lam'] and tc['lam'] do not necessarily need to
#			match as interpolateLambda will interpolate qe['throughput']
#			so that they do.
#		tc: Transmission Curve as a dict of two numpy arrays. 
#			One must be called 'lam', and the other must be called 
#			must be called 'throughput'. tc['lam'] must be an array of 
#			wavelength bins, and tc['throughput'] must be an array of 
#			percentages transmitted at the given 'lam'. tc['lam'] must
#			be in nm, and tc['throughput'] should be a number between 0 and
#			1. qe['lam'] and tc['lam'] do not necessarily need to
#			match as interpolateLambda will interpolate tc['throughput']
#			so that they do.
#		lambdaBinSize: the width of bins in wavelength space to be
#			returned from interpolateLambda().
#		effectiveArea: effective area of camera in m^2
#		darkCurrent: signal due to dark current in electron/s. Added to
#			every pixel in image, plus poisson photon noise
#		readSigma: standard deviation for gaussian read noise
#		dnBinSize: number of photons per intensity bin
#		dnDepthMax: saturation level of detector
#		sc: Spacecraft object this camera is attached to
#		msg: debug message used to turn features of camera on or off
#		hotDark: NOT YET IMPLEMENTED!
#
#	Computed Variables:
#		angularHeight: Size of detector field of view on the sky.
#		angularWidth: 
#		RA: Right ascension of all stars loaded from db
#		DE: Declination of all stars loaded from db
#		n1: First unit vector component in inertial space of all stars loaded 
#			from db via loadAllStars()
#		n2: Second unit vector component in inertial space of all stars loaded 
#			from db via loadAllStars()
#		n3: Third unit vector component in inertial space of all stars loaded 
#			from db via loadAllStars()
#		VT: Tycho Visual band magnitude of all stars loaded from db
#		BVT: Tycho (B-V) color index of all stars loaded from db
#		starID: DINO C-REx catalog number loaded from db
#		T: Estimated stellar temperature loaded from db
#		solidAngleSubtended: estimated solid angle subtended by each star.
#			loaded from db.
#		angularHeight: Angular height of the detector field of view. 
#			Calculated using calculateFOV()
#		angularWidth: Angular width of the detector field of view. 
#			Calculated using calculateFOV()
#		angularDiagonal: Angular size of maximum dimension of detector
#			(the diagonal since it's a rectangular detector)
#		solarBB: Blackbody intensity of the sun measured at the wavelengths
#			in lambdaSet. Used by image.updateState() to calculate beacon
#			intenisty incident on detector
#		lambdaSet: unified set of wavelength values for blackbody,
#			transmission, and quantum efficiency curves
#		qe: qe curve interpolated at wavelengths of lambdaSet
#		tc: transmission curve interpolated at wavelengths of lambdaSet
#		sensitivityCurve: Full throughput curve (qe*tc)
#		images: dictionary to hold images
#
#	Keyword Arguments:
#		db: File path to a stellar database file. If not provided, 
#			db/tycho_small.db is used.
# 
###############################################################################
class camera:
	def __init__(
		self, 
		detectorHeight, 
		detectorWidth, 
		focalLength,
		resolutionHeight,
		resolutionWidth,
		body2cameraDCM,
		maxMag,
		minMag,
		qe,
		tc,
		lambdaBinSize,
		effectiveArea,
		darkCurrent,
		readSigma,
		dnBinSize,
		dnDepthMax,
		psfSigma,
		dt,
		scState,
		scDCM,
		bodies,
		takeImage,
		**kwargs
		):

		from em import planck
		from numpy import pi, array

		try: 
			db = kwargs['db']
		except:
			db = 'db/tycho.db'
		
		try:
			msg = kwargs['debug']
		except:
			import bodies as bod
			msg = {
			'addStars': 1,'rmOcc': 1, 'addBod': 1, 'psf': 1, 
			'raster': 1, 'photon': 1, 'dark': 1, 'read': 1, 
			'verbose': 0}


		FOV = self.calculateFOV(
			focalLength,
			detectorHeight, 
			detectorWidth,
			)

		if msg['addStars']:
			allstars = self.loadAllStars(db, maxMag, minMag)
		lambdaSet = self.findLambdaSet(qe, tc, lambdaBinSize)
		qe = self.interpolateLambdaDependent(qe,lambdaSet)
		tc = self.interpolateLambdaDependent(tc,lambdaSet) 
		sensitivityCurve = qe['throughput']*tc['throughput']
		lambdaSet = lambdaSet[sensitivityCurve != 0]

		if msg['addStars']:
			self.RA = allstars['RA']
			self.DE = allstars['DE']
			self.n1 = allstars['n1']
			self.n2 = allstars['n2']
			self.n3 = allstars['n3']
			self.VT = allstars['VT']
			self.BVT = allstars['BVT']
			self.starID = allstars['starID']
			self.T = allstars['T']
			self.solidAngleSubtended = allstars['solidAngleSubtended']
		else:
			## Numpy array of right ascensions of all stars loaded into camera
			self.RA = array([]) 
			## Numpy array of declinations of all stars loaded into camera
			self.DE = array([]) 
			## Numpy array of 1st inertial unit vector coordinate of all stars loaded into camera
			self.n1 = array([]) 
			## Numpy array of 2nd inertial unit vector coordinate of all stars loaded into camera
			self.n2 = array([]) 
			## Numpy array of 3rd inertial unit vector coordinate of all stars loaded into camera
			self.n3 = array([]) 
			## Numpy array of Tycho visual magnitude of all stars loaded into camera
			self.VT = array([])
			## Numpy array of Tycho color index of all stars loaded into camera
			self.BVT = array([])
			## Numpy array of DINO db star id of all stars loaded into camera
			self.starID = array([])
			## Numpy array of computed temperatures of all stars loaded into camera
			self.T = array([])
			## Numpy array of solid angle subtended of all stars loaded into camera
			self.solidAngleSubtended = array([])

		## Numpy array solar flux values (W/m^2/sr) evaluated at each
		## wavelength in lambdaSet
		self.solarBB = planck(5778,lambdaSet*1e-9)
		self.lambdaSet = lambdaSet
		self.sensitivityCurve = sensitivityCurve[sensitivityCurve != 0]
		self.qe = qe
		self.tc = tc
		self.lambdaBinSize = lambdaBinSize
		self.resolutionHeight = resolutionHeight
		self.resolutionWidth = resolutionWidth		
		self.body2cameraDCM = body2cameraDCM
		self.maxMag = maxMag
		self.minMag = minMag
		self.angularHeight = FOV['alpha']
		self.angularWidth = FOV['beta']
		self.angularDiagonal = FOV['gamma']
		self.detectorHeight = detectorHeight
		self.detectorWidth = detectorWidth
		self.focalLength = focalLength
		self.readSigma = readSigma
		self.effectiveArea = effectiveArea
		self.darkCurrent = darkCurrent
		self.images = {}
		self.msg = msg
		self.dnBinSize = dnBinSize
		self.dnDepthMax = dnDepthMax
		self.psfSigma = psfSigma
		self.dt = dt
		self.bodies = bodies
		self.scState = scState
		self.scDCM = scDCM
		self.takeImage = takeImage
		self.imgTime = []

	###########################################################################
	#	loadAllStars() is used to load stellar data from a database (nominally
	# 	db/tycho.db into a camera instance).
	###########################################################################
	def loadAllStars(self,db, maxMag, minMag):
		"""!
		@param db: database to load stars from. If camera object is initialized 
		nominally, this will be db/tycho.db. db can be changed by initializing
		camera with the db kwarg.
		@param maxMag: maximum magnitude to load into camera. Can be used when
		modeling a noisy camera. Images are created more quickly when fewer stars
		are present, so cutting out dimmer stars may be desirable
		@param minMag: minum magnitude to load into camera. WARNING: Setting this
		value to anything greater than -2 is physically unrealistic as it will remove
		all of the brightest stars! Parameter is included for debugging purposes
		@return RA: Right ascension of all stars in db within magnitude bounds
		@return DE: Declination of all stars in db within magnitude bounds
		@return VT: Tycho visual magnitude of all stars in db within magnitude bounds
		@return BVT: Tycho color index of all stars in db within magnitude bounds
		@return n1: 1st coordinate of intertial unit vector of all stars in db within magnitude bounds
		@return n2: 2nd coordinate of intertial unit vector of all stars in db within magnitude bounds
		@return n3: 3rd coordinate of intertial unit vector of all stars in db within magnitude bounds
		@return starID: Right ascension of all stars in db within magnitude bounds
		@return solidAngleSubtended: Right ascension of all stars in db within magnitude bounds
		@return T: Right ascension of all stars in db within magnitude bounds
		"""
		import sqlite3
		from numpy import array, sin, cos, deg2rad, logical_and

		conn = sqlite3.connect(db)
		c = conn.cursor()

		selectString = "SELECT RA, DE, VTmag, BTmag, id, reduction_term, " + \
			"computed_temperature from tycho_data"
		c.execute(selectString)


		RA = []
		DE = []
		VT = []
		BT = []
		phi = []
		theta = []
		starID = []
		solidAngleSubtended = []
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
			starID.append(float(row[4]))
			#this IS a cheat. force the flux from the star
			#to be zero if we haven't gotten a computed T in for
			#it yet. The temp doesn't matter if the reduction
			#term is zero
			try:
				solidAngleSubtended.append(float(row[5]))
			except:
				solidAngleSubtended.append(float(0))
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
		starID = array(starID)
		solidAngleSubtended = array(solidAngleSubtended)
		T = array(T)

		ind = VT < maxMag
		ind = logical_and(ind, VT > minMag)
		RA = RA[ind]
		DE = DE[ind]
		VT = VT[ind]
		BVT = BVT[ind]
		theta = theta[ind]
		phi = phi[ind]
		starID = starID[ind]
		T = T[ind]
		solidAngleSubtended = solidAngleSubtended[ind]
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
			'starID': starID,
			'solidAngleSubtended': solidAngleSubtended,
			'T': T	
			}

	###############################################################################
	#
	# calculateFOV() finds the angular size of a camera FOV using freshman-level
	#	optics. The method uses a thin lens approximation.
	#
	#	Notes:
	#		Assumes lens is thin. Units for detector dimensions and focal length
	#		can be anything, but they must all match.
	#
	###############################################################################
	def calculateFOV(self, focalLength, detectorHeight, detectorWidth, **kwargs):
		"""!
		@param focalLength: Focal length of modeled lens 
		@param detectorHeight: y-direction dimension of the lens
		@param detectorWidth: x-direction dimension of the lens
		@return A dict with the following entries:
		@return alhpa: angular size of FOV along y-direction (height) 
		@return beta: angular size of FOV along x-direction (width) 
		"""

		from numpy import sqrt, arctan2, rad2deg
		f = float(focalLength)
		a = float(detectorHeight)
		b = float(detectorWidth)
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
	# findLambdaSet() Caclculates a set of wavelengths at which the QE and
	#	transmission curves will be interpolated and the planck function will
	#	be evaluated. It takes the smallest and largest wavelength values 
	#	between qe and tc and calculates a numpy arange between them with a
	#	spacing of lambdaBinSize
	###############################################################################
	def findLambdaSet(self, qe, tc, lambdaBinSize):
		"""!
		@param qe: qe dictionary input into camera object by user
		@param tc: tc dictionary input into camera object by user
		@param lambdaBinSize: lambdaBinSize input into camera object by user
		@return lambdaSet: unified set of wavelengths for tc/qe/planck
		"""		
		from numpy import concatenate, arange
		allLambdas = concatenate([qe['lam'],tc['lam']])
		minLambda = min(allLambdas)
		maxLambda = max(allLambdas)
		lambdaSet = arange(minLambda,maxLambda + lambdaBinSize,lambdaBinSize)
		return lambdaSet

	###############################################################################
	#
	# interpolateLambdaDependent() this function creates unified tc and qe curves
	#	that can be multiplied simply to create a full sensitivity curve
	#	Notes: I moved the actual funciton to util within the python framework
	#		because it was useful elsewhere, too.
	###############################################################################
	def interpolateLambdaDependent(self, ex,lambdaSet):
		"""!
		@param ex: tc or qe dictionary input into camera object by user
		@param lambdaSet: unified lambda set (output of camera.findLambdaSet)
		@return interpolatedSet: qe or tc dictionary evaluated at wavelength values
			of lambdSet
		"""
		from util import interpolateLambdaDependent
		interpolatedSet = interpolateLambdaDependent(ex,lambdaSet)
		return interpolatedSet

	###############################################################################
	#
	# calculateHotDark()
	#
	###############################################################################
	def calculateHotDark(self, detectorArray, resolutionHeight, resolutionWidth):
		"""!
		@param detector_array, 
		@param resolution_height, 
		@param resolution_width
		@return Impure_Array (Normalized) - 1-D resolution array with corrupted pixel factors
		"""
		# Module Imports
		import random
		import numpy as np
		##################

		# Gaussian Knob
		mu = 0  		# Mean
		sigma = 1e-6  	# Std Dev
		knob = abs(np.random.normal(mu, sigma, 1))  # Initial Number of Corrupt Pixels

		# Image Dimensions
		numRows = resolution_width
		numCols = resolution_height

		imageRes = np.full((numRows, numCols), 125)		# 2D Array of Pixel Value
		numPix = resolution_height * resolution_width	# Total Number of Pixels

		# Completely Impure Image
		impureRes = np.floor(np.random.rand(*imageRes.shape) * 255)

		# Integration Time (days)
		simulationTime = 10 	# Seconds (Basilisk User Input)
		timeStep = 1  			# Time Step should be per day of Basilisk Time

		# Integration Loop
		while timeStep < simulationTime:

			knob += knob + abs(np.random.normal(mu, sigma, 1))  # Knob
			if knob >= 1:
				knob = 1.0

			numImp = int(np.ceil(abs(knob) * numPix))  # Impure Pixels
			good = (numPix - numImp)  # Good Pixels

			# Boolean Mask of Impure Elements
			boolMask = np.array([0] * good + [1] * numImp).astype(np.bool)
			np.random.shuffle(boolMask)  # Shuffle
			np.random.shuffle(boolMask)  # Shuffle
			np.random.shuffle(boolMask)  # Shuffle
			mask = boolMask.reshape(resolution_width, resolution_height)  # Shuffled Mask

			# Altered Resolution
			imageRes[mask] = impureRes[mask]

			# Pixel Threshold Values
			avg 		= 125.
			hotThresh 	= 225.
			warmThresh 	= 200.
			coldThresh 	= 50.
			darkThresh 	= 25.

			# Hot Pixels
			hot = (imageRes >= hotThresh)
			imageRes[hot] = random.randint(hotThresh, 255)
			# Warm Pixels
			warm = np.logical_and(imageRes < hotThresh, imageRes > avg)
			imageRes[warm] = random.randint(warmThresh, hotThresh + 25)
			# Average Pixels
			nice = (imageRes == avg)
			imageRes[nice] = avg
			# Cold Pixels
			cold = np.logical_and(imageRes < avg, imageRes > darkThresh)
			imageRes[cold] = random.randint(darkThresh, coldThresh)
			# Dark Pixels
			dark = (imageRes <= darkThresh)
			imageRes[dark] = random.randint(0, darkThresh)

			# Integration Time Step
			timeStep += 0 + timeStep


		# Reshape 'imageRes' into 1-D Array
		impureNorm = [x / avg for x in imageRes]

		self.hotDarkArray = impureNorm

		return self.hotDarkArray

	###############################################################################
	# updateState() is the master function of the camera object. It reads
	#		takeImage from the BSK messaging system and opens, closes, or updates
	#		images appropraitely 
	###############################################################################

	def updateState(self):
		# self.hotDark = self.calculateHotDark(dt)

		#initialize a counter so we can check if we have open images
		openImages = 0
		takeImage = self.takeImage
		#loop through all images counting the ones that are open.
		#if you find an open one, set its key to openImageKey.
		#if there's more than one open, the var will be overwritten, but
		#that's OK since the code will error out later.
		for key in self.images:
			if self.images[key].imageOpen == 1: openImages +=1
			openImageKey = key

		if openImages > 1.:
			#we should never have more than one image open at a time, so
			#if openImages > 1, we have a problem. Report an error and
			#return
			print('ERROR: Camera Object Has Multiple Open Images!!!')
			return
		elif openImages == 1.:
			#if there is exactly one open image and we are still imaging
			#we need to update the state.
			if takeImage == 1:
				self.images[openImageKey].updateState()
			#if there is exactly one open image but the exec has told the
			#camera to stop imaging, just close it. We still neet do run
			#the image's updateState() method in order to get it to actually
			#create the image. It only collects attitude information until
			#this step.
			else:
				self.images[openImageKey].imageOpen = 0
				self.images[openImageKey].updateState()
				try:
					self.images[openImageKey].imgTime = self.imgTime
				except:
					self.images[openImageKey].imgTime = -1
				self.images[openImageKey].beaconPos = []
				self.images[openImageKey].beaconID = []
				self.images[openImageKey].beaconRad = []
				self.images[openImageKey].beaconAlbedo = []
				self.images[openImageKey].cameraParam = {
				        'resolution': 
				        	(self.resolutionHeight,self.resolutionWidth),
				        'focalLength': 
				        	self.focalLength,
				        'sensorSize': 
				        	(self.detectorHeight,self.detectorWidth),
				        'FOV': 
				        	(self.angularHeight,self.angularWidth),
				        'pixelSize': 
				        	(self.detectorHeight/self.resolutionHeight,
				            self.detectorWidth/self.resolutionWidth)
				    }

				self.images[openImageKey].imgBeaconPos = []
				for each in self.bodies:
					self.images[openImageKey].imgBeaconPos.append(
						each.state)

				self.images[openImageKey].imgPos = self.scState
				self.images[openImageKey].imgDCM = \
					self.body2cameraDCM.dot(self.scDCM)
				self.images[openImageKey].imgTime = self.imgTime

				try: 
					self.msg['verbose'] == 1
				except:
					self.msg['verbose'] = 0

				if not(self.msg['verbose']):
					delattr(self.images[openImageKey],'RA')
					delattr(self.images[openImageKey],'DE')
					delattr(self.images[openImageKey],'DCM')
					delattr(self.images[openImageKey],'alpha')
					delattr(self.images[openImageKey],'beta')
					delattr(self.images[openImageKey],'gamma')
					delattr(self.images[openImageKey],'scenes')
					delattr(self.images[openImageKey],'c1')
					delattr(self.images[openImageKey],'c2')
					delattr(self.images[openImageKey],'c3')
					delattr(self.images[openImageKey],'n1')
					delattr(self.images[openImageKey],'n2')
					delattr(self.images[openImageKey],'n3')
					delattr(self.images[openImageKey],'BVT')
					delattr(self.images[openImageKey],'VT')
					delattr(self.images[openImageKey],'starID')
					delattr(self.images[openImageKey],'camera')
					delattr(self.images[openImageKey],'solidAngleSubtended')
					delattr(self.images[openImageKey],'I')
					delattr(self.images[openImageKey],'T')
		else:
			if takeImage == 1.:
				#if we have no open image, count the number of open images
				#and use that as the key for the next one (python indexes at
				#zero, so if we want to index using autoincrementing numbers, 
				#this will give us the next number).
				newKey = len(self.images.keys())
				#initialize the image.
				self.images[newKey] = image(self,self.msg)
				#give the image its first update.
				self.images[newKey].updateState()
				
				#a bunch of attributes are only ever saved to the camera
				#as debugging artifacts. If verbose == 0 or it is not
				#set, then we go into the except section and destroy
				#all the unneeded artifacts. That way we don't need to
				#unnessarily blow up RAM when saving a lot of images.

		return

class image:
	"""!
	Image is a class for the image data product produced by a spacecraft
		camera.

	methods:

	User Variables:
		None. All variables in the image object are added dynamically

	Computed Variables:
		imageOpen: boolean that keeps track of if the image is currently open
		alpha: list of euler angle rotatons about axis 3 at each timestep
		beta: list of euler angle rotatons about axis 2 at each timestep
		gamma: list of euler angle rotatons about axis 1 at each timestep
		camera: camera object representing the camera this was taken by
		DCM: list of DCMs from intertial to camera frame at each timestep
		scenes: list of scene objects
		starID: list of IDs for each star that is visible at any time
			during exposure
		RA: Similar to the RA varable that is held in the camera object,
			but with only stars/facets that are in the FOC of this image. 
		DE: Similar to the DE varable that is held in the camera object,
			but with only stars/facets that are in the FOC of this image.
		n1: Similar to the n1 varable that is held in the camera object,
			but with only stars/facets that are in the FOC of this image.
		n2: Similar to the n2 varable that is held in the camera object,
			but with only stars/facets that are in the FOC of this image.
		n3: Similar to the n3 varable that is held in the camera object,
			but with only stars/facets that are in the FOC of this image.
		VT: Similar to the VT varable that is held in the camera object,
			but with only stars/facets that are in the FOC of this image.
		BVT: Similar to the BVT varable that is held in the camera object,
			but with only stars/facets that are in the FOC of this image.
		c1: 1st component of Camera frame unit vector of all stars visible
			at any time in this image
		c2: 2nd component of Camera frame unit vector of all stars visible
			at any time in this image
		c3: 3rd component of Camera frame unit vector of all stars visible
			at any time in this image
		detectorArray: array with summed intensities of light incident on
			each pixel

	Keyword Arguments:
		None

	Notes:
		-alpha, beta, gamma really need to be removed here. Although the
		Euler angle transformation is a lot easier to grasp intuitively, it
		would make way a camera-to-body DCM that is stored in the camera
		object along with an inertial to spacecraft DCM that is stored in the
		spacecraft object. This way the image object will only need camera and
		msg passed in.
		-This class ultimately needs to be split into image and frame classes
		where frames are totaled in order to make the final image. The frames
		should be calculated with a margin added to the true FOV because once
		we start to add saturation effects, it will matter a lot what is right
		outside the FOV (i.e. if the sun is right next to the FOV, the 
		background counts will be very high, even though the sun is not imaged
		directly)

	"""
	def __init__(
		self, 
		camera,
		msg,
		**kwargs
		):

		self.imageOpen = 1
		self.alpha = []
		self.beta = []
		self.gamma = []
		self.camera = camera
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
		self.detectorArray = []
		self.starID = []

	###############################################################################
	# updateState() is the master function of the image object. It does one of
	#		two things. If self.camera.msg['takeImage'] == 1, it will save the
	#		current inertial to camera DCM to the DCM list. if  
	#		self.camera.msg['takeImage'] == 0, then it will find all stars and
	#		bodies that are visible at any time in the exposure and add them to 
	#		attirbutes to be used when making scenes. Next, it takes each DCM,
	#		calculates which stars and bodies are visible at each time step
	#		and creates a scene.
	###############################################################################
	def updateState(self):

		if self.camera.takeImage== 1:
			#if the image is still actively being taken, all we do is save
			#off a DCM for that scene. The image will only be created once
			#self.camera.msg['takeImage'] is set to zero and
			#image.updateState() is run.
			self.DCM.append(self.camera.body2cameraDCM.dot(
				self.camera.scDCM))
		else:
			from scipy.linalg import inv
			from numpy import arctan2, arcsin, array, zeros, pi, floor
			from em import planck
			import pdb

			#this finds a subset of stars that are in the FOV at some point
			#during the exposure. Some functionality of findStarsInFOV()
			#is not needed in this step (like converting to pixel/line), 
			#but it is done anyway so we can reuse the function.

		
			#use the first attitude of the exposure as the central attitude
			#for calculating subset of stars.
			DCM_0 = self.DCM[0]
			DCM_0_inv = inv(DCM_0)
			alpha_0 = arctan2(DCM_0[0,1],DCM_0[0,0])
			beta_0 = -arcsin(DCM_0[0,2])
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

			#When we call findStarsInFOV()
			#to find all the stars in the exposure at any time, we want
			#to find all physical objects. This way we only render bodies
			#once per exposure. All the rest is just filtering them out
			#per scene and adding psf/noise characteristics of the image.
			fullExposureMsg = dict(self.camera.msg)
			fullExposureMsg['psf'] = 0
			fullExposureMsg['raster'] = 0
			fullExposureMsg['photon'] = 0
			fullExposureMsg['dark'] = 0
			fullExposureMsg['read'] = 0

			FOV = self.findStarsInFOV(
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
				self.camera.angularDiagonal/2,
				self.camera.angularDiagonal/2,
				self.camera.resolutionHeight,
				self.camera.resolutionWidth,
				self.camera.RA, self.camera.DE, 
				self.camera.n1, self.camera.n2, self.camera.n3,
				self.camera.VT, self.camera.BVT, self.camera.T,
				self.camera.T, #spoof so the fcn won't break :-(
				self.camera.solidAngleSubtended,
				self.camera.maxMag, 
				self.camera.starID,
				fullExposureMsg
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
			self.starID = FOV['starID']
			self.solidAngleSubtended = FOV['solidAngleSubtended']
			I = []
			import matplotlib.pyplot as plt


			for i in range(0,len(self.T)): 

				T = self.T[i]
				if T == 5778: 
					flux_per_m2_per_nm_per_sr_at_star = self.camera.solarBB
				else:
					flux_per_m2_per_nm_per_sr_at_star = planck(T,self.camera.lambdaSet*1e-9)
				solidAngleSubtended = self.solidAngleSubtended[i]
				
				flux_per_m2_per_nm_per_sr_at_obs = flux_per_m2_per_nm_per_sr_at_star*solidAngleSubtended
				flux_per_m2_per_nm = pi*flux_per_m2_per_nm_per_sr_at_obs
				flux_per_m2 = flux_per_m2_per_nm*self.camera.lambdaBinSize
				flux = flux_per_m2*self.camera.effectiveArea
				photons_per_sec = flux/self.photonEnergy(self.camera.lambdaSet*1e-9)
				electrons_per_sec = photons_per_sec*self.camera.sensitivityCurve
				I.append(sum(electrons_per_sec))
			
			self.I = array(I)
			self.starID = FOV['starID']


			#since the camera object removes occulted objects
			#and adds body renderings, we don't need to do that
			#all again here, so we force those parts of the msg
			# to be zero.
			sceneMsg = dict(self.camera.msg) 
			sceneMsg['rmOcc'] = 0
			sceneMsg['addBod'] = 0

			if self.camera.msg['psf']:
				psf = self.psf(self.camera.psfSigma)


			#create one scene per DCM we collected above
			for i in range(0,len(self.DCM)):
				FOV = self.findStarsInFOV(
					self.alpha[i] + alpha_0, 
					self.beta[i] + beta_0, 
					self.gamma[i] + gamma_0, 
					0,
					0,
					0,
					0,
					self.camera.angularHeight/2,
					self.camera.angularWidth/2,
					self.camera.resolutionHeight,
					self.camera.resolutionWidth,
					self.RA, self.DE, 
					self.n1, self.n2, self.n3,
					self.VT, self.BVT, self.T,
					self.I,
					self.solidAngleSubtended,
					self.camera.maxMag, 
					self.starID,
					sceneMsg
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
						FOV['starID']
						)
					)
			i = 0

			for eachScene in self.scenes:
				i+=1
				if self.camera.msg['psf']:
					pixel = psf['x'].reshape(len(psf['x']),1) + eachScene.pixel
					line = psf['y'].reshape(len(psf['y']),1) + eachScene.line
					I = psf['I'].reshape(len(psf['I']),1)*eachScene.I

					eachScene.psfPixel = pixel.reshape(len(eachScene.pixel)*len(psf['x']))
					eachScene.psfLine = line.reshape(len(eachScene.line)*len(psf['y']))

					eachScene.psfI = I.reshape(len(psf['I'])*len(eachScene.I))


				#self.addBkgdStarLight()
				if self.camera.msg['photon']:
					eachScene.psfI = self.addPoissonNoise(eachScene.psfI)

				if self.camera.msg['raster']:
					eachScene.detectorArray = \
						self.rasterize(
							self.camera.resolutionWidth,
							self.camera.resolutionHeight,
							eachScene.psfPixel,
							eachScene.psfLine,
							eachScene.psfI
							)*self.camera.dt

					if self.camera.msg['dark']:
						eachScene.detectorArray = \
							eachScene.detectorArray + \
							self.addPoissonNoise(
								zeros(len(eachScene.detectorArray)) + \
								self.camera.darkCurrent
								)

			self.detectorArray = 0
			for eachScene in self.scenes:
				self.detectorArray += eachScene.detectorArray

			if self.camera.msg['read']:
				self.detectorArray = self.addReadNoise(
						self.detectorArray,
						self.camera.readSigma
						)

			self.detectorArray[self.detectorArray < 0] = 0
			
			#convert phton counts to DN
			self.detectorArray = floor(self.detectorArray/self.camera.dnBinSize)

			#cut off pixels that are saturated
			self.detectorArray[self.detectorArray > self.camera.dnDepthMax] = self.camera.dnDepthMax

	########################################################################################################################
	########################################################################################################################
			if self.camera.msg['hotDark']:
				hd = self.camera.calculateHotDark(
													self.camera.hotDarkArray,
													self.camera.resolutionHeight,
													self.camera.resolutionWidth
													)
				import numpy as np
				hot_dark = np.asarray(hd).flatten()
				self.camera.hotDarkArray = hotDark
				self.detectorArray = self.detectorArray * hotDark


			self.detectorArray = self.detectorArray.reshape(
				self.camera.resolutionWidth,
				self.camera.resolutionHeight
				)

	###########################################################################
	# findStarsInFOV() finds the pixel and line coordinates of a 
	###########################################################################
	def findStarsInFOV(
		self, 
		alpha, 
		beta, 
		gamma, 
		alphaMax,
		betaMax,
		alphaMin,
		betaMin,
		halfAlpha,
		halfBeta,
		alphaResolution, 
		betaResolution,
		RA, DE, 
		n1, n2, n3,
		VT, BVT,T,I,
		solidAngleSubtended,
		maxMag, 
		starIDs,
		msg,
		**kwargs):

		"""!
		@param alpha: euler angle rotation about axis 3 at t0
		@param beta: euler angle rotation about axis 2 at t0
		@param gamma: euler angle rotation about axis 1 at t0
		@param alphaMax: largest  alpha at any time
		@param betaMax: largest beta any time
		@param alphaMin: smallest  alpha at any time
		@param betaMin: smallest beta any time
		@param halfAlpha: Half the angular height of the camera
		@param halfBeta: Half the angular width of the camera
		@param alphaResolution: resolution in the line dimension of the camera
		@param betaResolution: resolution in the pixel dimension of the camera
		@param RA: RA of all stars that may be in FOV
		@param DE: DE of all stars that may be in FOV
		@param n1: inertial coord n1 of all stars that may be in FOV
		@param n2: inertial coord n2 of all stars that may be in FOV
		@param n3: inertial coord n3 of all stars that may be in FOV
		@param VT: Tycho visual magnitude of all stars that may be in FOV
		@param BVT: Tycho color index of all stars that may be in FOV
		@param T: Computed temperature of all stars that may be in FOV
		@param I: Computed incident intensity of all stars that may be in FOV
		@param solidAngleSubtended: solidAngleSubtended of all stars that may be in FOV
		@param maxMag: maximum magnitude of camera
		@param starIDs: starIDs of all stars that may be in FOV
		@param msg: debug message passed from camera
		@return pixel: pixel coord of all stars in FOV
		@return line: line cood of all stars in FOV
		@return RA: RA of all stars in FOV
		@return DE: DE of all stars in FOV
		@return n1: inertial coord n1 of all stars in FOV
		@return n2: inertial coord n2 of all stars in FOV
		@return n3: inertial coord n3 of all stars in FOV
		@return VT: Tycho visual magnitude of all stars in FOV
		@return BVT: Tycho color index magnitude of all stars in FOV
		@return c1: camera frame coord c1 of all stars that may be in FOV
		@return c2: camera frame coord c2 of all stars that may be in FOV
		@return c3: camera frame coord c3 of all stars that may be in FOV
		@return starID : starID of all stars in FOV
		@return I: Computed incident intensity of all stars in FOV
		@return T: T of all stars in FOV
		@return solidAngleSubtended': solidAngleSubtended of all stars in FOV
		"""
		from numpy import array, deg2rad, sin, cos, append, sqrt, zeros, ones, logical_and
		from numpy.linalg import norm

		if len(self.camera.bodies) != 0:

			n = 1	
			bodies = self.camera.bodies

			#calculate how far each body is from the sc
			for body in bodies: body.distFromSc = norm(body.state[0:3] - self.camera.scState[0:3])
			#sort bodies by how far they are from the sc
			#this needs to be done 

			sortedBodies = list(bodies)
			sortedBodies.sort(key=lambda x:x.distFromSc, reverse=True)

			for body in sortedBodies:
				n+=1
				if msg['rmOcc']:
					occCheck = self.removeOccultations(body,n1,n2,n3)
					n1 = n1[occCheck]
					n2 = n2[occCheck]
					n3 = n3[occCheck]
					RA = RA[occCheck]
					DE = DE[occCheck]
					starIDs = starIDs[occCheck]
					T = T[occCheck]
					solidAngleSubtended = solidAngleSubtended[occCheck]
					I = I[occCheck]


				if msg['addBod']:
					from lightSimFunctions import lightSim
					DCM = self.camera.body2cameraDCM.dot(self.camera.scDCM)
					facets = lightSim(
						DCM,
						self.camera.scState[0:3],
						body.state[0:3],
						(
							self.camera.angularHeight,
							self.camera.angularWidth
							),
						100,
						100,
						False,
						body.albedo,
						body.r_eq,
						body.id
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
					surf_n1 -= self.camera.scState[0]
					surf_n2 -= self.camera.scState[1]
					surf_n3 -= self.camera.scState[2]

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
					starIDs = append(starIDs, zeros(len(surf_n1)))

					if len(surf_n1) > 1:
						solidAngleSubtended = append(
							solidAngleSubtended,
							facets['netAlbedo']*\
							facets['facetArea'])
						I = append(I,
							facets['netAlbedo']*\
							facets['facetArea'])
					else:
						solidAngleSubtended = append(
							solidAngleSubtended,
							sum(facets['netAlbedo']*\
							facets['facetArea']))
						I = append(I,
							sum(facets['netAlbedo']*\
							facets['facetArea']))

		c2_max = alphaMax + sin(deg2rad(halfAlpha))
		c2_min = alphaMin - sin(deg2rad(halfAlpha))
		c3_max = betaMax + sin(deg2rad(halfBeta))
		c3_min = betaMin - sin(deg2rad(halfBeta))

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
		ind = abs(c2/c1*self.camera.focalLength) < self.camera.detectorWidth/2
		#Remove stars outside the FOV in the c3 direction
		ind = logical_and(ind,abs(c3/c1*self.camera.focalLength) < self.camera.detectorHeight/2)
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
		starIDs = starIDs[ind]
		T = T[ind]
		solidAngleSubtended = solidAngleSubtended[ind]
		I = I[ind]

		#using similar triangles
		pixel = -self.camera.focalLength*c2/c1*\
			betaResolution/self.camera.detectorWidth + \
			betaResolution/2

		line = -self.camera.focalLength*c3/c1*\
			alphaResolution/self.camera.detectorHeight + \
			alphaResolution/2

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
			'starID' : starIDs,
			'I': I,
			'T': T,
			'solidAngleSubtended': solidAngleSubtended
		}
	###########################################################################
	#
	# removeOccultations() is a method to remove stars and beacon facets
	#	that are behind beacons
	#
	###########################################################################

	def removeOccultations(self,body,n1,n2,n3):
		"""!
		@param body: Object. Instantiation of the body class from bodies.py.
			Variables used from that class include equatorial radius,
			polar radius, and state (for the position to the body).

		@param n1, n2, n3: Numpy Float Arrays. All should be the same size. Each
			star/body fact that is still in the FOV at this point should
			have an entry in each, and each index MUST correspond to data
			from the same star/facet!
		@return occCheck: Numpy Array. A list of booleans, each corresponding to 
		"""
		from numpy import array, stack, einsum, logical_or
		from numpy.linalg import norm
		#this needs to be multiplied by a transformation matrix in order
		#to account for planetary attitude, but I haven't thought through
		#that just yet.
		try:
			r_eq = body.r_eq
			r_pole = body.r_pole
		except:
			r_eq = body.r_eq
			r_pole = body.r_eq
		A = array([
			[r_eq**-2,0,0],
			[0,r_eq**-2,0],
			[0,0,r_pole**-2]
			])

		v = stack([n1,n2,n3])
		#agnostic to coordinate frame origin. Just needs same axis
		#directions as v.
		y_0 = self.camera.scState[0:3] - body.state[0:3]
		#1xn
		vTAy_0 = (v.T.dot(A).dot(y_0)).T
		#1xn
		vTAv = einsum('ij,ji->i',v.T,A.dot(v))
		#scalar
		y_0TAy_0 = y_0.T.dot(A).dot(y_0)


		discriminant = vTAy_0**2 - vTAv*(y_0TAy_0 - 1)
		occCheck = logical_or(discriminant < 0,v.T.dot(y_0) > 0)
		return occCheck


	###########################################################################
	#
	# gaussian() calculates the value of a gaussian r away from the center
	#	with standard deviation sigma. Assumes covariance is a multiple of
	#	identity
	###########################################################################
	def gaussian(self,r,sigma):
		"""!
		@param r: distance of point in question from center of gaussian
		@param sigma: standard deviation of gaussian
		@return I: Intensity of gaussian
		"""

		from numpy import exp
		I = exp(-r**2/(2*sigma**2))
		return I

	###########################################################################
	#
	# psf() computes the point spread function for point sources in an image
	#
	###########################################################################
	def psf(self,sigma):
		"""!
		@param sigma: standard deviation of gaussian point spread.
		@return psfDict: dictionary with x, y, and I values for PSF. x and y
			convert to pixel/line, and I is intensity.
		"""
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
			I = append(I,ones(len(theta))*self.gaussian(r,sigma))
		I = I/sum(I)

		psfDict = { "x":x, "y":y, "I":I }

		return psfDict

	###########################################################################
	#
	# rasterize() floors the pixel and line coordinates and the uses pandas
	#		to sum all intensity that falls in the same bin.
	#
	###########################################################################
	def rasterize(self,pixelResolution,lineResolution,pixelCoord,lineCoord,intensity):
		"""!
		@param pixelResolution: number of pixels in the width dimension of the
			detector array
		@param lineResolution: number of pixels in the height dimension of the
			detector array
		@param pixelCoord: x (pixel) coordinate of every point source in scene
		@param lineCoord: y (line) coordinate of every point source in scene
		@param intensity: incident intensity of every point source in scene

		@return detectorArray: array with summed intenisty for every pixel in
			the detector array
		"""

		from numpy import floor, zeros, array, arange, append 
		from numpy import concatenate, logical_and
		from pandas import DataFrame

		#adding PSF introduces some values that are not on the detector. Remove them here

		positiveCoords = logical_and(pixelCoord > 0, lineCoord > 0)
		pixelCoord = pixelCoord[positiveCoords]
		lineCoord = lineCoord[positiveCoords]
		intensity = intensity[positiveCoords]

		notTooBig = logical_and(pixelCoord < pixelResolution, lineCoord < lineResolution)
		pixelCoord = pixelCoord[notTooBig]
		lineCoord = lineCoord[notTooBig]
		intensity = intensity[notTooBig]

		intPixCoord = floor(pixelCoord).astype(int)
		intLineCoord = floor(lineCoord).astype(int)

		detectorPosition = (lineResolution*intLineCoord + intPixCoord)
		detectorPosition = append(detectorPosition,arange(pixelResolution*lineResolution))
		intensity = append(intensity,zeros(pixelResolution*lineResolution))

		data = concatenate([detectorPosition,intensity])
		data = data.reshape(2,int(len(data)/2)).T
		df = DataFrame(data,columns = ["Position","Intensity"])
		detectorArray = df.groupby("Position").sum().values.T[0]

		return detectorArray

	###########################################################################
	#
	# addReadNoise() adds gaussian read noise. It is called once per image.
	#
	###########################################################################
	def addReadNoise(self,detectorArray,sigma):
		"""!	
		@param detectorArray: array that read noise is to be added to
		@param sigma: standard deviation of the gaussian noise to be added
		@return detectorArray: detectorArray with read noise added
		"""
		from numpy.random import randn
		detectorArray = detectorArray + sigma*randn(len(detectorArray))
		return detectorArray

	###########################################################################
	#
	# addPoissonNoise() adds poisson noise
	#
	###########################################################################
	def addPoissonNoise(self, I):
		"""!
		@param I: array of intensities to add poisson nouse to
		@return I: intensity array with poisson noise added
		"""
		from numpy.random import poisson
		I = I + poisson(I)
		return I

	###########################################################################
	#	planck() is a function that calculates a planck blackbody function
	#
	###########################################################################
	def planck(self, T,lam):
		"""!
		@param T: temperature of the star in Kelvin
		@param lam: wavelength bins at which to calculate in METERS(!)
		@return I: intensity of light at each wavelengh bin in W/m^2/nm/sr
		"""
		from constants import h, c, k_B
		from numpy import exp

		topPart = 2*h*c**2
		bottomPart = lam**5*(exp(h*c/(lam*k_B*T))-1)

		I = topPart/bottomPart*1e-9 #1e-9 to convert to per nm

		return I

	###########################################################################
	#	stefanBoltzmann() is a function that total flux from a star given 
	#
	#	Notes:
	#		This function isn't used in DINO C-REx except as part of the test
	#		for planck(). The Planck function integrated over all wavelength
	#		space should be identically equal to the Stephan-Boltzmann function
	#
	###########################################################################
	def stefanBoltzmann(self, T):
		"""!
		@param T: temperature of the star in Kelvin
		@return F: total flux at stellar surface in W/m^2
		"""
		from constants import sigma
		return sigma*T**4

	###########################################################################
	#	photonEnergy() calculates the energy per photon at a wavelength lam.
	#		it is used to convert how many photons are incident on the 
	#		detector per joule of energy
	#
	###########################################################################

	def photonEnergy(self, lam):
		"""!
		@param lam: array of wavelength values to evaluate
		@return photonEnergy: energy of a photon at each lam value
		"""
		from em import photonEnergy
		return photonEnergy(lam)
		from constants import h, c
		photonEnergy = h*c/lam
		return photonEnergy

	###########################################################################
	#	addZodiacalLight() 
	#
	#	Notes: the data represented in pioneer10BkdgStarLight comes from
	#	Lienert's 1997 reference on diffuse night sky brightness on p39.
	# 
	#	Units for the array are S10_sun, given by leinart to be 1.28e-8*.156
	#	W/m^2/sr
	#
	#	Rows of the array represent bins of declination.
	#
	###########################################################################
	def addZodiacalLight(self):
		"""!
		@param attitudeDCM: DCM for the inertial to camera attitude.
		@return zodiacalFlux: 
		"""
		from numpy import array, identity, flip, vstack, hstack, transpose
		from numpy import arctan2, pi, rad2deg, arcsin, argmin
		sunAngleRad = arctan2(self.camera.scDCM[1],self.camera.scDCM[0]) + pi
		sunAngleDeg = rad2deg(sunAngleRad)

		raBins = array([
			0,5,10,15,20,25,30,35,40,45,60,75,90,105,120,135,150,165,180
			])
		decBins = array([
			0,5,10,15,20,25,30,45,60,75
			])
		C12 = self.camera.scDCM[0,1]
		C11 = self.camera.scDCM[0,0]
		C13 = self.camera.scDCM[0,2]

		camRArad = arctan2(C12,C11)
		camDErad = -arcsin(C13)

		camRAdeg = abs(rad2deg(camRArad) - sunAngleDeg)
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

		zodiacalLight = zodiacalLight*1.18e-8*.156

		return zodiacalLight[camRAbin,camDEbin]

	###########################################################################
	#	addBkgdStarLight() 
	#	Notes: the data represented in pioneer10BkdgStarLight comes from
	#	Lienert's 1997 reference on diffuse night sky brightness on p78 and 79.
	# 
	#	Units for the array are S10_sun, given by leinart to be 1.28e-8*.156
	#	W/m^2/sr
	#
	#	Rows of the array represent bins of declination, and are spaced by 10
	#	degrees. Columns are binned by right ascension and are spaced every 10
	#	degrees.
	###########################################################################
	def addBkgdStarLight(self):
		"""!	
		@param attitudeDCM: DCM for the inertial to camera attitude.
		
		@return bkgdFlux: flux to be applied to every pixel in the detector array
		@return representing the brightness of unresolved stars
		"""
	

		from numpy import array, arcsin, arctan2, rad2deg, deg2rad, flip

		decBins = array([-80, -70, -60, -50, -40, -30, -24, -20, -16, -12, -8, 
			-4, 0,4, 8, 12, 16, 20, 24, 30, 40, 50, 60, 70, 80])
		raBins = array([0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,
			160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,
			320,330,340,350])

		C12 = self.camera.body2cameraDCM.dot(self.camera.scDCM)[0,1]
		C11 = self.camera.body2cameraDCM.dot(self.camera.scDCM)[0,0]
		C13 = self.camera.body2cameraDCM.dot(self.camera.scDCM)[0,2]

		camRArad = arctan2(C12,C11)
		camDErad = -arcsin(C13)

		camRAdeg = rad2deg(camRArad)
		camDEdeg = rad2deg(camDErad)
		if camDEdeg < 0: camDEdeg
		if camRAdeg < 0: camRAdeg += 360


		#the place where the difference between camera DE and the decliation bins
		#given by Leinert is the smallest will be the delination bin that we want.
		camDEbin = min(abs(camDEdeg-decBins)) == abs(camDEdeg-decBins)
		camRAbin = min(abs(camRAdeg-raBins)) == abs(camRAdeg-raBins)



		pioneer10BkdgStarLight = array([
			[ 68,  55,  55,  56,  49,   35,  39,  40,  38,  42,  49,  46,  50,
			  56,  43,  51,  59,  56,  61,  71, 116, 192, 284, 162, 124],                  
			[ 74,  74,  42,  38,  36,   35,  34,  43,  39,  34,  44,  45,  35,
			  48,  47,  51,  55,  52,  63,  64, 108, 181, 261, 185, 121],
			[ 78,  67,  51,  44,  40,   31,  36,  40,  39,  39,  43,  38,  32,
			  53,  45,  40,  56,  58,  57,  79,  99, 158, 250, 176, 106],
			[ 78,  53,  48,  44,  39,   38,  34,  46,  45,  41,  40,  51,  43,
			  47,  64,  56,  61,  56,  66,  82,  97, 187, 265, 157, 110],
			[ 73,  69,  49,  42,  36,   36,  45,  42,  40,  44,  45,  48,  53,
			  49,  59,  54,  61,  62,  60,  74, 116, 177, 240, 150, 109],
			[ 72,  63,  56,  41,  48,   50,  44,  41,  44,  51,  49,  56,  60,
			  59,  61,  70,  74,  79,  78,  85, 141, 168, 205, 140, 111],
			[ 77,  71,  53,  52,  48,   50,  51,  52,  52,  57,  60,  68,  68,
			  73,  88,  85,  82,  90, 103, 108, 133, 183, 175, 128, 117],
			[ 81,  90,  72,  62,  54,   58,  66,  64,  64,  74,  76,  75,  87,
			  82,  93,  95, 112, 122, 102, 104, 189, 164, 161, 137, 103],
			[ 81, 341,  64,  72,  64,   63,  80,  84,  90, 115, 134, 131, 136,
			 135, 135, 125, 146, 143, 130, 179, 214, 284, 125, 115,  85],
			[ 88, 226,  88,  93,  83,   90, 117, 123, 138, 147, 155, 158, 174,
			 174, 220, 238, 258, 218, 246, 227, 169, 137, 112,  91,  93],
			[ 98, 107, 114, 116, 135,  176, 193, 248, 240, 248, 225, 253, 283,
			 280, 284, 263, 229, 233, 207, 210, 145, 119, 108,  77,  84],
			[102, 111, 141, 182, 218,  315, 331, 321, 297, 312, 295, 267, 225,
			 206, 182, 184, 165, 135, 138, 105, 108,  93,  81,  73,  77],
			[110, 141, 176, 275, 308,  385, 327, 288, 249, 200, 166, 178, 173,
			 167, 155,  -1,  -1,  -1,  -1, 105,  75,  63,  75,  79,  73],
			[113, 142, 224, 337, 272,  274, 201, 156, 137,  -1,  -1,  -1,  -1,
			  -1,  -1,  -1,  -1,  -1,  -1,  -1,  77,  52,  57,  74,  70],
			[119, 206, 372, 294, 398,  175, 114, 116, 126, 109, 101,  -1,  -1,
			  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  53,  57,  80,  65],
			[116, 224, 490, 334, 193,  120, 104,  85,  72,  68,  65,  65,  66,
			  67,  68,  -1,  -1,  -1,  -1,  -1,  -1,  52,  55,  57,  55],
			[116, 265, 846, 250, 156,  100,  75,  74,  62,  65,  57,  68,  71,
			  71,  65,  62,  59,  60,  58,  55,  52,  46,  52,  46,  49],
			[141, 308, 744, 229, 116,   88,  76,  59,  63,  61,  56,  52,  55,
			  60,  68,  66,  58,  52,  52,  50,  60,  36,  51,  54,  52],
			[145, 340, 632, 230,  99,   78,  61,  61,  57,  58,  54,  56,  54,
			  54,  41,  47,  50,  31,  45,  38,  37,  43,  45,  38,  52],
			[137, 335, 450, 214, 120,   72,  67,  59,  63,  54,  51,  42,  43,
			  55,  45,  41,  41,  36,  41,  33,  37,  35,  42,  50,  51],
			[145, 410, 565, 262, 115,   91,  77,  59,  53,  58,  51,  49,  44,
			  48,  40,  34,  39,  42,  38,  38,  34,  37,  45,  50,  58],
			[154, 351, 516, 256, 144,   95,  73,  68,  66,  71,  58,  52,  47,
			  50,  49,  44,  57,  44,  40,  40,  40,  43,  43,  48,  59],
			[153, 337, 559, 284, 173,  102,  87,  81,  80,  77,  54,  53,  58,
			  48,  53,  47,  42,  47,  54,  38,  35,  44,  40,  54,  61],
			[156, 310, 538, 481, 249,  134, 110,  88,  97,  79,  77,  62,  71,
			  60,  54,  50,  52,  43,  46,  41,  55,  44,  49,  61,  69],
			[146, 262, 537, 357, 242,  200, 152, 137, 106, 101,  82,  83,  77,
			  77,  64,  64,  54,  54,  55,  53,  47,  51,  53,  60,  68],
			[128, 213, 373, 571, 485,  366, 217, 145, 141, 145, 108, 110, 101,
			  89,  87,  79,  75,  68,  67,  59,  57,  59,  56,  65,  71],
			[127, 190, 284, 494, 550,  571, 344, 356, 301, 172, 144, 119, 136,
			 120, 123, 119, 116,  90,  85,  89,  81,  66,  64,  78,  75],
			[115, 155, 205, 284, 626, 1393, 649, 498, 474, 261, 199, 162, 204,
			 225, 201, 166, 157, 148, 132, 120,  94,  96,  73,  82,  76],
			[101, 106, 150, 135, 262,  461, 481, 534, 542, 515, 509, 343, 210,
			 254, 382, 341, 237, 218, 206, 176, 145, 111,  96,  91,  81],
			[103, 113, 129, 191, 141,  183, 211, 246, 272, 274, 260, 257, 331,
			 366, 374, 271, 283, 250, 327, 387, 228, 168, 117,  90,  89],
			[102, 107,  93,  99, 112,  112, 117, 134, 144, 143, 162, 149, 178,
			 215, 237, 323, 366, 407, 316, 421, 403, 245, 163, 102,  93],
			[ 93,  88,  78,  80,  77,   81,  89,  85,  89,  91,  91, 109, 128,
			 138, 145, 174, 168, 227, 240, 288, 280, 292, 187, 133,  94],
			[ 87,  81,  72,  65,  64,   73,  70,  73,  72,  77, 100,  75,  85,
			  93, 109, 107, 131, 119, 153, 192, 279, 286, 225, 136, 105],
			[ 75,  76,  56,  54,  60,   42,  66,  60,  56,  67,  67,  63,  81,
			  78,  78, 103,  85,  93, 105, 142, 188, 384, 259, 167, 119],
			[ 66,  69,  57,  52,  51,   53,  43,  48,  58,  64,  57,  47,  49,
			  67,  65,  58,  63,  76,  76,  97, 171, 268, 241, 178, 122],
			[ 66,  64,  59,  43,  46,   56,  48,  49,  43,  50,  44,  47,  62,
			  48,  55,  59,  61,  69,  72,  84, 104, 203, 256, 166, 110] 
			])


		bkgdFlux = pioneer10BkdgStarLight[camRAbin,camDEbin][0]*1.28e-8*.156

		return bkgdFlux
class scene:
	"""!
	scene is a class for collecting each integration step during an exposure.
	These will then be summed to make frames, which are the physically relevant
	unit of an image. frames are then summed to create the full images.

	User Variables:
		None. All variables in the image object are added dynamically

	Computed Variables:
		starIDs: Similar to the starID varable that is held in the camera 
			object, but with only stars/facets that are in the FOC of 
			this scene. 
		RA: Similar to the RA varable that is held in the camera object,
			but with only stars/facets that are in the FOC of this scene. 
		DE: Similar to the DE varable that is held in the camera object,
			but with only stars/facets that are in the FOC of this scene.
		n1: Similar to the n1 varable that is held in the camera object,
			but with only stars/facets that are in the FOC of this scene.
		n2: Similar to the n2 varable that is held in the camera object,
			but with only stars/facets that are in the FOC of this scene.
		n3: Similar to the n3 varable that is held in the camera object,
			but with only stars/facets that are in the FOC of this scene.
		VT: Similar to the VT varable that is held in the camera object,
			but with only stars/facets that are in the FOC of this scene.
		BVT: Similar to the BVT varable that is held in the camera object,
			but with only stars/facets that are in the FOC of this scene.
		c1: Similar to the c1 varable that is held in the image object,
			but with only stars/facets that are in the FOC of this scene.
		c2: Similar to the c2 varable that is held in the image object,
			but with only stars/facets that are in the FOC of this scene.
		c3: Similar to the c3 varable that is held in the image object,
			but with only stars/facets that are in the FOC of this scene.
		I: Similar to the I varable that is held in the image object,
			but with only stars/facets that are in the FOC of this scene.
		pixel: Similar to the pixel varable that is held in the image object,
			but with only stars/facets that are in the FOC of this scene.
		line: Similar to the line varable that is held in the image object,
			but with only stars/facets that are in the FOC of this scene.
		detectorArray: array of summed intensities on each pixel during
			the integration timestep of this scene. Does not include noise.
			Added by image.updateState()
	"""
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
		starIDs,
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
		self.starIDs = starIDs

class beacon:
	"""!
	Beacon is a class to carry information about opnav beacons.
	All data is set by the user, ane instances of the class are really
	just carried around to help organize data.

	User Variables:
		None. All variables in the image object are added dynamically

	Computed Variables:
		None.
	"""
	def __init__(self):
		self.state = -1
		self.r_eq = -1
		self.id = -1
		self.albedo = -1


