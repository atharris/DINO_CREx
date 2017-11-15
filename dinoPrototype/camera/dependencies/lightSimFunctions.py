#	Title   : lightSimFunctions.py
#	Author  : Joe Park
#	Date    : 03/19/17
#	Synopsis: Functions for lighting simulation module

import math
import numpy as np
# import lightSimPlots as lSP
import matplotlib.pyplot as plt
import csv


################################################
################################################


def getCBparam(ind):
    """Return value of albedo and radius given index of celestial body
    Depends on CELESTIAL_BODIES_PARAMETERS.csv file in same directory.
    :param index: index number of celestial body
    :return Geometric Albedo Value, Equatorial Radius [m]"""

    indchk = 0
    with open("CELESTIAL_BODY_PARAMETERS.csv", 'rb') as f:
        reader = csv.reader(f)
        for row in reader:
            if indchk == ind:
                nameCBparam = row[0]
                albedoCBparam = float(row[1])
                radiusCBparam = float(row[2])
            indchk += 1

    return albedoCBparam, radiusCBparam, nameCBparam


##################################################
##################################################


def mapSphere(latRes, longRes, radCB):
    """Generate latitude and longitude coordinates for semi-spherical mapping. Method assumes equal divisions in
    latitude and equidistant divisions in longitude (number of longitudinal pts at each latitude dependent on horizontal
    radius of celestial body at a given latitude).
    Assumptions: surface mapping is fine enough to approximate surface area as 2D rectangle
    :param latRes: number of latitude points
    :param longRes: number of surface map pts along zero latitude
    :param radCB: radius of celestial body [m]
    :return N x 2 array of latitude and longitude coordinate for semi-spherical mapping [degrees]
            2D Surface area of a single rectangular facet [m^2]"""
    from numpy import zeros, cos, deg2rad, rad2deg
    from datetime import datetime
    start_ms = datetime.now()
    # delta latitude
    deltLat = (180./latRes)

    # delta longitude at equator and distance between surface mapping pts
    deltLongEquator = (180./longRes)
    deltHoriz = radCB*math.radians(deltLongEquator)

    # only using semi-sphere
    ptsLat = np.arange(-90.,90.+deltLat,deltLat)
    ptsLatLong = np.empty((0,2),float)
    ptsLatLongC = np.empty((0,2),float)

    latC = ptsLat[abs(ptsLat) < 90] - .5 * deltLat
    rCB = radCB*cos(deg2rad(latC))
    deltLong = rad2deg((deltHoriz/rCB))
    nptsLong, longRemainder = divmod(180., deltLong)

    # pdb.set_trace()
    # loop through each latitude and calculate longitudinal pts
    # (builds array of latitude and longitude pts of facet vertices)
    for currentLat in np.nditer(ptsLat):

        currentLatC = currentLat - .5 * deltLat

        if abs(currentLatC) < 90:

            # horizontal semicircle perimeter at current latitude
            current_rCB = radCB*math.cos(math.radians(currentLatC))

            # calculate number of surface map pts at current longitude and remainder
            currentDeltLong = math.degrees((deltHoriz/current_rCB))
            currentNptsLong, longRemainder = divmod(180., currentDeltLong)
            longI = longRemainder/2.
            longF = longI + (currentNptsLong*currentDeltLong)
            # pdb.set_trace()
            if currentNptsLong == 0:
                currentPtsLongC = np.array([0.0])
            else:
                currentPtsLongC = np.arange(
                    longI+(currentDeltLong/2.), 
                    longF, currentDeltLong
                    )

            # create array of latlong for facet vertices

            tmpStack = np.vstack(
                (
                    zeros(len(currentPtsLongC)) + currentLat,
                    currentPtsLongC
                    )
                ).T

            ptsLatLong = np.vstack((ptsLatLong,tmpStack))
    
    ptsLatLongC = ptsLatLong
    ptsLatLongC[:,0] -= .5 * deltLat

    # assume each facet has equal areas totaling seim-sphere
    nptsLatLongC = ptsLatLongC.shape[0]
    facetArea = (4.*math.pi*radCB**2) / nptsLatLongC

    mkplots = False
    if mkplots:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(ptsLatLong[:, 0], ptsLatLong[:, 1], c='b', marker='.')
        ax.scatter(ptsLatLongC[:, 0], ptsLatLongC[:, 1], c='y', marker='*')
        ax.set_xlabel('Latitude (deg)')
        ax.set_ylabel('Longitude (deg)')
        plt.grid(True)
        plt.show()
    return ptsLatLongC, facetArea

###################################################
###################################################


def checkFoV(posCB, posObs, dcmBN, fov, radiusCB):
    """Check if celelstial body is in camera field of view
    :param posCB: position of celestial body in heliocentric coordinates [m]
    :param posObs: position of camera in heliocentric coordinates [m]
    :param dcmNB: direction cosine matrix of observer attitude [body to heliocentric coord. frame]
    :param field of view tuple (horizontal, vertical) [degrees]
    :return: true if in field of view, false otherwise
    """
    from datetime import datetime
    from numpy.linalg import norm
    start = datetime.now()

    # compute relative position of celestial body to observer in body coordinates
    posObs2cbHelio = (posCB-posObs)
    posObs2cbBody = np.matmul(dcmBN, posObs2cbHelio)
    eObs2cbBody = posObs2cbBody / np.linalg.norm(posObs2cbBody)

    #compute angular size of body
    posObs2cbHelioNorm = norm(posObs2cbHelio)
    bodyAngsize = math.degrees(
        math.atan2(radiusCB/2,posObs2cbHelioNorm))

    # take field of view angles from centerline
    horizFOVcenter = fov[0]/2.
    vertFOVcenter = fov[1]/2.

    # compute azimuth and elevation of celestial body
    az = math.degrees(math.atan2(eObs2cbBody[1], eObs2cbBody[0]))
    el = math.degrees(math.asin(eObs2cbBody[2]))

    # check if limits in camera specs
    chk = True
    if (abs(az) > horizFOVcenter + bodyAngsize) or (abs(el) > vertFOVcenter + bodyAngsize):
        chk = False

    if np.linalg.norm(posObs2cbBody) < radiusCB:
        chk = False

    return chk


#####################################
#####################################


def lumos(posCB, posObs, albedoCB, radCB, latRes, longRes):
    """Simulate illumination of a celestial body assuming constant geometric albedo on a spherical surface
    :param posCB: position of celestial body in heliocentric coordinates [m]
    :param posObs: position of observer in heliocentric coordinates [m]
    :param albedoCB: geometric albedo of celestial body (pure diffuse reflection)
    :param radCB: spherical radius of celestial body [km]
    :return:    array of flux reduction at each surface point (1 x N) [unitless],
                array of cartesian coordinates for surface mapping (N x 3) in heliocentric coordinates [m]
                square facetArea [m^2]    """

    from numpy import deg2rad, sin, cos, vstack, pi
    from datetime import datetime
    start = datetime.now()
    # Constants
    fluxRefDistance = 695700. # radius of the sun in km
    ########################   ############
    # Compute transformation matrix from body to heliocentric coordinates
    # celestial body coordinate frame: {+i towards sun, j (k cross i), +k upwards}

    # unit vector to sun from celestial body center in heliocentric coordinates
    distanceCB = np.linalg.norm(posCB)
    iBodyHelio = -posCB/ distanceCB

    # unit vector normal to i_body in the i-j body frame
    ijNormBodyHelio = np.array([-iBodyHelio[1], iBodyHelio[0], 0])
    ijNormBodyHelio = ijNormBodyHelio/np.linalg.norm(ijNormBodyHelio)
    kBodyHelio = np.cross(iBodyHelio,ijNormBodyHelio)
    kBodyHelio = kBodyHelio / np.linalg.norm(kBodyHelio)
    jBodyHelio = np.cross(kBodyHelio, iBodyHelio)
    jBodyHelio = jBodyHelio/ np.linalg.norm(jBodyHelio)

    # direction cosine matrix for body to helio coordinate frame
    dcmNB = (np.array([iBodyHelio,jBodyHelio,kBodyHelio])).T

    ####################################
    # Generate surface map of celestial body in body frame coordinates
    # create spherical grid using equal area facets
    ptsLatLong, facetArea = mapSphere(latRes, longRes, radCB)
    nLatLong = ptsLatLong.shape[0]


    # generate unit vector of surface normal vector for each surface point (in body frame)

    pts_eNormBody = np.zeros((nLatLong,3))

    latitudes = ptsLatLong[:,0]
    longitudes = ptsLatLong[:,1]

    eNorm_kBody = sin(deg2rad(latitudes))
    eNorm_ijNody = cos(deg2rad(latitudes))
    eNorm_iBody = sin(deg2rad(longitudes))*eNorm_ijNody
    eNorm_jBody = cos(deg2rad(longitudes))*eNorm_ijNody
    pts_eNormBody = vstack([eNorm_iBody,eNorm_jBody,eNorm_kBody]).T

    # generate surface map of ijk heliocentric coordinates from center of celestial body
    ptsHelio = np.zeros((nLatLong,3))

    zbody = radCB*sin(deg2rad(latitudes))
    xybody = radCB*cos(deg2rad(latitudes))
    ybody = xybody*cos(deg2rad(longitudes))
    xbody = xybody*sin(deg2rad(longitudes))
    xyzbody = vstack([xbody,ybody,zbody]).T
    ptsHelio = np.matmul(dcmNB,xyzbody.T).T
    indPts = 0

    ###################################################
    # Calculate net flux reduction due to diffuse reflection and sun to CB inverse square law

    # calculate geometric albedo using cosine law

    ptsAlbedo = np.zeros((nLatLong, 1))
    ptsAlbedo = np.zeros(nLatLong)
    ptsAlbedo[pts_eNormBody[:,0] > 0] = \
        pts_eNormBody[:,0][pts_eNormBody[:,0] > 0]
    ptsAlbedo = ptsAlbedo.reshape(len(ptsAlbedo),1)*albedoCB



    # calculate flux decay due to inverse square laws
    fluxDecaySun2cb = (fluxRefDistance/distanceCB)**2
    distanceCB2obs = np.linalg.norm(posCB-posObs)
    fluxDecayCB2obs = (radCB/distanceCB2obs)**2
    fluxDecayCB2obs = (1/distanceCB2obs)**2
    fluxDecayNet = fluxDecaySun2cb * ptsAlbedo * fluxDecayCB2obs
    return ptsHelio, fluxDecayNet, facetArea

###################################################
###################################################


def project2CamView(posCB, posObs, attdeCam, xyzHelio, fluxDecay, fov, facetArea):
    """Remove celestial body surface map points that are out of the camera view and convert heliocentric cartesian
    coordinates to azimuth and elevation from camera point of view
    :param posCB: position of celestial body in heliocentric coordinates [m]
    :param posObs: position of camera in heliocentric coordinates [m]
    :param dcmNB: direction cosine matrix of observer attitude [body to heliocentric coord. frame]
    :param field of view tuple (horizontal, vertical) [degrees]
    :param facetArea: area of single facet on spherical approximation [m^2]
    :return: array of azimuth elevation for each surface point (2 x N) [degrees]
            (azimuth is angle from i-k plane in i-j plane, elevation is angle from i-j plane)
            array of xyz for each visible surface point (helio coord. frame) [m]
            array of flux values at each surface point (1 x N) [W/m^2]
    """
    from numpy import sqrt, arccos, arcsin, arctan2, rad2deg, vstack, einsum
    from datetime import datetime
    import pdb

    start = datetime.now()
    # compute position of celestial body relative to observer in body coordinates
    dcmBN = attdeCam
    posObs2cbHelio = (posCB-posObs)
    eObs2cbHelio = posObs2cbHelio / np.linalg.norm(posObs2cbHelio)
    posObs2facetHelio = posObs2cbHelio + xyzHelio
    rOb2facetHelio =  sqrt(
        posObs2facetHelio[:,0]**2+posObs2facetHelio[:,1]**2+posObs2facetHelio[:,2]**2
        )

    rOb2facetHelio = rOb2facetHelio.reshape(len(rOb2facetHelio),1)
    eObs2facetHelio = posObs2facetHelio/rOb2facetHelio
    posObs2cbBody = np.matmul(dcmBN, posObs2cbHelio)
    # cycle through surface points and calculate dot product with unit vector from camera to CB
    # (in helio coordinates), input points with phase angle > 90 deg (visible to camera) in new arrays

    horizCamfov = fov[0] / 2.
    vertCamfov = fov[1] / 2.

    r = sqrt(xyzHelio[:,0]**2+xyzHelio[:,1]**2+xyzHelio[:,2]**2).reshape(
        len(xyzHelio),1)
    eXyzHelio = xyzHelio/r
    # cosCamPhase2 = np.dot(eObs2cbHelio, eXyzHelio.T)
    cosCamPhase = einsum('ij,ji->i',eObs2facetHelio,eXyzHelio.T)

    camPhase = arccos(cosCamPhase)

    # check if current facet is visible to observer
    ind = camPhase > math.pi/2
    eXyzHelio = eXyzHelio[ind]
    cosCamPhase = cosCamPhase[ind]
    camPhase = camPhase[ind]
    xyzHelio = xyzHelio[ind]
    xyzBody = np.matmul(dcmBN, xyzHelio.T)
    r_pt = posObs2cbHelio + xyzBody.T
    r_ptNorm = sqrt(r_pt[:,0]**2+r_pt[:,1]**2+r_pt[:,2]**2)
    az = rad2deg(arctan2(r_pt[:,1], r_pt[:,0]))
    el = rad2deg(arcsin(r_pt[:,2] / r_ptNorm))
    xyzCamview = xyzBody
    xyzCamviewHelio = xyzHelio
    azelPts = vstack((az,el)).T
    fluxDecayOut = fluxDecay[ind]
    xyzVisCam = xyzBody.T
    facetAreaCamview = facetArea * -cosCamPhase
    facetAreaCamview = facetAreaCamview.reshape(
        len(facetAreaCamview),1)

    xyzCamviewHelio = np.matmul(np.linalg.inv(dcmBN),xyzCamviewHelio.T).T
    xyzVisCam = np.matmul(np.linalg.inv(dcmBN),xyzVisCam.T).T


    return azelPts, xyzCamviewHelio, xyzVisCam, fluxDecayOut, facetAreaCamview



###################################################
###################################################

def xyz2RADec(xyzCB):
    """Calculate right ascension and declination from xyz coordinates
        :param xyzCB: 1x3 array of position
        :return raDec: 1x2 array of center of CB in right ascension and declination
                ra [deg] is angle from [1,0,0] to point in the xy horizontal plane
                dec [deg] is angle above xy plane to point"""

    from numpy import arctan2, arcsin, sqrt, rad2deg
    from datetime import datetime
    start = datetime.now()
    if len(xyzCB.shape) == 1:

        r = np.linalg.norm(xyzCB)
        ra = math.degrees(math.atan2(xyzCB[1], xyzCB[0]))
        dec = math.degrees(math.asin(xyzCB[2]/r))

    else:


        r = sqrt(xyzCB[:,0]**2+xyzCB[:,1]**2+xyzCB[:,2]**2)
        ra = rad2deg(arctan2(xyzCB[:,1],xyzCB[:,0]))
        dec = rad2deg(arcsin(xyzCB[:,2]/r))
        ra = ra.reshape(len(ra),1)
        dec = dec.reshape(len(dec),1)

    raDec = (ra, dec)
    return raDec


###################################################
###################################################



def lightSim(attdeCam, posCam, posCB, fov, latRes, longRes, doPtSource,
    albedoCB, radiusCB, nameCB):
    """Conduct lighting simulation of multiple bodies.
    :param posCB: Nx3 array of positions of celestial bodies in heliocentric coordinates [m]
    :param posObs: position of camera in heliocentric coordinates [m]
    :param attdeCam: direction cosine matrix of observer attitude [body to heliocentric coord. frame]
    :param field of view tuple (horizontal, vertical) [degrees]
    :param doPtSource: bool whether to treat CB as a point source
    :return: dict of summary output for each visible CB with nameCB as key
             summary output dict has the following values
            {'bodypos': (npts x 3) numpy array [m]
            'facetRA' : right ascension of facet in heliocentric coord. frame [deg]
            'facetDec': declination of facet in heliocentric coord. frame [deg]
            'fluxDecay': flux decay due to inverse square laws and diffuse reflection
            'facetArea': projected facet area normal to observer [m^2]
    """
    from numpy import array

    # check if CB is in field of view
    fovChk = checkFoV(posCB, posCam, attdeCam, fov, radiusCB)
    currentCbDict = -1

    # compute lighting simulation
    if fovChk:

        # illuminate sphere
        ptsPosHelio, ptsAlbedoHelio, facetArea = lumos(
            posCB, posCam, albedoCB, radiusCB, latRes, longRes)

        # limit points to those in camera field of view and compute azel
        ptsAzelCam, ptsXyzHelio, ptsXyzCam, ptsAlbedoCam, facetAreaCam = project2CamView(
            posCB, posCam, attdeCam, ptsPosHelio, ptsAlbedoHelio, fov, facetArea)

        # check to see if there are any illuminated points visible for the celestial body
        # (celestial body may be in field of view with no illuminated points (observer on dark side of body)
        if ptsAlbedoCam.shape[0] > 0:
            currentCbDict = {}
            if doPtSource:

                # compute position of CB center in camera view
                dcmBN = attdeCam
                posObs2cbHelio = (posCB - posCam)
                posObs2cbBody = np.matmul(dcmBN, posObs2cbHelio.T)

                raDec = xyz2RADec(posObs2cbBody)

                #bodypos needs to be the positon of this facet relative
                #to the body. Because a point source facet will always
                #be in the same place as the center of its body,
                #currentCbDict['bodypos'] should always be zero for
                #a point source 
                currentCbDict['bodypos'] = array([[0,0,0]])
                currentCbDict['facetRA'] = array([raDec[0]])
                currentCbDict['facetDec'] = array([raDec[1]])
                currentCbDict['netAlbedo'] = ptsAlbedoCam
                currentCbDict['facetArea'] = facetAreaCam

            else:

                raDec = xyz2RADec(ptsXyzHelio)

                currentCbDict['bodypos'] = ptsXyzCam
                currentCbDict['facetRA'] = raDec[0]
                currentCbDict['facetDec'] = raDec[1]
                currentCbDict['netAlbedo'] = ptsAlbedoCam
                currentCbDict['facetArea'] = facetAreaCam

    return currentCbDict

