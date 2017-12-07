import math
import matplotlib.pyplot as plt
import numpy as np
# from scipy import misc

import imageProcessingFunctions as imfunc
import searchLocationFunctions as locfunc
import objectIDFunctions as idfunc
import dynamicFunctions as dyn



def imageProcessing(imageMap, cameraParameters, r_N_cam, sigma_BN_est,
                    r_N_beacons, beaconIDs, beaconRadius, makePlts = False):
    """
    Identifies objects in an image using a reference catalog.
    @param  imageMap            1xN numpy array of image where N is image horizontal resolution x image vertical
                                resolution.
    @param  cameraParameters:    python dict with the following entries
                                ['resolution'] tuple of horizontal x vertical camera sensor resolution
                                ['focal length'] camera sensor effective focal length [m]
                                ['pixel size'] tuple of horizontal x vertical camera sensor pixel size [m]
                                ['field of view'] tuple of horizontal x vertical field of view [deg]
    @param  r_N_cam
    @param  sigma_BN_est        Initial estimate of s/c attitude in modified rodriguez parameters
    @param  beaconIDs
    @param  beaconRadius
    @return objectID            List of reference catalog ID values for identified objects (NaN if none found)
    @return pixelLine           Tuple of pixel, line coordinates
    @return sigma_BN            S/C attitude in modified rodriguez parameters
    """

    ##################################################
    ##################################################
    # Parameter Definition

    pixelSizeMin = 15.  # number of pixels required for a beacon to be considered a resolved body
    pixelDimenMin = min(cameraParameters['pixel size'])
    singlePixelAngSize = math.degrees(math.atan(pixelDimenMin / cameraParameters['focal length']))
    angularSizeMin = singlePixelAngSize * pixelSizeMin

    minRes = min(cameraParameters['resolution'])
    ROI_parameters = {}
    ROI_parameters['signal_threshold'] = 4096/1000.
    ROI_parameters['ROI_border_width'] = 1
    ROI_parameters['max_search_dist'] = minRes/100.

    imageProcessingParam ={}
    imageProcessingParam['voteCountMinRatio'] = .2     # minimum ratio of positive matches out of possible matches
    imageProcessingParam['dthetaMax'] = 12.     #[deg] dependent on object ID reference catalog
    imageProcessingParam['filenameSearchCatalog'] = '../../../../external/tycho_mag_cutoff.db'
    imageProcessingParam['filenameObjectIDCatalog'] = '../../../../external/objectID_catalog.db'
    imageProcessingParam['dthetaError'] = 1E-4

    maxInitialEstimates = 20


    ##################################################
    ##################################################
    # Process Inputs

    # normalize to a 12 bit image and reshape 1D pixel map to 2D pixel map
    imageMap = (imageMap/np.amax(imageMap)) * 4096
    # imageMap = imageMap.reshape(cameraParameters['resolution'][0], cameraParameters['resolution'][1])

    # convert modified rodriguez parameter to a direction cosine matrix
    BN_dcm_cam = dyn.mrp2dcm(sigma_BN_est)


    ##################################################
    ##################################################
    # Check Current status of target beacons
    # 0 for out of field of view, 1 for point source, 2 for resolved body

    BN_dcm_cam = dyn.mrp2dcm(sigma_BN_est)

    numBeacons = len(beaconIDs)
    beaconStatus = imfunc.checkBeaconStatus(r_N_beacons, r_N_cam, BN_dcm_cam, cameraParameters['field of view'],
                                            beaconRadius, angularSizeMin)

    # Separate visible beacons from non-visible
    N_r_beaconsVisible = []
    N_r_beaconsPtSource = []
    N_r_beaconsResolved = []
    beaconIDsVisible = []
    beaconIDsPtSource = []
    beaconIDsResolved = []

    for ind in range(numBeacons):

        if beaconStatus[ind] != 0:
            N_r_beaconsVisible.append(r_N_beacons[ind])
            beaconIDsVisible.append(beaconIDs[ind])

        if beaconStatus[ind] == 1:
            N_r_beaconsPtSource.append(r_N_beacons[ind])
            beaconIDsPtSource.append(beaconIDs[ind])

        if beaconStatus[ind] == 2:
            N_r_beaconsResolved.append(r_N_beacons[ind])
            beaconIDsResolved.append(beaconIDs[ind])

    numBeaconsVisible = len(N_r_beaconsVisible)
    numBeaconsPtSource = len(N_r_beaconsPtSource)
    numBeaconsResolved = len(N_r_beaconsResolved)

    print '\nBeacon Status'
    print 'Visible Beacons: ', beaconIDsVisible
    print 'Point Source Beacons: ', beaconIDsPtSource
    print 'Resolved Body Beacons: ', beaconIDsResolved


    ##################################################
    ##################################################
    # Conduct initial search location for beacons and stars in field of view
    print '\nInitial Search Location ...'

    # search for stars in field of view
    pixel_line_stars_i, catalogIDs = locfunc.initial_stars_estimate(
        BN_dcm_cam, cameraParameters, imageProcessingParam['filenameSearchCatalog'])

    pixel_line_ptSource_i = pixel_line_stars_i

    if len(pixel_line_ptSource_i) > 20:
        pixel_line_ptSource_i = pixel_line_ptSource_i[0:maxInitialEstimates]


    if numBeaconsPtSource > 0:

        pixel_line_beacon_ptSource_i = locfunc.initial_beacons_estimate(
            N_r_beaconsPtSource, r_N_cam, BN_dcm_cam, cameraParameters)

        pixel_line_ptSource_i = np.vstack((pixel_line_ptSource_i,
                                           pixel_line_beacon_ptSource_i))

    pixelLineInitialEstimates = pixel_line_ptSource_i


    if numBeaconsResolved > 0:

        pixel_line_resolved_i = locfunc.initial_beacons_estimate(
            N_r_beaconsResolved, r_N_cam, BN_dcm_cam, cameraParameters)

        pixelLineInitialEstimates = np.vstack((pixel_line_resolved_i,
                                               pixelLineInitialEstimates))


    print '\n# of Initial Estimates of Pt Sources: ', len(pixel_line_ptSource_i)
    # for currentPL in pixel_line_ptSource_i:
    #     print currentPL

    if numBeaconsResolved > 0:
        print '\nInitial Estimate of Resolved Bodies'

        for currentPL in pixel_line_resolved_i:
            print currentPL


    ##################################################
    ##################################################
    # Determine Center and Centroid Locations
    print '\nCentroiding Center-finding ...'

    pixelLineCenter, DN = imfunc.findCentroidPointSource(
        imageMap, pixel_line_ptSource_i, ROI_parameters)

    if numBeaconsResolved > 0:

        pixelLineCenterResolved, DN = imfunc.findCenterResolvedBody(
            imageMap, pixel_line_resolved_i, ROI_parameters)

        pixelLineCenter = np.vstack(pixelLineCenterResolved, pixelLineCenter)

    pixelLineCenterFound = []
    for currentPL in pixelLineCenter:
        if currentPL is not None:
            pixelLineCenterFound.append(currentPL)

    numObjectsFound = len(pixelLineCenterFound)

    print '\n# of Identified Pixel and Line Center Locations: ', numObjectsFound
    for ind in range(numObjectsFound):
        print pixelLineCenterFound[ind]


    ##################################################
    ##################################################
    # Conduct Object ID
    print '\nObject ID ...'

    beaconIDsFound = []
    pixelLineBeaconFound = []

    if numBeaconsResolved > 0:
        for ind in range(numBeacons):
            if beaconStatus[ind] == 2:
                beaconIDsFound.append(beaconIDs)
                pixelLineBeaconFound.append(pixelLineCenterResolved[ind])

    if numBeaconsVisible == 0:
        objectIDs = idfunc.objectIDStars(pixelLineCenterFound, imageProcessingParam, cameraParameters)

        print '\nStar #', '\t\t', \
            "%6s" % 'Pixel ', '\t\t', \
            "%6s" % 'Line  ', '\t\t', \
            "%6s" % 'Ref ID', '\t\t', \
            "%6s" % 'Obj ID', '\t\t'
        for ind in range(numObjectsFound):
            print 'Star #', ind, '\t:', \
                "%6s" % round(pixelLineCenter[ind][0], 2), '\t', \
                "%6s" % round(pixelLineCenter[ind][1], 2), '\t\t', \
                "%6s" % catalogIDs[ind], '\t\t', \
                "%6s" % objectIDs[ind]

        # remove unidentified object IDs from results
        objectIDsMatched = []
        pixelLineCenterFoundMatched = []

        for ind in range(numObjectsFound):
            if objectIDs[ind] is not None:
                objectIDsMatched.append(objectIDs[ind])
                pixelLineCenterFoundMatched.append(pixelLineCenterFound[ind])

        numObjectsFoundMatched = len(objectIDsMatched)

    ##################################################
    ##################################################
    # Attitude Determination

    print '\nAttitude Determination ...'

    if numBeaconsResolved == 0 and numObjectsFoundMatched >= 2:

        radec = idfunc.catalogIDsToRaDec(objectIDsMatched,
                                         imageProcessingParam['filenameSearchCatalog'])

        # List containers for unit vectors (1x3) numpy arrays
        B_ehat = []
        N_ehat = []
        weight = []

        for indStar in range(numObjectsFoundMatched):

            B_ehatTemp =idfunc.pixelline_to_ehat(
                pixelLineCenterFoundMatched[indStar],
                cameraParameters)

            B_ehat.append(B_ehatTemp)

            N_ehatTemp = dyn.radec_to_unitvector(radec[indStar])
            N_ehat.append(N_ehatTemp)

            # verify results with actual s/c attitude
            B_ehatCheck = np.matmul(BN_dcm_cam, N_ehatTemp)

            weight.append(1)

            # print B_ehatCheck - B_ehatTemp

        dcmQUEST = dyn.quest(B_ehat, N_ehat)
        mrpQUEST = dyn.dcm2mrp(dcmQUEST)

        print '\nQUEST Attitude Solution: '
        print dcmQUEST

        print '\nS/C Attitude Truth: '
        print BN_dcm_cam

        print '\nMRP Solution: ', mrpQUEST

    else:
        mrpQUEST = None


    ##################################################
    ##################################################
    # Optional Plot Generation

    if makePlts:

        # original image
        fig1 = plt.figure(1)
        plt.imshow(imageMap, interpolation = 'none', cmap='viridis')
        plt.title('Input Image')
        plt.savefig('input_image.png')

        # initial estimates
        fig2 = plt.figure(2)
        plt.imshow(imageMap, interpolation = 'none', cmap='viridis')
        for ind in range(len(pixelLineInitialEstimates)):
            plt.scatter(pixelLineInitialEstimates[ind][0], pixelLineInitialEstimates[ind][1],
                        color='r', marker='+', s=10)
        plt.title('Initial Estimates')
        plt.savefig('initial_estimates.png')

        # measured center locations
        fig3 = plt.figure(3)
        plt.imshow(imageMap, interpolation='none', cmap='viridis')
        for ind in range(numObjectsFound):
            # plt.scatter(pixelLineInitialEstimates[ind][0], pixelLineInitialEstimates[ind][1],
            #             color='r', marker='+', s=10)
            plt.scatter(pixelLineCenter[ind][0], pixelLineCenter[ind][1],
                        color='y', marker='+', s=10)
        plt.title('Measured Center Locations')
        plt.savefig('measured_centers.png')
        # plt.show()


    ##################################################
    ##################################################
    # Set Outputs

    idOut = objectIDs
    pixelLineBeaconOut = pixelLineCenterFound

    if numBeaconsVisible == 0:
        idOut = None
        pixelLineBeaconOut = None

    if numBeaconsVisible > 0:
        idOut = beaconIDsFound
        pixelLineBeaconOut = pixelLineBeaconFound

    sigmaOut = mrpQUEST

    return idOut, pixelLineBeaconOut, sigmaOut
