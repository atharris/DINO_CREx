import math
import matplotlib.pyplot as plt
import numpy as np
# from scipy import misc

import imageProcessingFunctions as imfunc
import searchLocationFunctions as locfunc
import objectIDFunctions as idfunc
import dynamicFunctions as dyn


def imageProcessing(imageMap, cameraParameters, r_N_cam, sigma_BN_est,
                    r_N_beacons, beaconIDs, beaconRadius, makePlots=False, debugMode=False):
    """
    Identifies objects in an image using a reference catalog.
    @param  imageMap            NxM numpy array of image where N is image horizontal resolution and M is vertical
                                resolution.
    @param  cameraParameters:    python dict with the following entries
                                ['resolution'] tuple of horizontal x vertical camera sensor resolution
                                ['focal length'] camera sensor effective focal length [m]
                                ['pixel size'] tuple of horizontal x vertical camera sensor pixel size [m]
                                ['field of view'] tuple of horizontal x vertical field of view [deg]
    @param  r_N_cam
    @param  sigma_BN_est        Initial estimate of s/c attitude in a direction cosine matrix
    @param  beaconIDs
    @param  beaconRadius
    @return objectID            List of reference catalog ID values for identified objects (NaN if none found)
    @return pixelLine           Tuple of pixel, line coordinates
    @return sigma_BN            S/C attitude in modified rodriguez parameters
    """


    # debugMode = True

    ##################################################
    ##################################################
    # Parameter Definition

    pixelSizeMin = 15.  # number of pixels required for a beacon to be considered a resolved body
    pixelDimenMin = min(cameraParameters['pixel size'])
    singlePixelAngSize = math.degrees(math.atan(pixelDimenMin / cameraParameters['focal length']))
    angularSizeMin = singlePixelAngSize * pixelSizeMin

    minRes = min(cameraParameters['resolution'])
    imageMax = np.amax(imageMap)
    imageMax = 2**32
    ROI_parameters = {}
    ROI_parameters['signal_threshold'] = 200.
    ROI_parameters['ROI_border_width'] = 1
    ROI_parameters['max_search_dist'] = minRes/10.

    imageProcessingParam ={}
    imageProcessingParam['voteCountMinRatio'] = .25      # minimum ratio of positive matches out of possible matches
    imageProcessingParam['dthetaMax'] = 12.             #[deg] dependent on object ID reference catalog

    # filepath to catalog files relative to image processing unit test locations
    # imageProcessingParam['filenameSearchCatalog'] = '../../../../external/tycho_mag_cutoff.db'
    # imageProcessingParam['filenameObjectIDCatalog'] = '../../../../external/objectID_catalog.db'

    # filepath to catalog files relative to dinoScenarios
    imageProcessingParam['filenameSearchCatalog'] = '../external/tycho_mag_cutoff.db'
    imageProcessingParam['filenameObjectIDCatalog'] = '../external/objectID_catalog.db'
    imageProcessingParam['dthetaError'] = 5E-4
    maxInitialEstimates = 10

    if debugMode:
        print '\nParameters:'
        print 'Signal Threshold: ', ROI_parameters['signal_threshold']
        print 'Max Search Distance: ', ROI_parameters['max_search_dist']
        print 'Pixel Size Min: ', pixelSizeMin

    ##################################################
    ##################################################
    # Process Inputs

    # normalize to a 12 bit image
    # imageMap = (imageMap/np.amax(imageMap)) * 4096

    # convert modified rodriguez parameter to a direction cosine matrix
    BN_dcm_cam = dyn.mrp2dcm(sigma_BN_est)
    # BN_dcm_cam = dcm_BN_est


    ##################################################
    ##################################################
    # Check Current status of target beacons
    # 0 for out of field of view, 1 for point source, 2 for resolved body

    numBeacons = len(beaconIDs)
    beaconStatus = imfunc.checkBeaconStatus(r_N_beacons, r_N_cam, BN_dcm_cam,
                                            cameraParameters['field of view'],
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

    if debugMode:
        print '\nBeacon Status'
        print 'Visible Beacons: ', beaconIDsVisible
        print 'Point Source Beacons: ', beaconIDsPtSource
        print 'Resolved Body Beacons: ', beaconIDsResolved


    ##################################################
    ##################################################
    # Conduct initial search location for beacons and stars in field of view
    # Note that the order is important as earlier estimates are detected in the image first
    # Order is given in priority to resolved beacons, pt source beacons, then stars


    # search for stars in field of view
    pixel_line_stars_i, catalogIDs = locfunc.initial_stars_estimate(
        BN_dcm_cam, cameraParameters, imageProcessingParam['filenameSearchCatalog'])

    pixel_line_ptSource_i = pixel_line_stars_i

    if len(pixel_line_ptSource_i) > maxInitialEstimates:
        pixel_line_ptSource_i = pixel_line_ptSource_i[0:maxInitialEstimates]

    if numBeaconsPtSource > 0:
        pixel_line_beacon_ptSource_i = locfunc.initial_beacons_estimate(
            N_r_beaconsPtSource, r_N_cam, BN_dcm_cam, cameraParameters)
        pixel_line_ptSource_i = np.vstack((pixel_line_beacon_ptSource_i,
                                           pixel_line_ptSource_i))

    pixelLineInitialEstimates = pixel_line_ptSource_i


    if numBeaconsResolved > 0:
        pixel_line_resolved_i = locfunc.initial_beacons_estimate(
            N_r_beaconsResolved, r_N_cam, BN_dcm_cam, cameraParameters)
        pixelLineInitialEstimates = np.vstack((pixel_line_resolved_i,
                                               pixelLineInitialEstimates))


    if debugMode:
        print '\n# of Initial Estimates of Pt Sources: ', len(pixel_line_ptSource_i)
        for currentPL in pixel_line_ptSource_i:
            print round(currentPL[0], 4), round(currentPL[1],4)

        if numBeaconsResolved > 0:
            print '\nInitial Estimate of Resolved Bodies'
            for currentPL in pixel_line_resolved_i:
                print currentPL


    ##################################################
    ##################################################
    # Determine Center and Centroid Locations

    pixelLineCenterFound = []           # all objects with centers detected
    beaconIDsFound = []
    pixelLineCenterResolvedFound = []
    pixelLineCenterBeaconFound = []
    pixelLineCenterStars = []
    refCatalogIDsMatched = []


    # find centers of resolved beacons
    if numBeaconsResolved > 0:
        pixelLineCenterResolved, DN = imfunc.findCenterResolvedBody(
            imageMap, pixel_line_resolved_i, ROI_parameters)

        for indBeacon in range(len(pixelLineCenterResolved)):
            currentPL = pixelLineCenterResolved[indBeacon]

            if currentPL is not None:
                beaconIDsFound.append(beaconIDsResolved[indBeacon])
                pixelLineCenterResolvedFound.append(currentPL)
                pixelLineCenterBeaconFound.append(currentPL)
                # pixelLineCenterFound = np.vstack((currentPL, pixelLineCenterFound))

        numBeaconsResolvedFound = len(beaconIDsFound)

    else:
        numBeaconsResolvedFound = 0


    # find centers of stars and pt source beacons
    pixelLineCenter, DN = imfunc.findCentroidPointSource(
        imageMap, pixel_line_ptSource_i, ROI_parameters)

    # separate pt source beacons from stars
    for indCentroid in range(len(pixelLineCenter)):
        currentPL = pixelLineCenter[indCentroid]

        if currentPL is not None:

            if indCentroid < numBeaconsPtSource:
                pixelLineCenterFound.append(currentPL)

            if indCentroid >= numBeaconsPtSource:
                pixelLineCenterStars.append(currentPL)
                refCatalogIDsMatched.append(catalogIDs[indCentroid])


            if indCentroid < numBeaconsPtSource:
                beaconIDsFound.append(beaconIDsPtSource[indCentroid])
                pixelLineCenterBeaconFound.append(currentPL)


    numStarsFound = len(pixelLineCenterStars)
    numBeaconsFound = len(pixelLineCenterBeaconFound)
    numObjectsFound = len(pixelLineCenterFound)

    if debugMode:

        print 'Centroiding and Center-finding:'
        # print '\n# of Measured Pixel and Line Center Locations: ', numObjectsFound
        print '# of Measured Pt Source Locations: ', numStarsFound
        if numBeaconsResolvedFound > 0:
            print '# of Measured Resolved Beacons: ', numBeaconsResolvedFound
        print '\nMeasured Beacon Locations'
        for ind in range(numBeaconsFound):
            print round(pixelLineCenterBeaconFound[ind][0], 4), round(pixelLineCenterBeaconFound[ind][1],4)
        print '\nMeasured Star Locations'
        for ind in range(numStarsFound):
            print round(pixelLineCenterStars[ind][0], 4), round(pixelLineCenterStars[ind][1],4)



    ##################################################
    ##################################################
    # Conduct Object ID for Stars

    # if numBeaconsResolved == 0:
    if numStarsFound > 3:
        objectIDs = idfunc.objectIDStars(pixelLineCenterStars, imageProcessingParam, cameraParameters)
        numObjectsID = len(objectIDs)

        # remove unidentified object IDs from results
        objectIDsMatched = []
        pixelLineCenterFoundMatched = []


        print numObjectsID
        print pixelLineCenterFound
        for ind in range(numObjectsID):
            if objectIDs[ind] is not None:
                objectIDsMatched.append(objectIDs[ind])
                pixelLineCenterFoundMatched.append(pixelLineCenterStars[ind])

        numObjectsFoundMatched = len(objectIDsMatched)

        if debugMode:
            print '\nObject ID ...'
            print 'ObjectID Output: '
            print objectIDsMatched

    else:
        numObjectsFoundMatched = 0


    ##################################################
    ##################################################
    # Attitude Determination

    print '\nAttitude Determination ...'

    if numObjectsFoundMatched >= 3:

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

        if debugMode:

            print '\nStar #', '\t\t', \
                "%6s" % 'Pixel ', '\t\t', \
                "%6s" % 'Line  ', '\t\t', \
                "%6s" % 'Ref ID', '\t\t', \
                "%6s" % 'Obj ID', '\t\t', \
                "%6s" % 'RA Dec (ref for Obj ID)'

            for ind in range(numObjectsFoundMatched):
                print 'Star #', ind, '\t:', \
                    "%6s" % round(pixelLineCenterStars[ind][0], 2), '\t', \
                    "%6s" % round(pixelLineCenterStars[ind][1], 2), '\t\t', \
                    "%6s" % refCatalogIDsMatched[ind], '\t\t', \
                    "%6s" % objectIDsMatched[ind], \
                    "%6s" % radec[ind][0], \
                    "%6s" % radec[ind][1]

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

    if makePlots:

        # original image
        fig1 = plt.figure(1)
        plt.imshow(imageMap, interpolation = 'none', cmap='viridis')
        plt.title('Input Image')
        plt.savefig('input_image.png')

        # initial estimates
        fig2 = plt.figure(2)
        plt.imshow(imageMap, interpolation = 'none', cmap='viridis')
        for ind in range(len(pixelLineInitialEstimates)):
            plt.scatter(pixelLineInitialEstimates[ind][0]-.5,
                        pixelLineInitialEstimates[ind][1]-.5,
                        color='r', marker='+', s=10)
        plt.title('Initial Estimates')
        plt.savefig('initial_estimates.png')

        # measured center locations
        fig3 = plt.figure(3)
        plt.imshow(imageMap, interpolation='none', cmap='viridis')
        for ind in range(numObjectsFound):
            # plt.scatter(pixelLineInitialEstimates[ind][0], pixelLineInitialEstimates[ind][1],
            #             color='r', marker='+', s=10)
            plt.scatter(pixelLineCenterFound[ind][0]-.5,
                        pixelLineCenterFound[ind][1]-.5,
                        color='y', marker='+', s=30)
        plt.title('Measured Center Locations')
        plt.savefig('measured_centers.png')
        # plt.show()
        plt.close()

        # numZoomPlots = 4
        # windowSize = 10
        # for ind in range(numZoomPlots):
        #     fig4 = plt.figure(4)
        #     x = int(round(pixelLineInitialEstimates[ind][0]))
        #     y = int(round(pixelLineInitialEstimates[ind][1]))
        #     print '\n', x ,y, pixelLineCenterFound[ind], pixelLineInitialEstimates[ind]
        #     windowXMin = x - windowSize
        #     windowXMax = x + windowSize
        #     windowYMin = y - windowSize
        #     windowYMax = y + windowSize
            # plt.imshow(imageMap[windowYMin: windowYMax,
            #            windowXMin:windowXMax],
            #            interpolation='none', cmap='viridis')

            # add half pixel to initial estimates to compensate for imshow plotting
            # plt.scatter(pixelLineInitialEstimates[ind][0]-windowXMin-.5,
            #             pixelLineInitialEstimates[ind][1]-windowYMin-.5,
            #             color='g', marker='+', s=50)
            # plt.scatter(pixelLineCenterFound[ind][0]-windowXMin-.5,
            #             pixelLineCenterFound[ind][1]-windowYMin-.5,
            #             color='r', marker='+', s=50)
            # plt.scatter(pixelLineCenterFound[ind][0]-windowXMin+.5,
                        # pixelLineCenterFound[ind][1]-windowYMin+.5,
                        # color='g', marker='+', s=30)
            # plt.scatter(0,0,color='w',marker='+',s=30)
            # plt.show()


    ##################################################
    ##################################################
    # Set Outputs

    idOut = beaconIDsFound
    pixelLineBeaconOut = pixelLineCenterBeaconFound
    sigmaOut = mrpQUEST

    return idOut, pixelLineBeaconOut, sigmaOut
