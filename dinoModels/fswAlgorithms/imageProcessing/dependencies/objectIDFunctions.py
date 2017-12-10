# Title   : find_centroid_functions.py
# Author  : Joe Park
# Date    : 08/22/2017
# Synopsis: Finds centroid of a point source from a pixel map image. DINO C-REx module.
# Dependencies:     Specified star catalog file in '/star_catalog' directory

import math         # common math functions
import numpy as np  # matrix algebra
import sqlite3 as sql
import dynamicFunctions as dyn


##################################################
##################################################

def pixelline_to_radec(pl_in, attde_sc, cam_param):
    # Convert pixel / line coordinates to inertial RA and Dec
    # P/L --> body unit vector --> inertial unit vector --> inertial RA and Dec


    cam_focal_length = cam_param['focal length']
    cam_res = cam_param['resolution']
    cam_pixel_size = cam_param['pixel size']

    # convert pixel/line coordinates to a unit vector in camera body coordinates
    ehat_cam = pixelline_to_ehat(pl_in, cam_param)

    # convert to inertial coordinates using DCM
    ehat_helio = np.matmul(attde_sc.T, ehat_cam)

    # connvert inertial unit vector to RA and Dec
    radec = dyn.unitvector_to_radec(ehat_helio)

    return radec


##################################################
##################################################

# Does not account for image inversion through focal point
# Output:       ehat        unit vector in camera body coordinates

def pixelline_to_ehat(pl_in, cam_param):

    cam_focal_length = cam_param['focal length']
    cam_res = cam_param['resolution']
    cam_pixel_size = cam_param['pixel size']

    pixel = pl_in[0]
    line = pl_in[1]

    x_upperleft = cam_res[0] / 2. * cam_pixel_size[0]
    y_upperleft = cam_res[1] / 2. * cam_pixel_size[1]

    # distance to point on focal plane from upper left corner
    x_focal = pixel * cam_pixel_size[0]
    y_focal = line * cam_pixel_size[1]

    # calculate vector from focal plane to focal point in camera body coordinates
    vector_cam = np.array([cam_focal_length,
                           -(x_focal - x_upperleft),
                           -(y_focal - y_upperleft)])
    ehat_cam = vector_cam/np.linalg.norm(vector_cam)

    return ehat_cam


##################################################
##################################################

# Input:    pl_pair             ((pixel1, line1), (pixel2, line2))
#           cam_res
# Output    angular distance

def pixelline_to_angular_separation(pl1, pl2, cam_param):

    # Find unit vector in camera coordinates for each pixel/line coordinate
    ehat1_cam = pixelline_to_ehat(pl1, cam_param)
    ehat2_cam = pixelline_to_ehat(pl2, cam_param)

    # Find angular separation between two unit vectors
    ang_separation = math.degrees(math.acos(np.dot(ehat1_cam, ehat2_cam)))

    return ang_separation


##################################################
##################################################

def observed_measurements(pl_in, dthetaMax, cam_param):
    """
    Generates necessary measurements for Object ID algorithm from pixel/line coordinates of visible objects in image.

    @param pl_in The pixel/line coordinate of the object in the image //tuple of 1-D numpy arrays
    @param dthetaMax
    @param cam_param
    @return dtheta
    """

    dtheta = []
    pairIndex = []

    for ind in range(len(pl_in)):

        ind2 = ind + 1

        while ind2 < len(pl_in):

            current_dtheta = pixelline_to_angular_separation(
                pl_in[ind], pl_in[ind2], cam_param)

            if current_dtheta <= dthetaMax:
                dtheta.append(current_dtheta)
                pairIndex.append((ind, ind2))

            ind2 = ind2+1

    # print '\nObserved Delta Theta Calculations'
    # for ind in range(len(dtheta)):
    #     print pairIndex[ind], dtheta[ind]

    return dtheta, pairIndex


##################################################
##################################################

def generate_uniqueID_and_counts(rawIDs):
    """
    Generates list of unique IDs and counts of each.
    """

    uniqueIDs = []
    counts = []

    for currentID in rawIDs:

        if currentID not in uniqueIDs:

            uniqueIDs.append(currentID)
            netCount = 0

            for currentID2 in rawIDs:

                if currentID2 == currentID:
                    netCount = netCount+1

            counts.append(netCount)

    return uniqueIDs, counts

##################################################
##################################################

def objectIDStars(pl_in, imageProcessingParam, cameraParameters):
    """
    Identifies objects in an image using a reference catalog.
    @param  pl_in               Tuple of 1xN arrays of right ascension and declination in HCI coordinates [deg]
                                Assumed to be center locations of objects in image.
    @param  attde_sc            Camera attitude BN direction cosine matrix (camera to heliocentric coord. frame)
    @param imageProcessingParam Dict of user specified parameters
                                'filenameSearchCatalog' : SQL database filename in XXXXX catalog directory
                                'filenameObjectIDCatalog' : SQL database filename in XXXX catalog directory
                                'dthetaMax' : max value of dtheta in object ID reference catalog [deg]
                                'dthetaError' : one-sided error bound for ref. dtheta values in object ID catalog [deg]
                                'voteCountMinRatio' : ratio of total votes required out of total possible to considered
                                a positive ID match.
    @param cameraParameters:    python dict with the following entries
                                ['resolution'] tuple of horizontal x vertical camera sensor resolution
                                ['focal length'] camera sensor effective focal length [m]
                                ['pixel size'] tuple of horizontal x vertical camera sensor pixel size [m]
    @return objectID            List of reference catalog ID values for identified objects (None if none found)
    @return observationWeight   List of observation weight based on vote count (None if none found)
    """

    voteMinRatio = imageProcessingParam['voteCountMinRatio']
    dthetaMax = imageProcessingParam['dthetaMax']
    fname_catalog =imageProcessingParam['filenameObjectIDCatalog']
    ERROR_DTHETA = imageProcessingParam['dthetaError']

    nStars = len(pl_in)

    # container for each star to append running vote counts
    # value of -1 indicates empty (no matches)
    runningVoteCount = {}
    for indStar in range(nStars):
        runningVoteCount.update({indStar: -1})

    # container for final votes
    starID = []

    # calculate intra-stellar angular distance between each pair in star coordinate list
    dtheta, pairIndex = observed_measurements(pl_in, dthetaMax, cameraParameters)
    numSearches = len(dtheta)

    # determine number of object ID catalog searches for each object
    searchesPerObject = np.zeros(nStars)
    for indStar in range(nStars):
        for indSearch in range(numSearches):
            if indStar == pairIndex[indSearch][0] or indStar == pairIndex[indSearch][1]:
                searchesPerObject[indStar] = searchesPerObject[indStar]+1

    print '\nNumber of Object ID Catalog Searches'
    for ind in range(nStars):
        print searchesPerObject[ind]

    # open object ID reference catalog
    ref_catalog = sql.connect(fname_catalog)
    with ref_catalog:

        cursor = ref_catalog.cursor()

        # loop through each angular pair and find matching entries
        # loops through all (0,N) values then (1,N) values , etc...
        print '\nFinding Matching dtheta entries in catalog'
        print 'Total Number of catalog searches: ', numSearches

        for ind in range(len(dtheta)):

            # list of matched ID's for current dtheta measurement
            matchedIDs = []

            star1 = pairIndex[ind][0]
            star2 = pairIndex[ind][1]

            cursor.execute("SELECT * FROM tycho_objectID WHERE dtheta BETWEEN (?) AND (?)",
                           (dtheta[ind]-ERROR_DTHETA, dtheta[ind]+ERROR_DTHETA))
            rows = cursor.fetchall()

            # add matched ID's to a running list
            if len(rows) > 0:

                for row in rows:
                    matchedIDs.append(row[1])
                    matchedIDs.append(row[2])

                # store running list of matched ID's to container for each star used in dtheta calculation
                if runningVoteCount[star1] == -1:
                    runningVoteCount[star1] = list(matchedIDs)
                else:
                    runningVoteCount[star1].extend(matchedIDs)

                if runningVoteCount[star2] == -1:
                    runningVoteCount[star2] = list(matchedIDs)
                else:
                    runningVoteCount[star2].extend(matchedIDs)

        print 'Number of Catalog Searches: ', len(dtheta)

    # transform running list of IDs into a net vote count for each possible ID
    print '\nSumming Votes'

    # DEBUGGING
    print
    for ind in range(nStars):
        print 'Star ID ', ind, ': ', runningVoteCount[ind]
    print

    # container for net vote results
    # netIDs[0] = [4,3,6,8] ... measured id0 star has votes for reference ID's 4,3,6,8
    # netVoteCount[0] = [2,2,3,3] ... measured id0 star has vote counts of 2 for reference ids 4,3 and 3 votes for 6,8
    netIDs = []
    netVoteCount = []
    voteResults = []

    for indStar in range(nStars):

        # print '\nCHECK: Final Running Vote Count Results for star ', indStar, ': ', runningVoteCount[indStar]

        if runningVoteCount[indStar] != -1:

            uniqueIDs, counts = generate_uniqueID_and_counts(runningVoteCount[indStar])

            # Determine minimum number of votes required to be considered a match
            # voteMinCount = voteMinRatio * (sum(counts)/2.)
            voteMinCount = voteMinRatio * searchesPerObject[indStar]

            print 'Initial ID ',indStar,', Vote Result Max Counts of ', max(counts), ' for ID: ', \
                uniqueIDs[np.argmax(counts)], ' (positive ID requires ', voteMinCount, '/', searchesPerObject[indStar], ' searches'

            # first column are reference IDs, second column are vote counts
            netVotes = np.empty((0, 2))
            for indVote in range(len(uniqueIDs)):
                np.vstack((netVotes, np.array([uniqueIDs[indVote], counts[indVote]])))

            # Store positive ID matches or None if not found
            if max(counts) >= voteMinCount and len(runningVoteCount[indStar]) > 2:

                # voteResults.append(netVotes)
                netIDs.append(uniqueIDs[np.argmax(counts)])
                # netVoteCount.append(max(counts))

            else:
                netIDs.append(None)

        else:
            netIDs.append(None)

    # process votes to determine likeliest candidate

    # list of max vote count each observed star
    # maxVote = []
    # for indStar in range(nStars):
    #     maxVote[indStar] = max(netVoteCount[indStar])

    # store max vote count
    # for indStar in range(nStars):
    #     if maxVote[indStar] < voteMin:

    objectID = netIDs

    return objectID


##################################################
##################################################

def catalogIDsToRaDec(refIDs, fnameCatalog):

    ref_catalog = sql.connect(fnameCatalog)
    with ref_catalog:

        radec = []
        cursor = ref_catalog.cursor()

        for indID in range(len(refIDs)):

            cursor.execute("SELECT * FROM tycho_data WHERE id IS (?)", (refIDs[indID],))

            catalogMatch = cursor.fetchall()

            radec.append((catalogMatch[0][1], catalogMatch[0][2]))

    return radec

