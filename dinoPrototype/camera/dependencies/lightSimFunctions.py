#	Title   : DINO_solar_model.py
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


def mapSphere(lat_res, long_res, rad_cb):
    """Generate latitude and longitude coordinates for semi-spherical mapping. Method assumes equal divisions in
    latitude and equidistant divisions in longitude (number of longitudinal pts at each latitude dependent on horizontal
    radius of celestial body at a given latitude).
    Assumptions: surface mapping is fine enough to approximate surface area as 2D rectangle
    :param lat_res: number of latitude points
    :param long_res: number of surface map pts along zero latitude
    :param rad_cb: radius of celestial body [m]
    :return N x 2 array of latitude and longitude coordinate for semi-spherical mapping [degrees]
            2D Surface area of a single rectangular facet [m^2]"""
    import pdb
    from numpy import zeros, cos, deg2rad, rad2deg
    from datetime import datetime
    start_ms = datetime.now()
    # delta latitude
    delt_lat = (180./lat_res)

    # delta longitude at equator and distance between surface mapping pts
    delt_long_equator = (180./long_res)
    delt_horiz = rad_cb*math.radians(delt_long_equator)

    # only using semi-sphere
    pts_lat = np.arange(-90.,90.+delt_lat,delt_lat)
    pts_latlong = np.empty((0,2),float)
    pts_latlong_c = np.empty((0,2),float)

    lat_c = pts_lat[abs(pts_lat) < 90] - .5 * delt_lat
    r_cb = rad_cb*cos(deg2rad(lat_c))
    delt_long = rad2deg((delt_horiz/r_cb))
    npts_long, long_remainder = divmod(180., delt_long)

    # pdb.set_trace()
    # loop through each latitude and calculate longitudinal pts
    # (builds array of latitude and longitude pts of facet vertices)
    for current_lat in np.nditer(pts_lat):

        current_lat_c = current_lat - .5 * delt_lat

        if abs(current_lat_c) < 90:

            # horizontal semicircle perimeter at current latitude
            current_r_cb = rad_cb*math.cos(math.radians(current_lat_c))

            # calculate number of surface map pts at current longitude and remainder
            current_delt_long = math.degrees((delt_horiz/current_r_cb))
            current_npts_long, long_remainder = divmod(180., current_delt_long)
            long_i = long_remainder/2.
            long_f = long_i + (current_npts_long*current_delt_long)
            # pdb.set_trace()
            if current_npts_long == 0:
                current_pts_long_c = np.array([0.0])
            else:
                current_pts_long_c = np.arange(
                    long_i+(current_delt_long/2.), 
                    long_f, current_delt_long
                    )

            # create array of latlong for facet vertices

            tmp_stack = np.vstack(
                (
                    zeros(len(current_pts_long_c)) + current_lat,
                    current_pts_long_c
                    )
                ).T

            pts_latlong = np.vstack((pts_latlong,tmp_stack))
            # for current_long_c in np.nditer(current_pts_long_c):
            #     pts_latlong = np.vstack((pts_latlong, np.array([current_lat, current_long_c])))

            # create array of latlong for facet center locations
            # latitude value for center of facet

            # tmp_stack = np.vstack(
            #     (
            #         zeros(len(current_pts_long_c)) + current_lat_c,
            #         current_pts_long_c
            #         )
            #     ).T
            # pts_latlong_c = np.vstack((pts_latlong_c,tmp_stack))
            # for current_long_c in np.nditer(current_pts_long_c):
            #     if current_long_c < 180:
            #         pts_latlong_c = np.vstack((pts_latlong_c, np.array([current_lat_c, current_long_c])))
    
    pts_latlong_c = pts_latlong
    pts_latlong_c[:,0] -= .5 * delt_lat

    print("mapSphere: " + str(datetime.now() - start_ms))
    # assume each facet has equal areas totaling seim-sphere
    npts_latlong_c = pts_latlong_c.shape[0]
    facet_area = (4.*math.pi*rad_cb**2) / npts_latlong_c

    mkplots = False
    if mkplots:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(pts_latlong[:, 0], pts_latlong[:, 1], c='b', marker='.')
        ax.scatter(pts_latlong_c[:, 0], pts_latlong_c[:, 1], c='y', marker='*')
        ax.set_xlabel('Latitude (deg)')
        ax.set_ylabel('Longitude (deg)')
        plt.grid(True)
        plt.show()
    print("mapSphere: " + str(datetime.now() - start_ms))
    return pts_latlong_c, facet_area

###################################################
###################################################


def checkFoV(pos_cb, pos_obs, DCM_BN, fov, radius_cb):
    """Check if celelstial body is in camera field of view
    :param pos_cb: position of celestial body in heliocentric coordinates [m]
    :param pos_obs: position of camera in heliocentric coordinates [m]
    :param DCM_NB: direction cosine matrix of observer attitude [body to heliocentric coord. frame]
    :param field of view tuple (horizontal, vertical) [degrees]
    :return: true if in field of view, false otherwise
    """
    from datetime import datetime
    start = datetime.now()

    # compute relative position of celestial body to observer in body coordinates
    pos_obs2cb_helio = (pos_cb-pos_obs)
    pos_obs2cb_body = np.matmul(DCM_BN, pos_obs2cb_helio)
    e_obs2cb_body = pos_obs2cb_body / np.linalg.norm(pos_obs2cb_body)

    # take field of view angles from centerline
    horiz_fovcenter = fov[0]/2.
    vert_fovcenter = fov[1]/2.

    # compute azimuth and elevation of celestial body
    az = math.degrees(math.atan2(e_obs2cb_body[1], e_obs2cb_body[0]))
    el = math.degrees(math.asin(e_obs2cb_body[2]))

    # check if limits in camera specs
    chk = True
    if (abs(az) > horiz_fovcenter) or (abs(el) > vert_fovcenter):
        chk = False

    if np.linalg.norm(pos_obs2cb_body) < radius_cb:
        chk = False
    print("checkFOV: " + str(datetime.now() - start))

    return chk


#####################################
#####################################


def lumos(pos_cb, pos_obs, albedo_cb, rad_cb, lat_res, long_res):
    """Simulate illumination of a celestial body assuming constant geometric albedo on a spherical surface
    :param pos_cb: position of celestial body in heliocentric coordinates [m]
    :param pos_obs: position of observer in heliocentric coordinates [m]
    :param albedo_cb: geometric albedo of celestial body (pure diffuse reflection)
    :param rad_cb: spherical radius of celestial body [km]
    :return:    array of flux reduction at each surface point (1 x N) [unitless],
                array of cartesian coordinates for surface mapping (N x 3) in heliocentric coordinates [m]
                square facet_area [m^2]    """

    from numpy import deg2rad, sin, cos, vstack, pi
    from datetime import datetime
    start = datetime.now()
    print(start)
    # Constants
    FLUX_ref_distance = 695700. # radius of the sun in km
    ########################   ############
    # Compute transformation matrix from body to heliocentric coordinates
    # celestial body coordinate frame: {+i towards sun, j (k cross i), +k upwards}

    # unit vector to sun from celestial body center in heliocentric coordinates
    distance_cb = np.linalg.norm(pos_cb)
    ibody_helio = -pos_cb/ distance_cb

    # unit vector normal to i_body in the i-j body frame
    ijnormbody_helio = np.array([-ibody_helio[1], ibody_helio[0], 0])
    ijnormbody_helio = ijnormbody_helio/np.linalg.norm(ijnormbody_helio)
    kbody_helio = np.cross(ibody_helio,ijnormbody_helio)
    kbody_helio = kbody_helio / np.linalg.norm(kbody_helio)
    jbody_helio = np.cross(kbody_helio, ibody_helio)
    jbody_helio = jbody_helio/ np.linalg.norm(jbody_helio)

    # direction cosine matrix for body to helio coordinate frame
    dcm_NB = (np.array([ibody_helio,jbody_helio,kbody_helio])).T

    ####################################
    # Generate surface map of celestial body in body frame coordinates
    # create spherical grid using equal area facets
    pts_latlong, facet_area = mapSphere(lat_res, long_res, rad_cb)
    n_latlong = pts_latlong.shape[0]


    # generate unit vector of surface normal vector for each surface point (in body frame)

    pts_e_norm_body = np.zeros((n_latlong,3))

    latitudes = pts_latlong[:,0]
    longitudes = pts_latlong[:,1]

    e_norm_kbody = sin(deg2rad(latitudes))
    e_norm_ijbody = cos(deg2rad(latitudes))
    e_norm_ibody = sin(deg2rad(longitudes))*e_norm_ijbody
    e_norm_jbody = cos(deg2rad(longitudes))*e_norm_ijbody
    pts_e_norm_body = vstack([e_norm_ibody,e_norm_jbody,e_norm_kbody]).T

    # for ind_latlong in range(n_latlong):
    #     current_lat = pts_latlong[ind_latlong,0]
    #     current_long = pts_latlong[ind_latlong,1]

    #     e_norm_kbody = math.sin(math.radians(current_lat))
    #     e_norm_ijbody = math.cos(math.radians(current_lat))

    #     e_norm_ibody = math.sin(math.radians(current_long))*e_norm_ijbody
    #     e_norm_jbody = math.cos(math.radians(current_long))*e_norm_ijbody
    #     pts_e_norm_body[ind_latlong,:] = [e_norm_ibody, e_norm_jbody, e_norm_kbody]

    # generate surface map of ijk heliocentric coordinates from center of celestial body
    pts_helio = np.zeros((n_latlong,3))

    zbody = rad_cb*sin(deg2rad(latitudes))
    xybody = rad_cb*cos(deg2rad(latitudes))
    ybody = xybody*cos(deg2rad(longitudes))
    xbody = xybody*sin(deg2rad(longitudes))
    xyzbody = vstack([xbody,ybody,zbody]).T
    pts_helio = np.matmul(dcm_NB,xyzbody.T).T
    ind_pts = 0

    # for ind_latlong in range(n_latlong):

    #     current_lat = pts_latlong[ind_latlong,0]
    #     current_long = pts_latlong[ind_latlong,1]

    #     # determine surface point body coordinates
    #     zbody = rad_cb * math.sin(math.radians(current_lat))
    #     xybody = rad_cb * math.cos(math.radians(current_lat))
    #     ybody = xybody*math.cos(math.radians(current_long))
    #     xbody = xybody*math.sin(math.radians(current_long))
    #     xyzbody = np.array([xbody,ybody,zbody])

    #     # convert to heliocentric coordinates
    #     pts_helio[ind_pts,:] = np.matmul(dcm_NB,xyzbody)

    #     ind_pts += 1

    ###################################################
    # Calculate net flux reduction due to diffuse reflection and sun to CB inverse square law

    # calculate geometric albedo using cosine law

    pts_albedo = np.zeros((n_latlong, 1))
    pts_albedo = np.zeros(n_latlong)
    pts_albedo[pts_e_norm_body[:,0] > 0] = \
        pts_e_norm_body[:,0][pts_e_norm_body[:,0] > 0]
    pts_albedo = pts_albedo.reshape(len(pts_albedo),1)*albedo_cb

    # ind_latlong = 0
    # for  current_e_norm in pts_e_norm_body:

    #     phase_solar = math.acos(np.dot([1,0,0],current_e_norm))

    #     if phase_solar > math.pi/2:
    #         pts_albedo[ind_latlong] = 0.
    #     else:
    #         pts_albedo[ind_latlong] = math.cos(phase_solar)*albedo_cb

    #     ind_latlong += 1

    # calculate flux decay due to inverse square laws
    flux_decay_sun2cb = (FLUX_ref_distance/distance_cb)**2
    distance_cb2obs = np.linalg.norm(pos_cb-pos_obs)
    flux_decay_cb2obs = (rad_cb/distance_cb2obs)**2
    flux_decay_cb2obs = (1/distance_cb2obs)**2
    flux_decay_net = flux_decay_sun2cb * pts_albedo * flux_decay_cb2obs
    print("lumos: " + str(datetime.now()-start))
    return pts_helio, flux_decay_net, facet_area

###################################################
###################################################


def project2CamView(pos_cb, pos_obs, attde_cam, xyz_helio, flux_decay, fov, facet_area):
    """Remove celestial body surface map points that are out of the camera view and convert heliocentric cartesian
    coordinates to azimuth and elevation from camera point of view
    :param pos_cb: position of celestial body in heliocentric coordinates [m]
    :param pos_obs: position of camera in heliocentric coordinates [m]
    :param DCM_NB: direction cosine matrix of observer attitude [body to heliocentric coord. frame]
    :param field of view tuple (horizontal, vertical) [degrees]
    :param facet_area: area of single facet on spherical approximation [m^2]
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
    dcm_BN = attde_cam
    pos_obs2cb_helio = (pos_cb-pos_obs)
    e_obs2cb_helio = pos_obs2cb_helio / np.linalg.norm(pos_obs2cb_helio)
    pos_obs2facet_helio = pos_obs2cb_helio + xyz_helio
    r_ob2facet_helio =  sqrt(
        pos_obs2facet_helio[:,0]**2+pos_obs2facet_helio[:,1]**2+pos_obs2facet_helio[:,2]**2
        )

    r_ob2facet_helio = r_ob2facet_helio.reshape(len(r_ob2facet_helio),1)
    e_obs2facet_helio = pos_obs2facet_helio/r_ob2facet_helio
    pos_obs2cb_body = np.matmul(dcm_BN, pos_obs2cb_helio)
    # cycle through surface points and calculate dot product with unit vector from camera to CB
    # (in helio coordinates), input points with phase angle > 90 deg (visible to camera) in new arrays

    horiz_camfov = fov[0] / 2.
    vert_camfov = fov[1] / 2.

    # xyz_camview = np.empty((0,3), float)
    # xyz_camview_helio = np.empty((0,3), float)
    # facet_area_camview = np.empty((0,1), float)

    # azel_pts = np.empty((0,2), float)
    # flux_decay_out = np.empty((0,1), float)
    # xyz_vis_cam = np.empty((0,3), float)


    r = sqrt(xyz_helio[:,0]**2+xyz_helio[:,1]**2+xyz_helio[:,2]**2).reshape(
        len(xyz_helio),1)
    e_xyz_helio = xyz_helio/r
    # cos_cam_phase2 = np.dot(e_obs2cb_helio, e_xyz_helio.T)
    cos_cam_phase = einsum('ij,ji->i',e_obs2facet_helio,e_xyz_helio.T)

    cam_phase = arccos(cos_cam_phase)

    # check if current facet is visible to observer
    ind = cam_phase > math.pi/2
    e_xyz_helio = e_xyz_helio[ind]
    cos_cam_phase = cos_cam_phase[ind]
    cam_phase = cam_phase[ind]
    xyz_helio = xyz_helio[ind]
    xyz_body = np.matmul(dcm_BN, xyz_helio.T)
    r_pt = pos_obs2cb_helio + xyz_body.T
    r_pt_norm = sqrt(r_pt[:,0]**2+r_pt[:,1]**2+r_pt[:,2]**2)
    az = rad2deg(arctan2(r_pt[:,1], r_pt[:,0]))
    el = rad2deg(arcsin(r_pt[:,2] / r_pt_norm))
    xyz_camview = xyz_body
    xyz_camview_helio = xyz_helio
    azel_pts = vstack((az,el)).T
    flux_decay_out = flux_decay[ind]
    xyz_vis_cam = xyz_body.T
    # current_facet_proj_area = 
    facet_area_camview = facet_area * -cos_cam_phase
    facet_area_camview = facet_area_camview.reshape(
        len(facet_area_camview),1)

    xyz_camview_helio = np.matmul(np.linalg.inv(dcm_BN),xyz_camview_helio.T).T
    xyz_vis_cam = np.matmul(np.linalg.inv(dcm_BN),xyz_vis_cam.T).T
    # plt.plot(xyz_camview_helio[:,0],xyz_camview_helio[:,1],'.')
    # plt.plot(-pos_obs2cb_helio[0],-pos_obs2cb_helio[1],'.')

    # ind_npts = 0
    # for current_xyz_helio in xyz_helio:

    #     e_current_xyz_helio = current_xyz_helio / np.linalg.norm(current_xyz_helio)
    #     cos_cam_phase = np.dot(e_obs2cb_helio, e_current_xyz_helio)
    #     cam_phase = math.acos(cos_cam_phase)

    #     # check if current facet is visible to observer
    #     if cam_phase > math.pi/2:

    #         xyz_body = np.matmul(dcm_BN, current_xyz_helio.T)

    #         r_pt = pos_obs2cb_body + xyz_body

    #         az = math.degrees(math.atan2(r_pt[1], r_pt[0]))
    #         el = math.degrees(math.asin(r_pt[2] / np.linalg.norm(r_pt)))

    #         # # check if current facet is in camera field of view
    #         # if (abs(az) < horiz_camfov) & (abs(el) < vert_camfov):

    #         xyz_camview = np.vstack((xyz_camview, xyz_body))
    #         xyz_camview_helio = np.vstack((xyz_camview_helio, current_xyz_helio))

    #         azel_pts = np.vstack((azel_pts, np.array([[az, el]])))
    #         flux_decay_out = np.vstack((flux_decay_out, flux_decay[ind_npts]))
    #         xyz_vis_cam = np.vstack((xyz_vis_cam, xyz_body))

    #         current_facet_proj_area = facet_area * -cos_cam_phase
    #         facet_area_camview = np.vstack((facet_area_camview, current_facet_proj_area))

    #     ind_npts +=1
    print("project2CamView: " + str(datetime.now() - start))

    return azel_pts, xyz_camview_helio, xyz_vis_cam, flux_decay_out, facet_area_camview



###################################################
###################################################

def xyz2RADec(xyz_cb):
    """Calculate right ascension and declination from xyz coordinates
        :param xyz_cb: 1x3 array of position
        :return ra_dec: 1x2 array of center of CB in right ascension and declination
                ra [deg] is angle from [1,0,0] to point in the xy horizontal plane
                dec [deg] is angle above xy plane to point"""

    from numpy import arctan2, arcsin, sqrt, rad2deg
    from datetime import datetime
    start = datetime.now()
    if len(xyz_cb.shape) == 1:

        r = np.linalg.norm(xyz_cb)
        ra = math.degrees(math.atan2(xyz_cb[1], xyz_cb[0]))
        dec = math.degrees(math.asin(xyz_cb[2]/r))

    else:


        r = sqrt(xyz_cb[:,0]**2+xyz_cb[:,1]**2+xyz_cb[:,2]**2)
        ra = rad2deg(arctan2(xyz_cb[:,1],xyz_cb[:,0]))
        dec = rad2deg(arcsin(xyz_cb[:,2]/r))
        ra = ra.reshape(len(ra),1)
        dec = dec.reshape(len(dec),1)

        # ra = []
        # dec = []

        # for current_xyz in xyz_cb:

        #     current_r = np.linalg.norm(current_xyz)
        #     current_ra = math.degrees(math.atan2(current_xyz[1], current_xyz[0]))
        #     current_dec = math.degrees(math.asin(current_xyz[2]/current_r))

        #     ra.append(current_ra)
        #     dec.append(current_dec)

        # ra = np.array(ra).reshape(len(ra),1)
        # dec = np.array(dec).reshape(len(dec),1)

    ra_dec = (ra, dec)
    print("xyz2RADec: " + str(datetime.now() - start))
    return ra_dec


###################################################
###################################################



def lightSim(attde_cam, pos_cam, pos_cb, fov, lat_res, long_res, do_pt_source,
    albedo_cb, radius_cb, name_cb):
    """Conduct lighting simulation of multiple bodies.
    :param pos_cb: Nx3 array of positions of celestial bodies in heliocentric coordinates [m]
    :param pos_obs: position of camera in heliocentric coordinates [m]
    :param attde_cam: direction cosine matrix of observer attitude [body to heliocentric coord. frame]
    :param field of view tuple (horizontal, vertical) [degrees]
    :param do_pt_source: bool whether to treat CB as a point source
    :return: dict of summary output for each visible CB with cb_name as key
             summary output dict has the following values
            {'bodypos': (npts x 3) numpy array [m]
            'facetRA' : right ascension of facet in heliocentric coord. frame [deg]
            'facetDec': declination of facet in heliocentric coord. frame [deg]
            'fluxDecay': flux decay due to inverse square laws and diffuse reflection
            'facetArea': projected facet area normal to observer [m^2]
    """
    from numpy import array
    lightSim_output = {}

    # cycle through each CB
    ind_cb = 0


    # albedo_cb, radius_cb, name_cb = getCBparam(ind_cb)

    # check if CB is in field of view
    fov_chk = checkFoV(pos_cb, pos_cam, attde_cam, fov, radius_cb)
    # compute lighting simulation
    if fov_chk:

        current_CB_dict = {}

        # illuminate sphere
        pts_pos_helio, pts_albedo_helio, facet_area = lumos(
            pos_cb, pos_cam, albedo_cb, radius_cb, lat_res, long_res)

        # limit points to those in camera field of view and compute azel
        pts_azel_cam, pts_xyz_helio, pts_xyz_cam, pts_albedo_cam, facet_area_cam = project2CamView(
            pos_cb, pos_cam, attde_cam, pts_pos_helio, pts_albedo_helio, fov, facet_area)

        # check to see if there are any illuminated points visible for the celestial body
        # (celestial body may be in field of view with no illuminated points (observer on dark side of body)
        if pts_albedo_cam.shape[0] > 0:

            if do_pt_source:

                # compute position of CB center in camera view
                dcm_BN = attde_cam
                pos_obs2cb_helio = (pos_cb - pos_cam)
                print(pos_obs2cb_helio)
                pos_obs2cb_body = np.matmul(dcm_BN, pos_obs2cb_helio.T)

                ra_dec = xyz2RADec(pos_obs2cb_body)

                #bodypos needs to be the positon of this facet relative
                #to the body. Because a point source facet will always
                #be in the same place as the center of its body,
                #current_CB_dict['bodypos'] should always be zero for
                #a point source 
                current_CB_dict['bodypos'] = array([[0,0,0]])
                current_CB_dict['facetRA'] = array([ra_dec[0]])
                current_CB_dict['facetDec'] = array([ra_dec[1]])
                current_CB_dict['netAlbedo'] = pts_albedo_cam
                current_CB_dict['facetArea'] = facet_area_cam

            else:

                ra_dec = xyz2RADec(pts_xyz_helio)

                current_CB_dict['bodypos'] = pts_xyz_cam
                current_CB_dict['facetRA'] = ra_dec[0]
                current_CB_dict['facetDec'] = ra_dec[1]
                current_CB_dict['netAlbedo'] = pts_albedo_cam
                current_CB_dict['facetArea'] = facet_area_cam


        import matplotlib.pyplot as plt
    else: 
        current_CB_dict = -1

    return current_CB_dict

