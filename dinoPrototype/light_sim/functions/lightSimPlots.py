#	Title   : lightSimPlots.py
#	Author  : Joe Park
#	Date    : 03/28/2017
#	Synopsis: Plotting functions for the lighting simulation module.

from __future__ import division

import math
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick


#######################################
#######################################

def axisEqual3D(ax):
    extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = extents[:,1] - extents[:,0]
    centers = np.mean(extents, axis=1)
    maxsize = max(abs(sz))
    r = maxsize/2
    for ctr, dim in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)

    ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))
    ax.zaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))

    ax.xaxis.labelpad = 20
    ax.yaxis.labelpad = 20
    ax.zaxis.labelpad = 20

    ax.tick_params(axis='x', pad=5)
    ax.tick_params(axis='y', pad=5)
    ax.tick_params(axis='z', pad=5)



######################################
######################################

def plotCB_cart(xyz,flux):

    x = xyz[:,0]
    y = xyz[:,1]
    z = xyz[:,2]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x,y,z,c='r',marker='.')
    ax.set_xlabel('I-helio')
    ax.set_ylabel('J-helio')
    ax.set_zlabel('K-helio')
    axisEqual3D(ax)
    plt.title('CB in Helio Coordinates')
    plt.show()



#######################################
#######################################

def plotCB_azel(azel,flux):

    az = azel[:,0]
    el = azel[:,1]
    xyz = np.empty((0,3), float)

    r = 1
    ind_el = 0
    for current_az in az:
        current_el = el[ind_el]
        xy = r*math.acos(math.radians(current_el))
        x = xy*math.atan(math.radians(current_az))
        y = xy*math.atan(math.radians(current_az))
        z = r*math.asin(math.radians(current_el))
        xyz = np.vstack((xyz,np.array([x,y,z])))
        ind_el += 1

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(az, el, flux, c='r',marker='.')
    ax.set_xlabel('azimuth')
    ax.set_ylabel('elevation')
    ax.set_zlabel('flux')
    plt.title('Flux vs. Camera Azimuth Elevation')
    plt.show()



###########################################
###########################################

def plotCB_flux(xyz, flux):

    AU = 149597870.7E3  # kxm

    normflux = 1-(flux/max(flux))

    x = xyz[:,0]
    y = xyz[:,1]
    z = xyz[:,2]
    # fig = plt.figure()
    plt.style.use('classic')
    fig = plt.figure(facecolor="1")
    ax = fig.add_subplot(111, projection='3d', axisbg='w')
    ax.scatter(x,y,z,c=normflux, marker='.', cmap='gray', lw=0.001)
    ax.set_axis_bgcolor('white')
    ax.set_xlabel('I-helio [km]')
    ax.set_ylabel('J-helio [km]')
    ax.set_zlabel('K-helio [km]')
    axisEqual3D(ax)
    ax.view_init(20, -60)
    #plt.gray()
    plt.title('CB in Helio Coordinates')
    plt.show()



def plotCB_fluxIJ(xyz, flux):

    AU = 149597870.7  # km

    normflux = 1-(flux/max(flux))

    x = xyz[:,0]
    y = xyz[:,1]
    z = xyz[:,2]
    # fig = plt.figure()
    plt.style.use('classic')
    fig = plt.figure(facecolor="1")
    ax = fig.add_subplot(111, axisbg='w')
    #ax = fig.add_subplot(111, projection='3d', axisbg='w')
    ax.scatter(x,y,c=normflux, marker='.', cmap='gray', lw=0.001)
    ax.set_axis_bgcolor('white')
    ax.set_xlabel('I-helio [km]')
    ax.set_ylabel('J-helio [km]')
    #ax.set_zlabel('K-helio [km]')
    #axisEqual3D(ax)
    #ax.view_init(20, -60)
    #plt.gray()
    plt.title('CB in Helio Coordinates')
    plt.show()



def plotCBandObs(pos_obs, pos_cb, dcm_obs):

    AU = 149597870.7  # km

    pos_obs = pos_obs/AU
    pos_cb = pos_cb/AU

    plt.style.use('classic')
    fig = plt.figure(facecolor="1")
    ax = fig.add_subplot(111, projection='3d', axisbg='w')
    ax.scatter(x,y,z,c=normflux, marker='.', cmap='gray', lw=0.001)
    ax.scatter(pos_obs[0], pos_obs[1], pos_obs[2],c='r', marker='o', lw=1)
    ax.scatter(pos_cb[0], pos_cb[1], pos_cb[2],c='b', marker='o', lw=1)

    ax.set_axis_bgcolor('white')
    ax.set_xlabel('I-helio [AU]')
    ax.set_ylabel('J-helio [AU]')
    ax.set_zlabel('K-helio [AU]')
    axisEqual3D(ax)
    ax.view_init(90, 0)
    #plt.gray()
    plt.title('Target Beacon and Observer Location')
    plt.show()


