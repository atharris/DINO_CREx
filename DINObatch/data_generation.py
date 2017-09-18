# data_generation.py

import os
import pickle
import sys
import time
import optparse
import socket
import pdb
import scipy.integrate as integ
import scipy.io as io
import spiceypy as SP
import matplotlib.pyplot as plt
import numpy as np


def generate_data(sc_ephem_file, planet_beacons,
                  beacon_ids, n_observations=10,
                  start_et=None, end_et=None, ref='J2000'):
    # Load Ephemeris Files
    SP.furnsh(sc_ephem_file)
    SP.furnsh('SPICE/de430.bsp')
    SP.furnsh('SPICE/epoch_11JUL2020.bsp')
    # SP.furnsh('SPICE/mar097.bsp')
    SP.furnsh('SPICE/naif0011.tls')

    #for beacon_id in beacon_ids:
    #    SP.furnsh('SPICE/' + str(beacon_id) + '.bsp')

    # Identify Spacecraft Body
    body_id = -100  # SP.spkobj(sc_ephem_file)
    body_id_str = str(body_id)
    # Create Time Array
    observation_times = np.linspace(start_et, end_et, num=n_observations,
                                    endpoint=True)  # range(start_et, end_et, n_observations - 1)

    # Create Beacon and Spacecraft states at Observation Times
    # sc_states = []
    # earth_states = []
    # mars_states = []
    # jupiter_states = []
    # venus_states = []
    # for t in observation_times:
    #     earth_states.append(SP.spkezr(targ='3', et=t, ref=ref, abcorr='None', obs='SUN')[0])
    #     mars_states.append(SP.spkezr(targ='4', et=t, ref=ref, abcorr='None', obs='SUN')[0])
    #     jupiter_states.append(SP.spkezr(targ='5', et=t, ref=ref, abcorr='None', obs='SUN')[0])
    #     sc_states.append(SP.spkezr(targ=body_id_str, et=t, ref=ref, abcorr='None', obs='SUN')[0])
    #     venus_states.append(SP.spkezr(targ='2', et=t, ref=ref, abcorr='None', obs='SUN')[0])
    #
    # ephemerides = {'spacecraft': np.array(sc_states).T,
    # 'earth': np.array(earth_states).T,
    # 'mars': np.array(mars_states).T,
    # 'jupiter': np.array(jupiter_states).T,
    # 'venus': np.array(venus_states).T
    # }

    sc_states = []
    for t in observation_times:
        sc_states.append(SP.spkezr(targ=body_id_str, et=t, ref=ref, abcorr='None', obs='SUN')[0])
    ephemerides = {'spacecraft': np.array(sc_states).T}

    planet_ids = {'mercury': '1', 'venus': '2', 'earth': '3', 'mars':'4', 'jupiter':'5'}
    for planet in planet_beacons:
        planet = planet.lower()

        states = []
        for t in observation_times:
            planet_id = planet_ids[planet]
            states.append(SP.spkezr(targ=planet_id, et=t, ref=ref, abcorr='None', obs='SUN')[0])

        ephemerides[planet] = np.array(states).T

    for beacon_id in beacon_ids:
        beacon_states = []
        for t in observation_times:
            beacon_states.append(
                SP.spkezr(targ=str(beacon_id), et=t, ref=ref, abcorr='None', obs='SUN')[0])
        beacon_states = np.array(beacon_states).T
        ephemerides[str(beacon_id)] = np.array(beacon_states)

    return ephemerides, observation_times

# generate_data(sc_ephem_file='Ephem1.bsp', beacon_ephem_files=[], n_observations=100,
#               start_et=SP.utc2et(start_et), end_et=SP.utc2et(end_et))
