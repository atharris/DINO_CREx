# data_generation.py

import sys, os, inspect

filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))
path2 = os.path.dirname(os.path.abspath(filename))
bskName = 'Basilisk'
dinoName = 'DINO_CREx'
splitPath = path.split(dinoName)
dinoSpicePath = splitPath[0] + dinoName + '/DINObatch/SPICE/'
bskSpicePath = splitPath[0] + bskName + '/External/EphemerisData/'
bskPath = splitPath[0] + bskName + '/'
sys.path.append(bskPath + 'modules')
sys.path.append(bskPath + 'PythonModules')
sys.path.append(dinoSpicePath)

import pyswice
import numpy as np
from batchFilter import runRef




def generate_data(sc_ephem_file, planet_beacons,
                  beacon_ids, n_observations=10,
                  start_et=None, end_et=None, extras = {}, realData = 'OFF', ref='J2000'):
    # Load Ephemeris Files
    pyswice.furnsh_c(sc_ephem_file)
    pyswice.furnsh_c(bskSpicePath + 'de430.bsp')
    pyswice.furnsh_c(dinoSpicePath + 'epoch_11JUL2020.bsp')
    pyswice.furnsh_c(dinoSpicePath + 'naif0011.tls')

    for beacon in planet_beacons:
        found = pyswice.new_intArray(1)
        pyswice.intArray_setitem(found, 0, 0)
        beaconid = pyswice.new_intArray(1)
        pyswice.bodn2c_c(beacon, beaconid, found)
        beacon_ids.append(pyswice.intArray_getitem(beaconid, 0))

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
        state = pyswice.new_doubleArray(6)
        lt = pyswice.new_doubleArray(1)
        stateArray = np.zeros(6)
        pyswice.spkezr_c(body_id_str, t, ref, 'None', 'SUN', state, lt)
        for i in range(6):
            stateArray[i] = pyswice.doubleArray_getitem(state, i)
        sc_states.append(stateArray)
    ephemerides = {'spacecraft': np.array(sc_states).T}

    dyn_ref_state = []
    if realData == 'OFF':
        IC0 = sc_states[0]
        n_state = IC0.shape[0]
        phi0 = np.identity(n_state)
        # input to the propagator takes the ref_state and STM at t0, as well as the list of times
        prop_input = (IC0, phi0, observation_times, extras)

        # execute propagation
        state = runRef(prop_input)
        for jj in range(len(observation_times)):
            dyn_ref_state.append(state[jj][0:n_state])
        ephemerides = {'spacecraft': np.array(dyn_ref_state).T}


    # for planet in planet_beacons:
    #     planet = planet.lower()
    #
    #     states = []
    #     for t in observation_times:
    #         stateArray = np.zeros(6)
    #         state = pyswice.new_doubleArray(6)
    #         lt = pyswice.new_doubleArray(1)
    #         pyswice.spkezr_c(body_id_str, t, ref, 'None', 'SUN', state, lt)
    #         for i in range(6):
    #             stateArray[i] = pyswice.doubleArray_getitem(state, i)
    #         sc_states.append(stateArray)
    #
    #     ephemerides[planet] = np.array(states).T

    for beacon_id in beacon_ids:
        beacon_states = []
        for t in observation_times:
            stateArray = np.zeros(6)
            state = pyswice.new_doubleArray(6)
            lt = pyswice.new_doubleArray(1)
            pyswice.spkezr_c(str(beacon_id), t, ref, 'None', 'SUN', state, lt)
            for i in range(6):
                stateArray[i] = pyswice.doubleArray_getitem(state, i)
            beacon_states.append(stateArray)
        beacon_states = np.array(beacon_states).T
        ephemerides[str(beacon_id)] = np.array(beacon_states)
    return ephemerides, observation_times

# generate_data(sc_ephem_file='Ephem1.bsp', beacon_ephem_files=[], n_observations=100,
#               start_et=SP.utc2et(start_et), end_et=SP.utc2et(end_et))
