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


try:
    import pyswice
except ImportError:
    from Basilisk import pyswice
    bskSpicePath = splitPath[0] + bskName + '/supportData/EphemerisData/'
import numpy as np
from batchFilter import runRef




def generate_data(sc_ephem_file, planet_beacons,
                  beaconIDs, n_observations=10,
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
        beaconIDs.append(pyswice.intArray_getitem(beaconid, 0))

    # Identify Spacecraft Body
    bodyID = -100  # SP.spkobj(sc_ephem_file)
    bodyIDStr = str(bodyID)
    # Create Time Array
    observationTimes = np.linspace(start_et, end_et, num=n_observations,
                                    endpoint=True)  # range(start_et, end_et, n_observations - 1)

    # Create Beacon and Spacecraft states at Observation Times
    # stateSpacecraft = []
    # earth_states = []
    # mars_states = []
    # jupiter_states = []
    # venus_states = []
    # for t in observationTimes:
    #     earth_states.append(SP.spkezr(targ='3', et=t, ref=ref, abcorr='None', obs='SUN')[0])
    #     mars_states.append(SP.spkezr(targ='4', et=t, ref=ref, abcorr='None', obs='SUN')[0])
    #     jupiter_states.append(SP.spkezr(targ='5', et=t, ref=ref, abcorr='None', obs='SUN')[0])
    #     stateSpacecraft.append(SP.spkezr(targ=bodyIDStr, et=t, ref=ref, abcorr='None', obs='SUN')[0])
    #     venus_states.append(SP.spkezr(targ='2', et=t, ref=ref, abcorr='None', obs='SUN')[0])
    #
    # ephemerides = {'spacecraft': np.array(stateSpacecraft).T,
    # 'earth': np.array(earth_states).T,
    # 'mars': np.array(mars_states).T,
    # 'jupiter': np.array(jupiter_states).T,
    # 'venus': np.array(venus_states).T
    # }

    stateSpacecraft = []
    for t in observationTimes:
        state = pyswice.new_doubleArray(6)
        lt = pyswice.new_doubleArray(1)
        stateArray = np.zeros(6)
        pyswice.spkezr_c(bodyIDStr, t, ref, 'None', 'SUN', state, lt)
        for i in range(6):
            stateArray[i] = pyswice.doubleArray_getitem(state, i)
        stateSpacecraft.append(stateArray)
    ephemerides = {'spacecraft': np.array(stateSpacecraft).T}

    dyn_ref_state = []
    if realData == 'OFF':
        IC0 = stateSpacecraft[0]
        stateDimension = IC0.shape[0]
        phi0 = np.identity(stateDimension)
        # input to the propagator takes the ref_state and STM at t0, as well as the list of times
        propagationInputs = (IC0, phi0, observationTimes, extras)
        # execute propagation
        state = runRef(propagationInputs)
        for jj in range(len(observationTimes)):
            dyn_ref_state.append(state[jj][0:stateDimension])
        ephemerides = {'spacecraft': np.array(dyn_ref_state).T}


    # for planet in planet_beacons:
    #     planet = planet.lower()
    #
    #     states = []
    #     for t in observationTimes:
    #         stateArray = np.zeros(6)
    #         state = pyswice.new_doubleArray(6)
    #         lt = pyswice.new_doubleArray(1)
    #         pyswice.spkezr_c(bodyIDStr, t, ref, 'None', 'SUN', state, lt)
    #         for i in range(6):
    #             stateArray[i] = pyswice.doubleArray_getitem(state, i)
    #         stateSpacecraft.append(stateArray)
    #
    #     ephemerides[planet] = np.array(states).T

    for beaconID in beaconIDs:
        beaconState = []
        for t in observationTimes:
            stateArray = np.zeros(6)
            state = pyswice.new_doubleArray(6)
            lt = pyswice.new_doubleArray(1)
            pyswice.spkezr_c(str(beaconID), t, ref, 'None', 'SUN', state, lt)
            for i in range(6):
                stateArray[i] = pyswice.doubleArray_getitem(state, i)
            beaconState.append(stateArray)
        beaconState = np.array(beaconState).T
        ephemerides[str(beaconID)] = np.array(beaconState)
    return ephemerides, observationTimes

# generate_data(sc_ephem_file='Ephem1.bsp', beacon_ephem_files=[], n_observations=100,
#               start_et=SP.utc2et(start_et), end_et=SP.utc2et(end_et))
