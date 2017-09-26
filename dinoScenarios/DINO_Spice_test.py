import sys, os, inspect
from numpy import linalg as la
import numpy as np
import ctypes
import math
import csv
import logging
from datetime import datetime

bskName = 'Basilisk'
bskPath = '../..' + '/' + bskName + '/'
sys.path.append(bskPath + 'modules')
sys.path.append(bskPath + 'PythonModules')

# import general simulation support files
import SimulationBaseClass
import unitTestSupport
import matplotlib.pyplot as plt
import macros as mc

# import simulation related support
import spacecraftPlus
import gravityEffector
import ephemeris_converter
import spice_interface
import pyswice
import BSK_plotting as BSKPlt



#########################################################
# setup simulation containers
#########################################################

simulationTime = mc.sec2nano(139643.532)
#   Set the number of data points to be logged, and therefore the sampling frequency
numDataPoints = 10000
samplingTime = simulationTime / (numDataPoints - 1)

# Create simulation variable names
simTaskName = "simTask"
simProcessName = "simProcess"

# Create a sim module as an empty container
scSim = SimulationBaseClass.SimBaseClass()
scSim.TotalSim.terminateSimulation()

# create the simulation process
dynProcess = scSim.CreateNewProcess(simProcessName)

# create the dynamics task and specify the integration update time
simulationTimeStep = mc.sec2nano(5.)
dynProcess.addTask(scSim.CreateNewTask(simTaskName, simulationTimeStep))



#########################################################
# initialize spacecraftPlus object and set properties
#########################################################

scObject = spacecraftPlus.SpacecraftPlus()
scObject.ModelTag = "spacecraftBody"
scObject.hub.useTranslation = True
scObject.hub.useRotation = False

# add spacecraftPlus object to the simulation process
scSim.AddModelToTask(simTaskName, scObject, None, 1)



#########################################################
# setup gravity effectors
#########################################################

ephemConvert = ephemeris_converter.EphemerisConverter()
ephemConvert.ModelTag = "ephemerisConverter"

spicePlanetNames = []
gravBodyList = []

sunGravBody = gravityEffector.GravBodyData()
sunGravBody.bodyInMsgName = "sun_planet_data"
sunGravBody.outputMsgName = "sun_display_frame_data"
sunGravBody.mu = 1.32712440018E20 # meters^3/s^2
sunGravBody.radEquator = 695508000.0 # meters
sunGravBody.isCentralBody = True
sunGravBody.useSphericalHarmParams = False
gravBodyList.append(sunGravBody)
spicePlanetNames.append(sunGravBody.bodyInMsgName[:-12])

earthGravBody = gravityEffector.GravBodyData()
earthGravBody.bodyInMsgName = "earth_planet_data"
earthGravBody.outputMsgName = "earth_display_frame_data"
earthGravBody.mu = 0.3986004415E+15  # meters^3/s^2
earthGravBody.radEquator = 6378136.6  # meters
earthGravBody.isCentralBody = False
earthGravBody.useSphericalHarmParams = False
gravBodyList.append(earthGravBody)
spicePlanetNames.append(earthGravBody.bodyInMsgName[:-12])

# Current SPICE(SPKINSUFFDATA) error for Mars SPICE file
#marsGravBody = gravityEffector.GravBodyData()
#marsGravBody.bodyInMsgName = "mars_planet_data"
#marsGravBody.outputMsgName = "mars_display_frame_data"
#marsGravBody.mu = 0.3986004415E+15  # meters^3/s^2  UPDATE TO ACTUAL INFO
#marsGravBody.radEquator = 6378136.6  # meters  UPDATE TO ACTUAL INFO
#marsGravBody.isCentralBody = False
#marsGravBody.useSphericalHarmParams = False
#gravBodyList.append(marsGravBody)
#spicePlanetNames.append(marsGravBody.bodyInMsgName[:-12])

moonGravBody = gravityEffector.GravBodyData()
moonGravBody.bodyInMsgName = "moon_planet_data"
moonGravBody.outputMsgName = "moon_display_frame_data"
moonGravBody.mu = 0.3986004415E+15  # meters^3/s^2     UPDATE TO ACTUAL INFO
moonGravBody.radEquator = 6378136.6  # meters          UPDATE TO ACTUAL INFO
moonGravBody.isCentralBody = False
moonGravBody.useSphericalHarmParams = False
gravBodyList.append(moonGravBody)
spicePlanetNames.append(moonGravBody.bodyInMsgName[:-12])

beacon1GravBody = gravityEffector.GravBodyData()
beacon1GravBody.bodyInMsgName = "ceres_planet_data"
beacon1GravBody.outputMsgName = "ceres_display_frame_data"
beacon1GravBody.mu = 0.000  # meters^3/s^2
beacon1GravBody.radEquator = 1000.0  # meters
beacon1GravBody.isCentralBody = False
beacon1GravBody.useSphericalHarmParams = False
gravBodyList.append(beacon1GravBody)

beacon2GravBody = gravityEffector.GravBodyData()
beacon2GravBody.bodyInMsgName = "vesta_planet_data"
beacon2GravBody.outputMsgName = "vesta_display_frame_data"
beacon2GravBody.mu = 0.000  # meters^3/s^2
beacon2GravBody.radEquator = 5000.0  # meters
beacon2GravBody.isCentralBody = False
beacon2GravBody.useSphericalHarmParams = False
gravBodyList.append(beacon2GravBody)

scObject.gravField.gravBodies = spacecraftPlus.GravBodyVector(gravBodyList)



#########################################################
# setup Spice Objects
#########################################################

# setup simulation start date/time
timeInitString = "2020 JULY 11 00:28:30.0"
spiceTimeStringFormat = '%Y %B %d %H:%M:%S.%f'

timeInit = datetime.strptime(timeInitString, spiceTimeStringFormat)

# setup SPICE interface
spiceObject = spice_interface.SpiceInterface()
spiceObject.ModelTag = "SpiceInterfaceData"
spiceObject.SPICEDataPath = bskPath + 'External/EphemerisData/'
spiceObject.OutputBufferCount = 10000

# load beacon SPICE ephemeris data
# Ceres used as a test case <- issues with actual beacon spice files
beaconEphemerisName = 'Beacons/2000001.bsp'
beaconSpiceName = 'ceres'
spicePlanetNames.append(beaconSpiceName)
spiceObject.loadSpiceKernel(beaconEphemerisName, spiceObject.SPICEDataPath)

# Vesta used as a test case <- issues with actual beacon spice files
beaconEphemerisName2 = 'Beacons/2000004.bsp'
beaconSpiceName2 = 'vesta'
spicePlanetNames.append(beaconSpiceName2)
spiceObject.loadSpiceKernel(beaconEphemerisName2, spiceObject.SPICEDataPath)

# 1990MF Asteroid used as Target Beacon
# currently not working ... SpiceName is recognized but "insufficient ephemeris data has been loaded to compute ..."
#beaconEphemerisName2 = '2008014.bsp'
#beaconSpiceName2 = '8014'      
#spicePlanetNames.append(beaconSpiceName2)
#spiceObject.loadSpiceKernel(beaconEphemerisName2, spiceObject.SPICEDataPath)

# 2001 FB44 Asteroid used as Target Beacon
# currently not working ... SpiceName is recognized but "insufficient ephemeris data has been loaded to compute ..."
#beaconEphemerisName2 = '2054509.bsp'
#beaconSpiceName2 = '54509'
#spicePlanetNames.append(beaconSpiceName2)
#spiceObject.loadSpiceKernel(beaconEphemerisName2, spiceObject.SPICEDataPath)

ephemNames = [beaconEphemerisName, beaconEphemerisName2]
spiceNames = [beaconSpiceName, beaconSpiceName2]

# define PlanetNames vector for SpiceInterface ... make sure external spice kernels are loaded as well
spiceObject.PlanetNames = spice_interface.StringVector(spicePlanetNames)

# pull in SPICE support libraries
pyswice.furnsh_c(spiceObject.SPICEDataPath + 'de430.bsp')  # solar system bodies
pyswice.furnsh_c(spiceObject.SPICEDataPath + 'naif0011.tls')  # leap second file
pyswice.furnsh_c(spiceObject.SPICEDataPath + 'de-403-masses.tpc')  # solar system masses
pyswice.furnsh_c(spiceObject.SPICEDataPath + 'pck00010.tpc')  # generic Planetary Constants Kernel

spiceObject.UTCCalInit = timeInitString

# add spice interface object to task list
scSim.AddModelToTask(simTaskName, spiceObject, None, 2)

print 'Beacon Spice Filepath'
print spiceObject.SPICEDataPath + ephemNames[0]

# Set zero base of SPICE output as the sun
spiceObject.zeroBase = 'sun'
print '\nSpice zero base: ',  spiceObject.zeroBase


##########################################################
# initial spacecraft position and velocity vector
##########################################################

# Initial position of spacecraft set to Beacon location at timeInitString in Earth Coordinate Frame
#pyswice.furnsh_c(spiceObject.SPICEDataPath + ephemNames[0]) # Beacon spice file
#scInitialState = 1000*pyswice.spkRead(beaconSpiceName, timeInitString, 'J2000', 'EARTH')
#rN = scInitialState[0:3]         # meters
#vN = scInitialState[3:6]         # m/s

# Alternatively, explicitly set the s/c initial position and velocity
rN = np.array([53850779.24415602, -141892420.9599383, -61159.29323671013]) * 1000.
vN = np.array([32.95818246003602, 14.73457852769852, -0.5045251942794007]) * 1000.



##########################################################
# setup simulation time
##########################################################

simulationTime = mc.sec2nano(139643.532)
numDataPoints = 100
samplingTime = simulationTime / (numDataPoints-1)



##########################################################
# setup message logging
##########################################################

print 'Grav Body List: ', gravBodyList
print 'Spice Planet Names: ', spicePlanetNames
print 'spiceObject Planet Names: ', spiceObject.PlanetNames

# log spacecraftPlus message
scSim.TotalSim.logThisMessage(scObject.scStateOutMsgName, samplingTime)

# log spice messages
sunConvertName = 'sun_ephemeris_data'
earthConvertName = 'earth_ephemeris_data'
beacon1ConvertName = 'beacon1_ephemeris_data'
beacon2ConvertName = 'beacon2_ephemeris_data'

messageMap = {}
messageMap[sunGravBody.bodyInMsgName] = sunConvertName
messageMap[earthGravBody.bodyInMsgName] = earthConvertName

ephemConvert.messageNameMap = ephemeris_converter.map_string_string(messageMap)



##########################################################
# initialize simulation and initial spacecraft state (must be done after sim initialization)
##########################################################

scSim.InitializeSimulation()

posRef = scObject.dynManager.getStateObject("hubPosition")
velRef = scObject.dynManager.getStateObject("hubVelocity")

posRef.setState(unitTestSupport.np2EigenVectorXd(rN))  # m - r_BN_N
velRef.setState(unitTestSupport.np2EigenVectorXd(vN))  # m - v_BN_N



##########################################################
# configure a simulation stop time and execute the simulation run
##########################################################

scSim.ConfigureStopTime(simulationTime)
scSim.ExecuteSimulation()



##########################################################
# pull messages and print results
##########################################################

pos_sc = scSim.pullMessageLogData(scObject.scStateOutMsgName+'.r_BN_N',range(3))
pos_sun = scSim.pullMessageLogData(sunGravBody.bodyInMsgName + '.PositionVector', range(3))
pos_earth = scSim.pullMessageLogData(earthGravBody.bodyInMsgName + '.PositionVector', range(3))
pos_moon = scSim.pullMessageLogData(moonGravBody.bodyInMsgName + '.PositionVector', range(3))
print ''
print 'Beacon pullMessageLogData Start'
pos_beacon1 = scSim.pullMessageLogData(beacon1GravBody.bodyInMsgName + '.PositionVector', range(3))
pos_beacon2 = scSim.pullMessageLogData(beacon2GravBody.bodyInMsgName + '.PositionVector', range(3))
print 'Beacon pullMessageLogData Finish'

np.set_printoptions(precision=8)

print ''
print 'SpacecraftPlus Logged Positions'
print pos_sc[-1:,1:]
print 'Earth Logged Position'
print pos_earth[-1:, 1:]
print 'Sun Logged Position'
print pos_sun[-1:, 1:]
print 'Moon Logged Position'
print pos_moon[-1:, 1:]
print 'Ceres Logged Position'
print pos_beacon1[-1:, 1:]
print 'Vesta Logged Position'
print pos_beacon2[-1:, 1:]

print ''
scSim.TotalSim.PrintSimulatedMessageData()



##########################################################
# Plotting
##########################################################

# generate plots of output

print '\n\n'
print 'CELESTIAL:'
print 'r_sc = ', la.norm(pos_sc[-1:, 1:])
print 'r_earth = ', la.norm(pos_earth[-1:, 1:])
print 'r_sun = ', la.norm(pos_sun[-1:, 1:])
print 'r_moon = ', la.norm(pos_moon[-1:, 1:])

dict_data_color = {
    'moon': [pos_moon, 'cyan'],
    'earth': [pos_earth, 'dodgerblue'],
    'sun': [pos_sun, 'orange'],
    'ceres': [pos_beacon1, 'g'],
    'vesta': [pos_beacon2, 'salmon'],
    }

BSKPlt.plot_multi_orbit_0(dict_data_color)
#BSKPlt.plot_spacecraft_orbit_0(dict_data_color, pos_sc)

sc_dict_data_color = {
    'moon': [pos_moon, 'cyan'],
    'earth': [pos_earth, 'dodgerblue'],
    'ceres': [pos_beacon1, 'g'],
    'vesta': [pos_beacon2, 'salmon'],
    }

BSKPlt.plot_spacecraft_orbit(sc_dict_data_color, pos_sc)

plt.show()



##########################################################
# Simulation cleanup
##########################################################

pyswice.unload_c(spiceObject.SPICEDataPath + 'de430.bsp')
pyswice.unload_c(spiceObject.SPICEDataPath + 'naif0011.tls')
pyswice.unload_c(spiceObject.SPICEDataPath + 'de-403-masses.tpc')
pyswice.unload_c(spiceObject.SPICEDataPath + 'pck00010.tpc')
pyswice.unload_c(spiceObject.SPICEDataPath + beaconEphemerisName)
pyswice.unload_c(spiceObject.SPICEDataPath + beaconEphemerisName2)
