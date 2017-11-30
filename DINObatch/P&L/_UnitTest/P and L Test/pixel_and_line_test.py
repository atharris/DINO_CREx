

import numpy as np
import pdb
import pixel_and_line_symbolics
from pixelLineBatch import fncH as H_function
from pixelLineBatch import fncG as G_function



## \defgroup Tutorials_1_3
##   @{
## How to setup orbital simulation with multiple gravitational bodies.
#
# Orbit Setup to Simulate Translation with Multiple Gravitational Bodies {#scenarioOrbitMultiBody}
# ====
#
# Scenario Description
# -----
# This script sets up a 3-DOF spacecraft which is traveling in a multi-gravity environment.  The purpose
# is to illustrate how to attach a multiple gravity model, and compare the output to SPICE generated
# trajectories.  The scenarios can be run with the followings setups
# parameters:
# Setup | scCase
# ----- | -------------------
# 1     | Hubble
# 2     | New Horizon
#
# To run the default scenario 1, call the python script through
#
#       python test_scenarioOrbitMultiBody.py
#
# When the simulation completes 2-3 plots are shown for each case.  One plot always shows
# the inertial position vector components, while the third plot shows the inertial differences
# between the Basilisk simulation trajectory and the SPICE spacecraft trajectory.  Read
# [test_scenarioBasicOrbit.py](@ref scenarioBasicOrbit) to learn how to setup an
# orbit simulation.
#
# The simulation layout is shown in the following illustration.  The SPICE interface object keeps track of
# the selection celestial objects, and ensures the gravity body object has the correct locations at
# each time step.
# ![Simulation Flow Diagram](Images/doc/test_scenarioOrbitMultiBody.svg "Illustration")
#
# The spacecraftPlus() module is setup as before, except that we need to specify a priority to this task.
# ~~~~~~~~~~~~~~~~~{.py}
#     # initialize spacecraftPlus object and set properties
#     scObject = spacecraftPlus.SpacecraftPlus()
#     scObject.ModelTag = "spacecraftBody"
#     scObject.hub.useTranslation = True
#     scObject.hub.useRotation = False
#
#     # add spacecraftPlus object to the simulation process
#     scSim.AddModelToTask(simTaskName, scObject, None, 1)
# ~~~~~~~~~~~~~~~~~
# If BSK modules are added to the simulation task process, they are executed in the order that they are added
# However, we the execution order needs to be control, a priority can be assigned.  The model with a higher priority
# number is executed first.  Modules with unset priorities will be given a priority of -1 which
# puts them at the
# very end of the execution frame.  They will get executed in the order in which they were added.
# For this scenario scripts, it is critical that the Spice object task is evaluated
# before the spacecraftPlus() model.  Thus, below the Spice object is added with a higher priority task.
#
# The first step to create a fresh gravity body factor class through
# ~~~~~~~~~~~~~~~~~{.py}
#   gravFactory = simIncludeGravBody.gravBodyFactory()
# ~~~~~~~~~~~~~~~~~
# This clears out the list of gravitational bodies, especially ifthe script is
# run multiple times using 'py.test' or in Monte-Carlo runs.
# Next a series of gravitational bodies are included.  Note that it is convenient to include them as a
# list of SPICE names.  The Earth is included in this scenario with the
# spherical harmonics turned on.  Note that this is true for both spacecraft simulations.
# ~~~~~~~~~~~~~~~~~{.py}
#    gravBodies = gravFactory.createBodies(['earth', 'mars barycenter', 'sun', 'moon', "jupiter barycenter"])
#    gravBodies['earth'].isCentralBody = True
#    gravBodies['earth'].useSphericalHarmParams = True
#    imIncludeGravBody.loadGravFromFile(bskPath + 'External/LocalGravData/GGM03S.txt'
#                                     , gravBodies['earth'].spherHarm
#                                     , 100
#                                     )
# ~~~~~~~~~~~~~~~~~
# The configured gravitational bodies are addes to the spacecraft dynamics with the usual command:
# ~~~~~~~~~~~~~~~~~{.py}
#    scObject.gravField.gravBodies = spacecraftPlus.GravBodyVector(gravFactory.gravBodies.values())
# ~~~~~~~~~~~~~~~~~
#
# Next, the default SPICE support module is created and configured.  The first step is to store
# the date and time of the start of the simulation.
# ~~~~~~~~~~~~~~~~~{.py}
#       timeInitString = "2012 MAY 1 00:28:30.0"
#       spiceTimeStringFormat = '%Y %B %d %H:%M:%S.%f'
#       timeInit = datetime.strptime(timeInitString,spiceTimeStringFormat)
# ~~~~~~~~~~~~~~~~~
# The following is a support macro that creates a `spiceObject` instance, and fills in typical
# default parameters.
# ~~~~~~~~~~~~~~~~~{.py}
#       gravFactory.createSpiceInterface(bskPath + 'External/EphemerisData/', timeInitString)
# ~~~~~~~~~~~~~~~~~
# Next the SPICE module is costumized.  The first step is to specify the zeroBase.  This is the inertial
# origin relative to which all spacecraft message states are taken.  The simulation defaults to all
# planet or spacecraft ephemeris being given in the SPICE object default frame, which is the solar system barycenter
# or SSB for short.  The spacecraftPlus() state output message is relative to this SBB frame by default.  To change
# this behavior, the zero based point must be redefined from SBB to another body.  In this simulation we use the Earth.
# ~~~~~~~~~~~~~~~~~{.py}
#   gravFactory.spiceObject.zeroBase = 'Earth'
# ~~~~~~~~~~~~~~~~~
# Finally, the SPICE object is added to the simulation task list.
# ~~~~~~~~~~~~~~~~~{.py}
#       scSim.AddModelToTask(simTaskName, gravFactory.spiceObject, None, -1)
# ~~~~~~~~~~~~~~~~~
# To unload the loaded SPICE kernels, use
# ~~~~~~~~~~~~~~~~~{.py}
# gravFactory.unloadSpiceKernels()
# ~~~~~~~~~~~~~~~~~
# This will unload all the kernels that the gravital body factory loaded earlier.
#
# Next we would like to import spacecraft specific SPICE ephemeris data into the python enviroment.  This is done
# such that the BSK computed trajectories can be compared in python with the equivalent SPICE directories.
# Note that this python SPICE setup is different from the BSK SPICE setup that was just completed.  As a result
# it is required to load in all the required SPICE kernels.  The following code is used to load either
# spacecraft data.
# ~~~~~~~~~~~~~~~~~{.py}
#     if scCase is 'NewHorizons':
#        scEphemerisFileName = 'nh_pred_od077.bsp'
#         scSpiceName = 'NEW HORIZONS'
#         vizPlanetName = "sun"
#     else:  # default case
#         scEphemerisFileName = 'hst_edited.bsp'
#         scSpiceName = 'HUBBLE SPACE TELESCOPE'
#         vizPlanetName = "earth"
#     pyswice.furnsh_c(gravFactory.spiceObject.SPICEDataPath + scEphemerisFileName)  # Hubble Space Telescope data
#     pyswice.furnsh_c(gravFactory.spiceObject.SPICEDataPath + 'de430.bsp')  # solar system bodies
#     pyswice.furnsh_c(gravFactory.spiceObject.SPICEDataPath + 'naif0011.tls')  # leap second file
#     pyswice.furnsh_c(gravFactory.spiceObject.SPICEDataPath + 'de-403-masses.tpc')  # solar system masses
#     pyswice.furnsh_c(gravFactory.spiceObject.SPICEDataPath + 'pck00010.tpc')  # generic Planetary Constants Kernel
# ~~~~~~~~~~~~~~~~~
# To unload the SPICE kernels loaded into the Python environment, use
# ~~~~~~~~~~~~~~~~~{.py}
#     pyswice.unload_c(gravFactory.spiceObject.SPICEDataPath + 'de430.bsp')  # solar system bodies
#     pyswice.unload_c(gravFactory.spiceObject.SPICEDataPath + 'naif0011.tls')  # leap second file
#     pyswice.unload_c(gravFactory.spiceObject.SPICEDataPath + 'de-403-masses.tpc')  # solar system masses
#     pyswice.unload_c(gravFactory.spiceObject.SPICEDataPath + 'pck00010.tpc')  # generic Planetary Constants Kernel
# ~~~~~~~~~~~~~~~~~
#
#
# The initial spacecraft position and velocity vector is obtained via the SPICE function call:
# ~~~~~~~~~~~~~~~~~{.py}
#       scInitialState = 1000*pyswice.spkRead(scSpiceName, timeInitString, 'J2000', 'EARTH')
#       rN = scInitialState[0:3]         # meters
#       vN = scInitialState[3:6]         # m/s
# ~~~~~~~~~~~~~~~~~
# Note that these vectors are given here relative to the Earth frame.  When we set the spacecraftPlus()
# initial position and velocity vectors through before initialization
# ~~~~~~~~~~~~~~~~~{.py}
#     scObject.hub.r_CN_NInit = unitTestSupport.np2EigenVectorXd(rN)  # m - r_CN_N
#     scObject.hub.v_CN_NInit = unitTestSupport.np2EigenVectorXd(vN)  # m - v_CN_N
# ~~~~~~~~~~~~~~~~~
# the natural question arises, how does Basilisk know relative to what frame these states are defined?  This is
# actually setup above where we set `.isCentralBody = True` and mark the Earth as are central body.
# Without this statement, the code would assume the spacecraftPlus() states are relative to the default zeroBase frame.
# In the earlier basic orbital motion script (@ref scenarioBasicOrbit) this subtleties were not discussed.
# This is because there
# the planets ephemeris message is being set to the default messages which zero's both the position and orientation
# states.  However, if Spice is used to setup the bodies, the zeroBase state must be carefully considered.
#
#
# Setup 1
# -----
#
# Which scenario is run is controlled at the bottom of the file in the code
# ~~~~~~~~~~~~~{.py}
# if __name__ == "__main__":
#     run( False,       # do unit tests
#          True,        # show_plots
#          'Hubble'
#        )
# ~~~~~~~~~~~~~
# The first 2 arguments can be left as is.  The remaining argument(s) control the
# simulation scenario flags to turn on or off certain simulation conditions.  The default
# scenario simulates the Hubble Space Telescope (HST) spacecraft about the Earth in a LEO orbit.
# The resulting position coordinates and orbit illustration are shown below.  A 2000 second simulation is
# performed, and the Basilisk and SPICE generated orbits match up very well.
# ![Inertial Position Coordinates History](Images/Scenarios/scenarioOrbitMultiBody1Hubble.svg "Position history")
# ![Perifocal Orbit Illustration](Images/Scenarios/scenarioOrbitMultiBody2Hubble.svg "Orbit Illustration")
# ![Trajectory Differences](Images/Scenarios/scenarioOrbitMultiBody3Hubble.svg "Trajectory Differences")
#
# Setup 2
# -----
#
# The next scenario is run by changing the bottom of the file in the scenario code to read
# ~~~~~~~~~~~~~{.py}
# if __name__ == "__main__":
#     run( False,       # do unit tests
#          True,        # show_plots
#          'NewHorizons'
#        )
# ~~~~~~~~~~~~~
# This case illustrates a simulation of the New Horizons spacecraft.  Here the craft is already a very
# large distance from the sun.  The
# resulting position coordinates and trajectorie differences are shown below.
# ![Inertial Position Coordinates History](Images/Scenarios/scenarioOrbitMultiBody1NewHorizons.svg "Position history")
# ![Trajectory Difference](Images/Scenarios/scenarioOrbitMultiBody3NewHorizons.svg "Trajectory Difference")
#
## @}
def main():
    spacecraft = np.array([[1000.,0.,0.,0.,0.,0.]])
    beacon     = np.array([[1200,1000,450]])

    ##################################################################################
    #
    # Camera/P&L Parameters
    #
    ##################################################################################
    
    extras = {}
    ## \var extras
    # Focal Length (mm)
    # 
    extras['focal_length'] = 100.

    # Camera resolution (pixels)
    extras['resolution'] = [2048., 512.]
    
    # width and height of pixels in camera
    extras['pixel_width']  = 2.5
    extras['pixel_height'] = 10.
 
    # direction coefficient of pixel and line axes
    extras['pixel_direction'] = 1.
    extras['line_direction']  = 1.

    extras['obs_beacons'] = [1.]

    angles = np.array([[0,np.pi / 4,np.pi / 2.]])
    angles = np.array([[0,        0,np.pi / 2.]])

    args = ( spacecraft, beacon, angles, extras )

    G = G_function( args )

    H = H_function( args )

    symbolic_results = pixel_and_line_symbolics.main()
    
    G_diff = G - np.array(symbolic_results[0:2])

    H_diff = H[:,0:3] - np.array([symbolic_results[2:5],symbolic_results[5:8]])

    if np.any( np.greater( np.abs( G_diff ), 10**(-10) ) ):
       print 'P&L G Function did not pass unit test :('
    else:
       print 'P&L G Function passed unit test!'

    if np.any( np.greater( np.abs( H_diff ), 10**(-10) ) ):
       print 'P&L H Function did not pass unit test :('
    else:
       print 'P&L H Function passed unit test!'

    return


if __name__ == "__main__":
    main()


