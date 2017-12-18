
import sys, os, inspect

filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))
dinoName = 'DINO_CREx'
splitPath = path.split(dinoName)
dinoCommonPath = splitPath[0] + dinoName + '/DINObatch/P&L/common_functions/'
sys.path.append(dinoCommonPath)

import numpy as np
import pdb
import pixel_and_line_symbolics
from pixelLineBatch import fncH
from pixelLineBatch import fncG 

## \defgroup test_PL testPL_GH - unit test for fncG() and fncH() of pixelLineBatch.py
##   @{
## The unit test for exportable functions of `pixelLineBatch.py`.
#
# Overview {#overview}
# ====
#
# Purpose
# -----
# This script runs a unit test of the `fncG()` `fncH()` functons found within \ref pixel_and_line . Differences between two methodologies are computed, and the functions pass if the error is below a defined threshold. 
#
# Contents
# -----
# The `testPL_GH.py` script is run to execute a single main function:
#
# - `main`
#
# This main function defines parameters which are then passed to two different methods of computing measurements (G) and the observation-state mapping matrix (H). 
#
# The Code
# =====
#
# `main`
# -----
# Because `testPL_GH.py` is a script meant to be run from the terminal, there are no inputs. The operator may choose to tweak various mission settings or spacecraft parameters, many of which are found in \ref init_vanilla .
#
# The first code in the script creates a single state for a spacecraft and a beacon
# ~~~~~~~~~~~~~~~~{.py}
#    spacecraft = np.array([[1000.,0.,0.,0.,0.,0.]])
#    beacon     = np.array([[1200,1000,450]])
# ~~~~~~~~~~~~~~~~
#
# These states consitute the position and velocity of the spacecraft and beacon at the time of measurement. The next lines then go on to define some parameters for the camera, which are stored in the `extras` dictionary found on many other batch filter functions. 
#
# After these definitions, the `fncG()` and `fncH()` functions of \ref pixel_and_line are run
# ~~~~~~~~~~~~~~~~{.py}
#    args = ( spacecraft, beacon, angles, extras )
#
#    G = fncG( args )
#
#    H = fncH( args )
# ~~~~~~~~~~~~~~~~
#
# The next step is to calculate G and H using symbolic functions that have been coded separately in \ref pixel_and_line_symbolic . At the time of this documentation, the main function of \ref pixel_and_line_symbolic uses the same inputs as this script, but there are no inputs for that of the symbolic function. Therefore, if an operator wishes to change some mission paramters for the test, they must then edit \ref pixel_and_line_symobolic , or develop a input support. This symbolic call is found at
# ~~~~~~~~~~~~~~~~{.py}
#    symbolic_results = pixel_and_line_symbolics.main()
# ~~~~~~~~~~~~~~~~
#
# Differences between the outputs of the two methods are then made, and if any absolute difference is greater than 10^-10, the test fails.
# ~~~~~~~~~~~~~~~~{.py}
#    diffG = G - np.array(symbolic_results[0:2])
#
#    diffH = H[:,0:3] - np.array([symbolic_results[2:5],symbolic_results[5:8]])
#
#    if np.any( np.greater( np.abs( diffG ), 10**(-10) ) ):
#       print 'P&L G Function did not pass unit test :('
#    else:
#       print 'P&L G Function passed unit test!'
#
#    if np.any( np.greater( np.abs( diffH ), 10**(-10) ) ):
#       print 'P&L H Function did not pass unit test :('
#    else:
#       print 'P&L H Function passed unit test!'
# ~~~~~~~~~~~~~~~~ 
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
    extras['FoL'] = 100.

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
    angles = np.array([0,        0,np.pi / 2.])

    args = ( spacecraft, beacon, angles, extras )
    pdb.set_trace()
    G = fncG( args )

    H = fncH( args )

    symbolic_results = pixel_and_line_symbolics.main()
    
    diffG = G - np.array(symbolic_results[0:2])

    diffH = H[:,0:3] - np.array([symbolic_results[2:5],symbolic_results[5:8]])

    if np.any( np.greater( np.abs( diffG ), 10**(-10) ) ):
       print 'P&L G Function did not pass unit test :('
    else:
       print 'P&L G Function passed unit test!'

    if np.any( np.greater( np.abs( diffH ), 10**(-10) ) ):
       print 'P&L H Function did not pass unit test :('
    else:
       print 'P&L H Function passed unit test!'

    return


if __name__ == "__main__":
    main()


