## \file pixel_and_line_unit_test.py
# This code initiates and runs a unit test for the following functions:
# Pixel and Line G function (measurement generation)
# Pixel and Line H function (measurement partials wrt estimate state)
#

import numpy as np
import pdb
import pixel_and_line_symbolics
from pixelLineBatch import fncH as H_function
from pixelLineBatch import fncG as G_function


##
# I'm making a list
# - here we go
#     - deeper
#         - even deeper  
#         .
#     should be sub list
#     .
# should be list
# .
# should not be list
#
def main():
    spacecraft = np.array([[1000.,0.,0.,0.,0.,0.]])
    beacon     = np.array([[1200,1000,450]])

    ##################################################################################
    #
    # Camera/P&L Parameters
    #
    ##################################################################################
    
    extras = {}
    # Focal Length (mm)
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

