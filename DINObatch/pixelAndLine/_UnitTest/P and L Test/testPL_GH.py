
import sys, os, inspect

filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))
dinoName = 'DINO_CREx'
splitPath = path.split(dinoName)
dinoCommonPath = splitPath[0] + dinoName + '/DINObatch/pixelAndLine/commonFunctions/'
sys.path.append(dinoCommonPath)

import numpy as np
import pdb
import pixel_and_line_symbolics
from pixelLineBatch import fncH
from pixelLineBatch import fncG 


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
    angles = np.array([[0,        0,np.pi / 2.]])

    args = ( spacecraft, beacon, angles, extras )

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


