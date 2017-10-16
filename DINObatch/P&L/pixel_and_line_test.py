import numpy as np
import pdb
from pixelLineBatch import fncH as H_function
from pixelLineBatch import fncG as G_function

def main():
    spacecraft = np.array([[1000.,0.,0.,0.,0.,0.]])
    beacon     = np.array([[1200,1000,400]])

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

    args = ( spacecraft, beacon, angles, extras )

    G = G_function( args )

    print 'G', G[0,0], G[0,1]

    H = H_function( args )
    pdb.set_trace()
    print 'H', H[0,0], H[0,1]

    return


if __name__ == "__main__":
    main()

