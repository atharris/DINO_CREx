#!/usr/local/bin/python
'''

 <<Description>>


 <<Summary>>

'''

__author__ = 'Marc Balducci'
__version__ = '$Revision$'[11:-2]
__date__ = '$Date$'[7:26]

################################################################################
#                     I M P O R T     L I B R A R I E S
################################################################################
##############################################################################
# Log path in order to get pyswice from BSK

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


import numpy as np



################################################################################
#                  E X P O R T E D     F U N C T I O N S
################################################################################

# -------------------------------------------------------------------------------

def biasedTime(input):
    t_span   = input[0]
    extras   = input[1]

    # value of the last updated and corrected time. effectively, this is
    # t_0 in the bias equation
    t_0      = extras['last_t_update'] # s
    delta_F  = extras['measured_frequency_offset'] # hz
    dot_F    = extras['frequency_aging'] # hz/s
    allan    = extras['allan_parameters'] # tau (s), deviation (N/A)   

    t_length = t_span.shape[0]
    
    


    return dState


################################################################################
#                       M A I N     F U N C T I O N:
################################################################################

def main():
    return


if __name__ == "__main__":
    main()
