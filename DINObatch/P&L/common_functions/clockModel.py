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

import numpy as np



################################################################################
#                  E X P O R T E D     F U N C T I O N S
################################################################################

# -------------------------------------------------------------------------------

def referenceTime(input):
    t_span   = input[0]
    extras   = input[-1]

    # value of the last updated and corrected time. effectively, this is
    # t_0 in the bias equation
    t_0      = extras['last_t_update'] # s
    F_meas   = extras['measured_frequency'] # hz
    F_0      = extras['reference_frequency'] # hz
    dot_F    = extras['frequency_aging'] # hz/s
    allan    = extras['allan_parameters'] # [ tau (s), deviation (N/A) ] 

    # length of the time span in question
    t_length = t_span.shape[0]
 
    # difference between measured frequency and reference frequency
    delta_F  = F_meas - F_0

    # preallocate variable for the time bias at each measurement time
    t_bias   = np.zeros( (t_length,) )
 
    # loop through each measurement to populate the bias vector
    for tt in xrange(t_length-1):
      # difference in time between current observation and last 
      delta_t   = t_span[tt+1] - t_span[tt]
      allan_idx = np.argmax( allan[:,0] > delta_t )
      allan_t   = 

    return dState

def simulatedTime(input):
    t_span   = input[0]
    extras   = input[-1]

    # value of the last updated and corrected time. effectively, this is
    # t_0 in the bias equation
    t_0      = extras['last_t_update'] # s
    F_meas   = extras['measured_frequency'] # hz
    F_0      = extras['reference_frequency'] # hz
    dot_F    = extras['frequency_aging'] # hz/s
    allan    = extras['allan_parameters'] # [ tau (s), deviation (N/A) ] 

    # length of the time span in question
    t_length = t_span.shape[0]
 
    # difference between measured frequency and reference frequency
    delta_F  = F_meas - F_0

    # preallocate variable for the time bias at each measurement time
    t_bias   = np.zeros( (t_length,) )
 
    # loop through each measurement to populate the bias vector
    for tt in xrange(t_length-1):
      # difference in time between current observation and last 
      delta_t   = t_span[tt+1] - t_span[tt]
      allan_idx = np.argmax( allan[:,0] > delta_t )
      allan_t   = 

    return dState

def initiateTimeSim(input):

    extras   = input[-1]

    # initial values for some parameters. These values are to be randomized
    # in order to provide a more "true" behavior for the simulation. The initial
    # values will act as means for the randomizer 
    F_0      = extras['reference_frequency'] # hz
    dot_F    = extras['frequency_aging'] # hz/s
    allan    = extras['allan_parameters'] # [ tau (s), deviation (N/A) ] 

    # length of the time span in question
    t_length = t_span.shape[0]
 
    # difference between measured frequency and reference frequency
    delta_F  = F_meas - F_0

    # preallocate variable for the time bias at each measurement time
    t_bias   = np.zeros( (t_length,) )
 
    # loop through each measurement to populate the bias vector
    for tt in xrange(t_length-1):
      # difference in time between current observation and last 
      delta_t   = t_span[tt+1] - t_span[tt]
      allan_idx = np.argmax( allan[:,0] > delta_t )
      allan_t   = 

    return dState
################################################################################
#                       M A I N     F U N C T I O N:
################################################################################

def main():
    return


if __name__ == "__main__":
    main()
