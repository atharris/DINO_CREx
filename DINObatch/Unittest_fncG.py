''' '''
'''
 ISC License

 Copyright (c) 2016-2017, Autonomous Vehicle Systems Lab, University of Colorado at Boulder

 Permission to use, copy, modify, and/or distribute this software for any
 purpose with or without fee is hereby granted, provided that the above
 copyright notice and this permission notice appear in all copies.

 THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

'''
import pytest
import sys, os, inspect

#
# Spice Unit Test
#
# Purpose:  Test the proper function of the Spice Ephemeris module.
#           Proper function is tested by comparing Spice Ephermis to
#           JPL Horizons Database for different planets and times of year
# Author:   Thibaud Teil
# Creation Date:  Dec. 20, 2016
#

filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))

import rngRngRtBatch as functions
import numpy as np
import math
import batchFilter
import matplotlib.pyplot as plt

testFailCount = 0

SCState = np.array([0.,0.,0.,0.,0.,0.])

BeaconStates = np.zeros([2,6])
BeaconStates[0,:] = np.array([10.,20.,15.,1.,2.,3.])
BeaconStates[1,:] = np.array([100.,2.,-15.,5.,1.,2.])

extras = {}
extras['obs_beacons'] = BeaconStates[:,0]
print len(BeaconStates[:,0])

input = [SCState, BeaconStates, extras]


G = functions.fncG(input)
# H = functions.fncH(input)

print G

print 'Test fail count =', testFailCount