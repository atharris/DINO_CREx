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
import numpy
import math
import batchFilter
import matplotlib.pyplot as plt

# input = []
#
# G = functions.fncG(input)
# H = functions.fncH(input)

a = numpy.array([1.,0.])

input = [1., 1., 1.]

normtest = batchFilter.norm(input)
tol = 1E-10
testFailCount = 0

if numpy.abs(normtest - numpy.sqrt(3)) > tol:
    testFailCount+=1
    print 'Test Failed'
else:
    print 'Test Passed'


input2 = [-1., 0., 0.]

normtest2 = batchFilter.norm(input2)


if numpy.abs(normtest2 - 1) > tol:
    testFailCount+=1
    print 'Test Failed'
else:
    print 'Test Passed'

input3 = [-100000., 0.0001, numpy.sqrt(2)]

normtest3 = batchFilter.norm(input3)

normtrue = numpy.linalg.norm(input3)

if numpy.abs(normtest3 - normtrue) > tol:
    testFailCount+=1
    print 'Test Failed'
else:
    print 'Test Passed'

print 'Test fail count =', testFailCount