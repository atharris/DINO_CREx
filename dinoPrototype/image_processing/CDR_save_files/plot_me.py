#! /usr/bin/env python3
# H+
#	Title   : camera_demo.py
#	Author  : Matt Muszynski
#	Date    : 03/23/17
#	Synopsis: 
#
#
#	$Date$
#	$Source$
#  @(#) $Revision$
#	$Locker$
#
#	Revisions:
#
# H-
# U+
#	Usage   :
#	Example	:
#	Output  :
# U-
# D+
#
# D-
###############################################################################

import numpy as np
import pdb
import matplotlib.pyplot as plt

_0_deg = np.load('0_deg.npz')
_45_deg = np.load('45_deg.npz')
_90_deg = np.load('90_deg.npz')
_135_deg = np.load('135_deg.npz')
stars_only = np.load('stars_only.npz')

plt.figure()
plt.imshow(_0_deg['detector_array'].reshape(512,512))
plt.figure()
plt.imshow(_45_deg['detector_array'].reshape(512,512))
plt.figure()
plt.imshow(_90_deg['detector_array'].reshape(512,512))
plt.figure()
plt.imshow(_135_deg['detector_array'].reshape(512,512))
plt.figure()
plt.imshow(stars_only['detector_array'].reshape(512,512))

##########
pdb.set_trace()
