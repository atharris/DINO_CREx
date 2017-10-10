#! /usr/bin/env python3
import pdb
from numpy import array, arange, pi, \
average, log10, exp, sqrt, set_printoptions, get_printoptions
import matplotlib.pyplot as plt
import sys

def pix_per_sec(FOV,res,slew_rate):
	return slew_rate/(FOV/res)

def plot_it(res):
	plt.figure()
	slew_rate = arange(0.01,1.01,0.01)
	for FOV in range(2,11,2):
		label = str(FOV) + ' Degree FOV'
		plt.plot(slew_rate,pix_per_sec(FOV,res,slew_rate),label=label)

	plt.xlabel('Slew Rate (deg/sec)')
	plt.ylabel('Pixels/sec')
	plt.title('Pixels/Second, ' + str(res) + 'x' + str(res) + ' Detector')
	plt.ylim([0,1000])
	plt.legend()
	plt.figure()
	slew_rate = arange(0.05,1,0.1)
	for FOV in range(2,11,2):
		label = str(FOV) + ' Degree FOV'
		plt.plot(slew_rate,0.1/pix_per_sec(FOV,res,slew_rate),label=label)

	plt.xlabel('Slew Rate (deg/sec)')
	plt.ylabel('Time Step')
	plt.title('Timestep for <0.1 Pixel Motion, ' + str(res) + 'x' + str(res) + ' Detector')
	plt.ylim([0,.02])
	plt.legend()

def table_it(res):
	slew_rate = arange(0.25,2.1,0.25)
	print(res)
	print(slew_rate)
	for FOV in range(2,11,2):
		print(FOV, (0.1/pix_per_sec(FOV,res,slew_rate)))

set_printoptions(linewidth=200,formatter={'float': lambda x: format(x, '6.2e')})
table_it(512)
table_it(1024)
table_it(2048)
plt.show()
pdb.set_trace()