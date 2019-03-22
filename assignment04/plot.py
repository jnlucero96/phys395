#!/usr/bin/env python

# import libraries
import matplotlib
import sys

try:
	import astropy.io.fits as pyfits
except ImportError:
	import pyfits

import numpy as np

from pylab import *
from scipy.signal import medfilt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm

def question1(data_file):

	# load trajectory data
	hdu = pyfits.open(data_file)

	data = hdu[0].data

	ny,nx = data.shape

	# antialias using median filter
	osy = 1; osx = 1

	if (osy != 1 or osx != 1):
		aakernel = [(osy & (~1)) + 1, (osx & (~1)) + 1]
		data = medfilt(data, aakernel)[::osy,::osx]; ny,nx = data.shape
		print("Antialiased with median kernel %s; output size is (%i,%i) pixels" % (aakernel,nx,ny))

	# data extent and meshes
	x0 = hdu[0].header['CRVAL1']
	y0 = hdu[0].header['CRVAL2']
	dx = hdu[0].header['CDELT1']
	dy = hdu[0].header['CDELT2']

	x = x0 + arange(nx)*dx*osx
	y = y0 + arange(ny)*dy*osy

	X, Y = meshgrid(x, y)

	hdu.close()

	abs_data_max = data.__abs__().max()

	# create the figure
	fig = figure(figsize=(10,8), frameon=False)

	# plot image data
	c = matplotlib.colors.LinearSegmentedColormap.from_list("difference", ["blue", "white", "red"])
	im = imshow(
		data,
		extent=(x[0],x[-1],y[0],y[-1]), origin='lower',
		aspect=(osx*dx)/(osy*dy), cmap=cm.get_cmap('afmhot'), interpolation='none',
		vmin=-abs_data_max, vmax=abs_data_max # make the bounds symmetric
		)

	# plot contours
	#contour(X, Y, data, 32, cmap=cm.jet)

	# make colorbar match the plot
	# divider = make_axes_locatable(plt.gca())
	# cax = divider.append_axes("right", size="5%", pad=0.07)
	cax = plt.gca()
	fig.colorbar(im, ax=cax)

	fig.savefig("question1.pdf")

def question2(data_file):

	# load trajectory data
	t, pend2_traj, energy_deviation = np.loadtxt(
		data_file, usecols=(0,2,3), unpack=True
		)

	# create the figure
	fig, ax = subplots(2, 1, figsize=(10,8), sharex=True)

	# plot image data
	ax[0].plot(t, pend2_traj, lw=3.0)
	ax[0].tick_params(labelsize=18)
	ax[0].set_ylabel("Pend. 2 Pos.", fontsize=20)

	ax[1].plot(t, energy_deviation, lw=3.0)
	ax[1].tick_params(labelsize=18)
	ax[1].yaxis.offsetText.set_fontsize(18)
	ax[1].yaxis.get_major_formatter().set_powerlimits((0,0))
	ax[1].set_ylabel("Energy deviation", fontsize=20)
	ax[1].set_xlabel("Time", fontsize=20)

	for i in range(2):
		ax[i].set_xlim([t.min(), t.max()])

	fig.tight_layout()
	fig.savefig("question2.pdf")

def question3(data_file):

	# load trajectory data
	hdu = pyfits.open(data_file)

	data = hdu[0].data

	ny,nx = data.shape

	# antialias using median filter
	osy = 1; osx = 1

	if (osy != 1 or osx != 1):
		aakernel = [(osy & (~1)) + 1, (osx & (~1)) + 1]
		data = medfilt(data, aakernel)[::osy,::osx]; ny,nx = data.shape
		print("Antialiased with median kernel %s; output size is (%i,%i) pixels" % (aakernel,nx,ny))

	# data extent and meshes
	x0 = hdu[0].header['CRVAL1']
	y0 = hdu[0].header['CRVAL2']
	dx = hdu[0].header['CDELT1']
	dy = hdu[0].header['CDELT2']

	x = x0 + arange(nx)*dx*osx
	y = y0 + arange(ny)*dy*osy

	X, Y = meshgrid(x, y)

	hdu.close()

	# calculate the max of the data
	abs_data_max = data.__abs__().max()

	fig = figure(figsize=(10,8), frameon=False)
	c = matplotlib.colors.LinearSegmentedColormap.from_list("difference", ["blue", "white", "red"])
	im = imshow(
		data,
		extent=(x[0],x[-1],y[0],y[-1]), origin='lower',
		aspect=(osx*dx)/(osy*dy), cmap=cm.get_cmap('afmhot'), interpolation='none',
		vmin=0.0, vmax=abs_data_max
		)

	# plot contours
	#contour(X, Y, data, 32, cmap=cm.jet)

	# make colorbar match the plot
	# divider = make_axes_locatable(plt.gca())
	# cax = divider.append_axes("right", size="5%", pad=0.07)
	cax = plt.gca()
	cax.grid(True)
	fig.colorbar(im, ax=cax)

	fig.savefig("question3.pdf")

if __name__ == "__main__":
	data_file = sys.argv[1]
	problem = sys.argv[2]

	if (int(problem) == 1):
		question1(data_file)
	elif (int(problem) == 2):
		question2(data_file)
	elif (int(problem) == 3):
		question3(data_file)
	else:
		"Input for problem cannot be understood."
