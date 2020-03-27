#!/usr/bin/env python

from numpy import loadtxt
from matplotlib import rcParams, rc, ticker, colors, cm
from matplotlib.style import use
from matplotlib.pyplot import subplots, close
from matplotlib.cm import get_cmap

use('seaborn-paper')
rc('text', usetex=True)
rcParams['mathtext.fontset'] = 'cm'
rcParams['text.latex.preamble'] = [
    r"\usepackage{amsmath}", r"\usepackage{lmodern}",
    r"\usepackage{siunitx}", r"\usepackage{units}",
    r"\usepackage{physics}", r"\usepackage{bm}"
]

# read out time, position, velocity, and energy deviation from data file
t, x, v, delta = loadtxt("p3_data.dat", unpack=True)

# ready the figure
fig, ax = subplots(3, 1, figsize=(10,12), sharex=True)

# plot positions
ax[0].plot(t, x, lw=2.0, color='black')
ax[0].set_ylabel(r"$x$",fontsize=24)
ax[0].set_title("Position", fontsize=30)

#plot velocities
ax[1].plot(t, v, lw=2.0, color='red')
ax[1].set_ylabel(r"$\dot{x}$", fontsize=24)
ax[1].set_title("Velocity", fontsize=30)

# plot energy deviations
ax[2].plot(t, delta, lw=2.0)
ax[2].set_ylabel(r"$\delta$", fontsize=24)
ax[2].set_title("Energy Deviation", fontsize=30)
ax[2].set_xlabel(r"$t$", fontsize=24)


for axis in ax: 
    axis.tick_params(labelsize=14)
    axis.grid(True)

fig.tight_layout()
fig.savefig("q3_solution.pdf")