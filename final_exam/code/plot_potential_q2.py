#!/usr/bin/env python
from pylab import *

x = linspace(-4, 4, 100)
f = 4*x**3 + 9*x**2 - 8*x -3
axhline(0.0, ls='--', color='k', lw=2.0)
grid(True)

plot(x, f, lw=3.0)
savefig("plot_potential_q2_figure.pdf")