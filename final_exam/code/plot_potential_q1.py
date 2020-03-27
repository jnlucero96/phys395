#!/usr/bin/env python
from pylab import *

x = linspace(-10, 10, 100)
f = cos(x) - x/5.0
axhline(0.0, ls='--', color='k', lw=2.0)
grid(True)

plot(x, f, lw=3.0)
savefig("plot_potential_q1_figure.pdf")