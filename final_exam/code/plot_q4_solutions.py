#!/usr/bin/env python

from pylab import *

initE, period = loadtxt("p4_data.dat", unpack=True)

semilogy(initE, period, lw=3.0)
xlabel("Initial Energy")
ylabel("Period of Oscillation")

tight_layout()
savefig("q4_figure.pdf")