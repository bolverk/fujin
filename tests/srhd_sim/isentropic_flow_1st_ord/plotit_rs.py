#! /usr/bin/python

import matplotlib
matplotlib.use('Qt4Agg')
import pylab
import numpy
import sys

aux = sys.argv[0]
pref = aux[:aux.find('plotit_rs.py')]

fname = "plot_rs.txt"
rawd = numpy.loadtxt(pref+fname);
r = rawd[:,0]
p = rawd[:,1]
v = rawd[:,2]

pylab.subplot(211)
pylab.plot(r,p,'x')
pylab.xlabel('Radius')
pylab.ylabel('Pressure')

pylab.subplot(212)
pylab.plot(r,v,'x')
pylab.xlabel('Radius')
pylab.ylabel('Velocity')

pylab.show()
