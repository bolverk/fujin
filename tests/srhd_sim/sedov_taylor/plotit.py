#! /usr/bin/python

import matplotlib
matplotlib.use('Qt4Agg')
import pylab
import numpy
import sys

pref = sys.argv[0].replace('plotit.py','')
fname = pref+"plot.txt"
rawd = numpy.loadtxt(fname);
r = rawd[:,0]
d = rawd[:,1]
p = rawd[:,2]
v = rawd[:,3]
e = rawd[:,4]

pylab.subplot(221)
pylab.plot(r,d)
pylab.xlabel('Radius')
pylab.ylabel('Density')

pylab.subplot(222)
pylab.plot(r,p)
pylab.xlabel('Radius')
pylab.ylabel('Pressure')

pylab.subplot(223)
pylab.plot(r,v)
pylab.xlabel('Radius')
pylab.ylabel('Velocity')

pylab.subplot(224)
pylab.plot(r,e)
pylab.xlabel('Radius')
pylab.ylabel('Energy')

pylab.show()
