#! /usr/bin/python

import matplotlib
matplotlib.use('Qt4Agg')
import pylab
import numpy
import sys
import re

pref = re.sub(r'/[^/]*$','',sys.argv[0])

fname = pref+"/plot.txt"
rawd = numpy.loadtxt(fname);
r = rawd[:,0]
d = rawd[:,1]
p = rawd[:,2]
v = rawd[:,3]

pylab.subplot(311)
pylab.plot(r,d)
pylab.xlabel('Radius')
pylab.ylabel('Density')

pylab.subplot(312)
pylab.plot(r,p)
pylab.xlabel('Radius')
pylab.ylabel('Pressure')

pylab.subplot(313)
pylab.plot(r,v)
pylab.xlabel('Radius')
pylab.ylabel('Velocity')

pylab.show()
