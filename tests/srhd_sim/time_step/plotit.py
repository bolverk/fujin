#! /usr/bin/python

import matplotlib
matplotlib.use('Qt4Agg')
import pylab
import numpy

fname = "plot.txt"
rawd = numpy.loadtxt(fname);
pylab.plot(rawd[:,0],rawd[:,1])
pylab.xlabel('Radius')
pylab.ylabel('Velocity')
pylab.show()
