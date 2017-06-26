#! /usr/bin/python

import matplotlib
matplotlib.use('Qt4Agg')
import pylab
import numpy

fname = "plot.txt"
rawd = numpy.loadtxt(fname)
r = rawd[:,0]
d = rawd[:,1]
p = rawd[:,2]
v = rawd[:,3]
e = rawd[:,4]

fname = "plot2.txt"
rawd = numpy.loadtxt(fname)
r2 = rawd[:,0]
d2 = rawd[:,1]
p2 = rawd[:,2]
v2 = rawd[:,3]
e2 = rawd[:,4]

pylab.subplot(221)
pylab.plot(r,d, r2, d2)
pylab.xlabel('Radius')
pylab.ylabel('Density')

pylab.subplot(222)
pylab.plot(r,p,r2,p2)
pylab.xlabel('Radius')
pylab.ylabel('Pressure')

pylab.subplot(223)
pylab.plot(r,v,r2,v2,)
pylab.xlabel('Radius')
pylab.ylabel('Velocity')

pylab.subplot(224)
pylab.plot(r,e,r2,e2)
pylab.xlabel('Radius')
pylab.ylabel('Energy')

pylab.show()

pylab.clf

fname = "plot3.txt"
rawd = numpy.loadtxt(fname)
rv = rawd[:,0]
ps = rawd[:,1]
vs = rawd[:,2]

fname = "plot4.txt"
rawd = numpy.loadtxt(fname)
rv2 = rawd[:,0]
ps2 = rawd[:,1]
vs2 = rawd[:,2]

pylab.subplot(121)
pylab.plot(rv,ps,rv2,ps2)
pylab.xlabel('Radius')
pylab.ylabel('Pressure')

pylab.subplot(122)
pylab.plot(rv,vs,rv2,vs2)
pylab.xlabel('Radius')
pylab.ylabel('Velocity')

pylab.show()
