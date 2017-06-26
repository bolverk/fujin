#! /usr/bin/python

import matplotlib
matplotlib.use('Qt4Agg')
import pylab
import numpy
import sys

def calc_celerity(v):

    import math

    return v/math.sqrt(1-v**2)

aux = sys.argv[0]
pref =  aux[0:aux.find('plotit.py')]

fname = "plot.txt"
rawd = numpy.loadtxt(pref+fname);
r_list = rawd[:,0]
d_list = rawd[:,1]
p_list = rawd[:,2]
v_list = rawd[:,3]
c_list = [calc_celerity(v) for v in v_list]

pylab.subplot(311)
pylab.plot(r_list,d_list)
pylab.ylabel('Density')

pylab.subplot(312)
pylab.plot(r_list,p_list)
pylab.ylabel('Pressure')

pylab.subplot(313)
pylab.plot(r_list,v_list)
pylab.xlabel('Radius')
pylab.ylabel('Celerity')

pylab.show()
