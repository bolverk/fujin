#! /usr/bin/python

def relveladd(v1, v2):

    return (v1+v2)/(1+v1*v2)

def get_vars(fname):
    
    import numpy
    
    rawd = numpy.loadtxt(fname)
    x = rawd[:,0]
    d = rawd[:,1]
    p = rawd[:,2]
    v = rawd[:,3]
    return x, d, p, v

def main():

    import numpy
    import math

    x0, d0, p0, v0 = get_vars('plot0.txt')
    xn, dn, pn, vn = get_vars('plot.txt')

    g = 5./3.
    e0 = [p/(g-1)+d for d,p in zip(d0,p0)]
    t = 0.8
    cs = [math.sqrt(g*p0[i]/(e0[i]+p0[i])) for i in range(len(x0))]
    vp = [relveladd(cs[i],v0[i]) for i in range(len(x0))]
    xa = [x0[i]+t*vp[i] for i in range(len(x0))]

    xaf = [xa[i] for i in range(len(xa)) if (xa[i]>0 and xa[i]<1)]
    p0f = [p0[i] for i in range(len(xa)) if (xa[i]>0 and xa[i]<1)]
    pnf = numpy.interp(xaf, xn, pn)

    res = sum(((p0f-pnf)/max(p0f))**2)/len(p0f)
	
    return res<0.01

if __name__=='__main__':

    import os

    os.system('rm -rf test_passed.res test_failed.res')
    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
