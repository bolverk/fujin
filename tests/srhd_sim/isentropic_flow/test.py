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

def load_snapshot(fname):

    import h5py
    import numpy

    with h5py.File(fname,'r') as f:
        return {field:numpy.array(f[field])
                for field in f}

def consolidate(pattern):

    from glob import glob
    import re
    import numpy

    file_list = sorted(glob(pattern),
                       key=lambda fname:int(re.search('(\d+)',fname)[0]))
    partitions = [load_snapshot(fname) for fname in file_list]
    return {field:numpy.concatenate([part[field] for part in partitions])
            for field in partitions[0]}

def main():

    import numpy
    import math
    from glob import glob

    if len(glob('initial_*.h5'))>0:
        initial = consolidate('initial_*.h5')
        final = consolidate('final_*.h5')
    else:
        initial = load_snapshot('initial.h5')
        final = load_snapshot('final.h5')
    x0 = initial['position']
    d0 = initial['density']
    p0 = initial['pressure']
    v0 = initial['celerity']/numpy.sqrt(initial['celerity']**2+1)
    xn = final['position']
    dn = final['density']
    pn = final['pressure']
    vn = final['celerity']/numpy.sqrt(final['celerity']**2+1)

    g = 5./3.
    t = 0.8
    e0 = [p/(g-1)+d for d,p in zip(d0,p0)]
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
