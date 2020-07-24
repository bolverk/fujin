#! /usr/bin/python

def velocity2celerity(v):

    import math

    return v/math.sqrt(1-v**2)

def calc_l1(a1, a2):

    import numpy
    return sum(numpy.abs(numpy.array(a1)-numpy.array(a2)))/len(a1)

def consolidate(fname):

    import h5py
    import numpy

    with h5py.File(fname,'r') as f:
        return {field:numpy.array(f[field]) for field in f}

def consolidate_all(pattern):

    from glob import glob
    import re
    import numpy
    import logging

    LOGLEVEL = os.environ.get('LOGLEVEL', 'WARNING').upper()
    logging.basicConfig(level=LOGLEVEL)

    file_list = sorted(glob(pattern),
                       key=lambda fname:int(re.search('(\d+)',fname)[0]))
    logging.debug(file_list)
    partitions = [consolidate(fname) for fname in file_list]
    return {field:numpy.concatenate([part[field] for part in partitions])
            for field in partitions[0]}

def main():

    import numpy
    import os
    import sys
    from glob import glob
    assert('FUJIN_ROOT' in os.environ)
    sys.path.append(os.environ['FUJIN_ROOT']+'/analytic')
    from ideal_gas_riemann_solver import RiemannProfile, Primitive

    if len(glob('final_*.h5'))>1:
        final = consolidate_all('final_*.h5')
    else:
        final = consolidate('final.h5')
    exact_prof = RiemannProfile(Primitive(1,1,velocity2celerity(0.9)),
                                Primitive(1,10,0), 4./3.)
    ss_cord = (final['position']-0.5)/final['time'][0]
    exact = {'density': [exact_prof.calcPrimitive(ssc).d for ssc in ss_cord],
             'pressure': [exact_prof.calcPrimitive(ssc).p for ssc in ss_cord],
             'celerity': [exact_prof.calcPrimitive(ssc).w for ssc in ss_cord]}

    if False:
        import matplotlib
        matplotlib.use('Qt4Agg')
        import pylab
        for n, field in enumerate(['density','pressure','celerity']):
            pylab.subplot(3,1,n+1)
            pylab.plot(final['position'], final[field])
            pylab.plot(final['position'], exact[field])
        pylab.show()

    l1_data = dict((key,calc_l1(exact[key],final[key])) for key in exact)

    return (l1_data['density']<0.1 and
            l1_data['pressure']<0.22 and
            l1_data['celerity']<0.011)

if __name__=='__main__':

    import os

    os.system('rm -rf test_passed.res test_failed.res')
    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
