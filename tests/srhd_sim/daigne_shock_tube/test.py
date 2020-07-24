#! /usr/bin/python

def approx_equal(x1, x2, tol):
    """
    Checkes whether two values agree, within a given error margin
    """

    return abs((x1-x2)/(x1+x2))<tol

def l1_norm(a1, a2):

    import numpy

    return sum(numpy.abs(numpy.array(a1)-numpy.array(a2))/len(a1))

def calc_celerity(v):

    import math

    return v/math.sqrt(1-v**2)

def consolidate(fname):

    import h5py
    import numpy

    f = h5py.File(fname)
    res = dict((key,numpy.array(f[key])) for key in f)
    f.close()
    return res

def consolidate_all(pattern):

    from glob import glob
    import re
    import numpy

    file_list = sorted(glob(pattern),
                       key=lambda fname:int(re.search('(\d+)',fname)[0]))
    partitions = [consolidate(fname) for fname in file_list]
    return {field:numpy.concatenate([part[field] for part in partitions])
            for field in partitions[0]}

def main():

    import sys
    sys.path.append(os.environ['FUJIN_ROOT']+'/analytic')
    from ideal_gas_riemann_solver import RiemannProfile, Primitive
    from glob import glob

    if len(glob('final_*.h5'))>1:
        final = consolidate_all('final_*.h5')
    else:
        final = consolidate('final.h5')
    exact_prof = RiemannProfile(
        Primitive(1.,1000.,0),
        Primitive(1.,0.1,0), 5./3.)
    ss_cord = (final['position']-0.5)/(final['time'][0])
    analytic = {}
    analytic['density'] = [exact_prof.calcPrimitive(ssc).d for ssc in ss_cord]
    analytic['pressure'] = [exact_prof.calcPrimitive(ssc).p for ssc in ss_cord]
    analytic['celerity'] = [exact_prof.calcPrimitive(ssc).w for ssc in ss_cord]

    l1_data = dict((key,l1_norm(analytic[key],final[key])) for key in analytic)

    return l1_data['pressure'] < 22 and l1_data['celerity'] < 0.37 and l1_data['density']<0.83

if __name__=='__main__':

    import os

    os.system('rm -rf test_passed.res test_failed.res')
    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')

