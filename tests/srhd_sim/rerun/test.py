#! /usr/bin/python

def calc_l1(a1, a2):

    import numpy
    return sum(numpy.abs(numpy.array(a1)-numpy.array(a2)))/len(a1)

def consolidate(fname):

    import h5py
    import numpy

    with h5py.File(fname,'r') as f:
        data = dict((key,numpy.array(f[key])) for key in f)
    return data

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

    import numpy
    import os
    import sys
    assert('FUJIN_ROOT' in os.environ)
    sys.path.append(os.environ['FUJIN_ROOT']+'/analytic')
    from ideal_gas_riemann_solver import RiemannProfile, Primitive
    from glob import glob

    continuous = consolidate('continuous.h5');
    interrupted = consolidate('interrupted.h5')
    checkpoint = consolidate('checkpoint.h5')

    for field in ['edges', 'density', 'pressure', 'celerity']:
        if calc_l1(continuous[field],
                   interrupted[field])>1e-12:
            return False    

    return True

if __name__=='__main__':
    
    import os

    os.system('rm -rf test_passed.res test_failed.res')
    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
