#! /usr/bin/python

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
    from glob import glob
    import logging

    LOGLEVEL = os.environ.get('LOGLEVEL', 'WARNING').upper()
    logging.basicConfig(level=LOGLEVEL)

    ai = 4./3.
    v = 0.99999
    lf = (1.0-v**2)**-0.5
    d = 1.0
    p = (lf-1.0)*(lf*ai-1.0)*d
    df = d*((ai*lf+1)/(ai-1))**3/lf**2
    p = (ai-1.0)*(lf-1)*df
    ans = p
    
    if len(glob('final_*.h5'))>1:
        final = consolidate('final_*.h5')
    else:
        final = load_snapshot('final.h5')
    index = int(len(final['position'])/20)
    res = final['pressure'][index]

    return abs(res-ans)/ans<0.1

if __name__=='__main__':

    import os

    os.system('rm -rf test_passed.res test_failed.res')
    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
