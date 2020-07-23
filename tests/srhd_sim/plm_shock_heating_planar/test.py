#! /usr/bin/python

def load_file(fname):

    import numpy
    import h5py

    with h5py.File(fname,'r') as f:
        return {field:numpy.array(f[field]) for field in f}

def consolidate_files():

    import logging
    LOGLEVEL = os.environ.get('LOGLEVEL', 'WARNING').upper()
    logging.basicConfig(level=LOGLEVEL)
    from glob import glob
    import re
    import numpy

    file_list = glob('./final_*.h5')
    file_list = sorted(file_list,
                       key=lambda fname:int(re.search('(\d+)',fname)[0]))
    logging.debug(file_list)
    partitions = [load_file(fname) for fname in file_list]
    return {field:numpy.concatenate([part[field] for part in partitions])
            for field in partitions[0]}

def main():

    import numpy
    from scipy import interpolate
    import logging
    import os
    from glob import glob
    LOGLEVEL = os.environ.get('LOGLEVEL', 'WARNING').upper()
    logging.basicConfig(level=LOGLEVEL)

    ai = 5./3.
    v = 1.0-1e-10
    lf = (1.0-v**2)**-0.5
    d = 1.0
    p = (lf-1.0)*(lf*ai-1.0)*d
    ans = p

    if len(glob('final_*.h5')) > 0:
        data = consolidate_files()
    else:
        data = load_file('final.h5')
    logging.debug([field for field in data])
    res = interpolate.interp1d(data['position'],
                               data['pressure'])(0.5)
    logging.debug(ans)
    logging.debug(res)

    return abs(res-ans)/ans<0.01

if __name__=='__main__':
    
    import os

    os.system('rm -rf test_passed.res test_failed.res')
    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
