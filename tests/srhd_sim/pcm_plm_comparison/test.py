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

def celerity2velocity(w):

    import numpy

    return w/numpy.sqrt(w**2+1)

def main():

    import numpy
    import logging
    from glob import glob

    LOGLEVEL = os.environ.get('LOGLEVEL', 'WARNING').upper()
    logging.basicConfig(level=LOGLEVEL)

    if len(glob('pcm_*.h5')):
        pcm = consolidate('pcm_*.h5')
        plm = consolidate('plm_*.h5')
    else:
        pcm = load_snapshot('pcm.h5')
        plm = load_snapshot('plm.h5')
    index = int(len(pcm['position'])/2)
    rawd = [pcm['pressure'][index],
             celerity2velocity(pcm['celerity'][index]),
             plm['pressure'][index],
             celerity2velocity(plm['celerity'][index])]          
    logging.debug(rawd)

    return abs(rawd[0]-rawd[2])/rawd[2]<0.01 and\
        abs(rawd[1]-rawd[3])/rawd[3]<0.3

if __name__=='__main__':
    
    import os

    os.system('rm -rf test_passed.res test_failed.res')
    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
