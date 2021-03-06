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
    from glob import glob

    if len(glob('final_*.h5'))>1:
        final = consolidate('final_*.h5')
    else:
        final = load_snapshot('final.h5')

    index = int(len(final['position'])/2+1)
    res = [final['pressure'][index],
           celerity2velocity(final['celerity'][index])]
    pre = 0.146165
    vel = 0.089294
    return (res[0]-pre)/pre<0.01 and \
        (res[1]-vel)/vel<0.01

if __name__=='__main__':
    
    import os

    os.system('rm -rf test_passed.res test_failed.res')
    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
