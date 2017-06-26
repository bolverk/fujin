#! /usr/bin/python

def consolidate(fname):

    import h5py
    import numpy

    with h5py.File(fname) as f:
        return dict((key,numpy.array(f[key])) for key in f)

def all_approx_equal(a1, a2, thres):

    for x,y in zip(a1,a2):
        if (abs(x)-abs(y))/(abs(x)+abs(y))>thres:
            return False
    return True

def main():

    import numpy

    data = consolidate('results.h5')
    return all_approx_equal(data['pressure'], data['reconstructed'], 1e-2)

if __name__=='__main__':

    import os

    os.system('rm -rf test_passed.res test_failed.res')
    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
