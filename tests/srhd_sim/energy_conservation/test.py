#! /usr/bin/python

def all_the_same(a):

    import numpy

    num = numpy.diff(a)
    den = 0.5*(a[1:]+a[:-1])

    return (numpy.abs(num/den)<1e-3).all()

def consolidate():

    import numpy
    from glob import glob
    import re
    import logging

    LOGLEVEL = os.environ.get('LOGLEVEL', 'WARNING').upper()
    logging.basicConfig(level=LOGLEVEL)

    file_list = glob('res_*.txt')
    file_list = sorted(file_list,
                       key=lambda fname:int(re.search('(\d+)',fname)[0]))
    partitions = numpy.array([numpy.loadtxt(fname) for fname in file_list])
    logging.debug(partitions)
    logging.debug(numpy.sum(partitions,axis=0))
    res = numpy.sum(partitions,axis=0)
    return res

def main():

    import numpy
    from glob import glob

    if len(glob('res_*.txt'))>1:
        res = consolidate()
    else:
        res = numpy.loadtxt('res.txt')
    return all_the_same(res)

if __name__=='__main__':
    import os

    os.system('rm -rf test_passed.res test_failed.res')
    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
