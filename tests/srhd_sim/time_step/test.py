#! /usr/bin/python

def main():

    import numpy
    import logging
    from glob import glob

    LOGLEVEL = os.environ.get('LOGLEVEL', 'WARNING').upper()
    logging.basicConfig(level=LOGLEVEL)

    if len(glob('res_*.txt'))>1:
        logging.debug([numpy.loadtxt(fname) for fname in glob('res_*.txt')])
        res = numpy.min([numpy.loadtxt(fname) 
                         for fname in glob('res_*.txt')])
    else:
        res = numpy.loadtxt('res.txt')
    logging.debug(res)
    ans = 0.00874773
    return abs((res-ans)/ans)<0.01

if __name__=='__main__':

    import os

    os.system('rm -rf test_passed.res test_failed.res')
    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
