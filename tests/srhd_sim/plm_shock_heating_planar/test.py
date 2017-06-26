#! /usr/bin/python

def main():

    import numpy
    ai = 5./3.
    v = 1.0-1e-10
    lf = (1.0-v**2)**-0.5
    d = 1.0
    p = (lf-1.0)*(lf*ai-1.0)*d
    res = numpy.loadtxt('res.txt')
    ans = p
    return abs(res-ans)/ans<0.01

if __name__=='__main__':
    
    import os

    os.system('rm -rf test_passed.res test_failed.res')
    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
