#! /usr/bin/python

def approx_equal(x1, x2, tol):
    """
    Checkes whether two values agree, within a given error margin
    """

    return abs((x1-x2)/(x1+x2))<tol

def main():

    import numpy
    rawd = numpy.loadtxt('res.txt')
    p1 = rawd[0]
    v1 = rawd[1]
    p2 = rawd[2]
    v2 = rawd[3]
    return p1==p2 and v1==-v2

if __name__=='__main__':
    
    import os

    os.system('rm -rf test_passed.res test_failed.res')
    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')

