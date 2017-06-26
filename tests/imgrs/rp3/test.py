#! /usr/bin/python

def approx_equal(x1, x2, tol):
    """
    Checkes whether two values agree, within a given error margin
    """

    return abs((x1-x2)/(x1+x2))<tol

def main():

    import numpy
    rawd = numpy.loadtxt('res.txt')
    p = rawd[0]
    v = rawd[1]
    pa = 20.7445
    va = 0.956038
    return approx_equal(p,pa,0.01) and \
        approx_equal(v,va,0.01)

if __name__=='__main__':
    
    import os

    os.system('rm -rf test_passed.res test_failed.res')
    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')

