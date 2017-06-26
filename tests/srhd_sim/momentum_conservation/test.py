#! /usr/bin/python

def all_the_same(a):
    """
    Verifies that all the terms in an arra a are the same
    """

    res = True
    for i in a:
        res = res and (i==a[0])
    return res

def main():

    import numpy
    res = numpy.loadtxt('res.txt')
    return all_the_same(res)

if __name__=='__main__':

    import os

    os.system('rm -rf test_passed.res test_failed.res')
    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')

