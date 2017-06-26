#! /usr/bin/python

def main():

    import numpy
    res = numpy.loadtxt('res.txt')
    ans = 0.302898
    return abs(res-ans)/ans<0.01

if __name__=='__main__':
    
    import os

    os.system('rm -rf test_passed.res test_failed.res')
    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
