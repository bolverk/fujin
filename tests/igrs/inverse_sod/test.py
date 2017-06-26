#! /usr/bin/python

def main():

    import numpy
    res = numpy.loadtxt('res.txt')
    return res[0]<1e-15 and abs(res[1]-res[2])/(res[1]+res[2])<0.01 and abs(res[2]-res[3])/(res[2]+res[3])<0.01

if __name__=='__main__':

    import os

    os.system('rm -rf test_passed.res test_failed.res')
    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
