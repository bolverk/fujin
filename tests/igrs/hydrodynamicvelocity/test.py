#! /usr/bin/python

def main():

    import numpy
    res = numpy.loadtxt('res.txt')
    return res[0]==-0.899119 and res[1]==0.901362

if __name__=='__main__':

    import os

    os.system('rm -rf test_passed.res test_failed.res')
    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')

