#! /usr/bin/python

def main():

    import numpy

    rawd = numpy.loadtxt('rpmax.txt')
    istart = 1000
    t = rawd[istart:,0]
    r = rawd[istart:,1]
    fitd = numpy.polyfit(numpy.log(t), numpy.log(r),1)
    sdpl = 0.4;
    return abs(fitd[0]-sdpl)/sdpl < 0.01

if __name__=='__main__':

    import os

    os.system('rm -rf test_passed.res test_failed.res')
    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
