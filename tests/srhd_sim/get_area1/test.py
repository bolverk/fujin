#! /usr/bin/python

def main():

    import numpy
    rawd = numpy.loadtxt('res.txt')
    res = 1
    n = len(rawd[:,0])
    for i in range(n):
        if rawd[i,0]!=rawd[i,1]:
            res = 0
    return res

if __name__=='__main__':

    import os

    os.system('rm -rf test_passed.res test_failed.res')
    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
