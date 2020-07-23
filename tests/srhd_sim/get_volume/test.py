#! /usr/bin/python

def verify_file(fname):

    import numpy
    rawd = numpy.loadtxt(fname)
    return ((rawd.T[0]-rawd.T[1])==0).all()

def main():

    from glob import glob

    for fname in glob('res*.txt'):
        if not verify_file(fname):
            return False
    return True

if __name__=='__main__':

    import os

    os.system('rm -rf test_passed.res test_failed.res')
    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
