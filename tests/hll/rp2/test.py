#! /usr/bin/python

def main():

    import numpy
    res = numpy.loadtxt('res.txt')
    ans = [0.146165, 0.089294]
    return abs(res[0]-ans[0])/ans[0]<0.1 and \
        abs(res[1]-ans[1])/ans[1]<0.2

if __name__=='__main__':

    import os

    os.system('rm -rf test_passed.res test_failed.res')
    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
