#! /usr/bin/python

def main():

    import numpy
    rawd = numpy.loadtxt('res.txt')
    ans = [1.02062,0.285774,1.33089]
    res = 1
    for i in range(len(rawd)):
        temp = abs((rawd[i]-ans[i])/ans[i])<0.01
        res = res*temp
    return res

if __name__=='__main__':

    import os

    os.system('rm -rf test_passed.res test_failed.res')
    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
