#! /usr/bin/python

def main():

    import numpy
    res = numpy.loadtxt('res.txt')
    pre = 0.146165
    vel = 0.089294
    return (res[0]-pre)/pre<0.01 and \
        (res[1]-vel)/vel<0.01

if __name__=='__main__':
    
    import os

    os.system('rm -rf test_passed.res test_failed.res')
    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
