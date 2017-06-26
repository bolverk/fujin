def main():

    import numpy
    rawd = numpy.loadtxt('res.txt')
    return abs((rawd[0]-rawd[1])/rawd[0])<1e-3;

if __name__=='__main__':

    import os

    os.system('rm -rf test_passed.res test_failed.res')
    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
