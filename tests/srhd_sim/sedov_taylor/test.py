#! /usr/bin/python

def consolidate():

    import logging
    import numpy
    from glob import glob

    LOGLEVEL = os.environ.get('LOGLEVEL', 'WARNING').upper()
    logging.basicConfig(level=LOGLEVEL)

    partitions = [numpy.loadtxt(fname) 
                  for fname in glob('rpmax_*.txt')]
    combined_r = numpy.zeros(len(partitions[0].T[0]))
    combined_t = numpy.zeros(len(partitions[0].T[0]))
    combined_p = numpy.zeros(len(partitions[0].T[0]))
    for n in range(len(partitions[0].T[0])):
        p_list = numpy.array([part.T[2][n] for part in partitions])
        r_list = numpy.array([part.T[1][n] for part in partitions])
        t_list = numpy.array([part.T[0][n] for part in partitions])
        combined_p[n] = numpy.max(p_list)
        combined_r[n] = r_list[numpy.argmax(p_list)]
        combined_t[n] = t_list[0]
    return numpy.array([combined_t,combined_r,combined_p]).T

def main():

    import numpy
    from glob import glob
    import logging

    LOGLEVEL = os.environ.get('LOGLEVEL', 'WARNING').upper()
    logging.basicConfig(level=LOGLEVEL)

    if len(glob('rpmax_*.txt'))>1:
        rawd = consolidate()
    else:
        rawd = numpy.loadtxt('rpmax.txt')
    logging.debug(rawd)
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
