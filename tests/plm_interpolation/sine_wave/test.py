def main():

    if False:

        import pylab
        import numpy
        
        edges = numpy.loadtxt('edges.txt')
        pcm_data = numpy.loadtxt('pcm_data.txt')
        plm_data = numpy.loadtxt('plm_data.txt')
        analytic_data = numpy.loadtxt('analytic_data.txt')

        print(len(edges))
        print(len(pcm_data))
        pylab.plot(edges,pcm_data)
        pylab.plot(edges,plm_data)
        pylab.plot(edges,analytic_data)
        pylab.show()

    import numpy

    points = numpy.loadtxt('points.txt')
    l1_pcm = numpy.loadtxt('l1_pcm.txt')
    l1_plm = numpy.loadtxt('l1_plm.txt')
    pcm_fitd = numpy.polyfit(numpy.log(points),
                             numpy.log(l1_pcm),1)
    plm_fitd = numpy.polyfit(numpy.log(points),
                             numpy.log(l1_plm),1)

    if False:

        import pylab

        pylab.loglog(points,l1_pcm)
        pylab.loglog(points,l1_plm)
        pylab.show()

    return abs(pcm_fitd[0]+1)<1e-2 and \
        abs(plm_fitd[0]+2)<2e-2

if __name__=='__main__':

    import os

    os.system('rm -rf test_passed.res test_failed.res')
    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')

