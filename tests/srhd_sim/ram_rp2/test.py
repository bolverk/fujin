def velocity2celerity(v):

    import math

    return v/math.sqrt(1-v**2)

def calc_l1(a1, a2):

    import numpy
    return sum(numpy.abs(numpy.array(a1)-numpy.array(a2)))/len(a1)

def consolidate(fname):

    import h5py
    import numpy

    with h5py.File(fname) as f:
        res = dict((key,numpy.array(f[key])) for key in f)
    return res

def main():

    import numpy
    import os
    import sys
    sys.path.append(os.environ['FUJIN_ROOT']+'/analytic')
    from ideal_gas_riemann_solver import RiemannProfile, Primitive
    final = consolidate('final.h5')

    exact_prof = RiemannProfile(Primitive(1.,1000.,0.),
                                Primitive(1.,1e-2,0.),5./3.)
    ss_cord = (final['position']-0.5)/final['time'][0]
    exact = {'density': [exact_prof.calcPrimitive(ssc).d for ssc in ss_cord],
             'pressure': [exact_prof.calcPrimitive(ssc).p for ssc in ss_cord],
             'celerity': [exact_prof.calcPrimitive(ssc).w for ssc in ss_cord]}

    if False:
        import matplotlib
        matplotlib.use('Qt4Agg')
        import pylab
        for n, field in enumerate(['density','pressure','celerity']):
            pylab.subplot(3,1,n+1)
            pylab.plot(final['position'], final[field])
            pylab.plot(final['position'], exact[field])
        pylab.show()

    l1_data = dict((key,calc_l1(exact[key],final[key])) for key in exact)
    return (l1_data['density']<0.66 and
            l1_data['pressure']<36 and
            l1_data['celerity']<0.33)

if __name__=='__main__':
    
    import os

    os.system('rm -rf test_passed.res test_failed.res')
    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
