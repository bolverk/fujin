#! /usr/bin/python

def approx_equal(x1, x2, tol):
    """
    Checkes whether two values agree, within a given error margin
    """

    return abs((x1-x2)/(x1+x2))<tol

def l1_norm(a1, a2):

    import numpy

    return sum(numpy.abs(numpy.array(a1)-numpy.array(a2))/len(a1))

def calc_celerity(v):

    import math

    return v/math.sqrt(1-v**2)

def consolidate(fname):

    import h5py
    import numpy

    f = h5py.File(fname)
    res = dict((key,numpy.array(f[key])) for key in f)
    f.close()
    return res

def main():

    import sys
    sys.path.append(os.environ['FUJIN_ROOT']+'/analytic')
    from ideal_gas_riemann_solver import RiemannProfile, Primitive

    #import numpy
    #rawd = numpy.loadtxt('res.txt')
    #p_f = rawd[0]
    #v_f = rawd[1]
    #c_f = calc_celerity(v_f)
    #p_dm = 20.7445
    #v_dm = 0.956038
    #c_dm = calc_celerity(v_dm)

    #f = open('gradesheet.txt','w')
    #f.write('Numeric pressure = '+str(p_f)+'\n')
    #f.write('Numeric celerity = '+str(c_f)+'\n')
    #f.write('Exact pressure = '+str(p_dm)+'\n')
    #f.write('Exact celerity = '+str(c_dm)+'\n')
    #f.write('Relative pressure difference = '+str((p_f-p_dm)/(p_f+p_dm))+'\n')
    #f.write('Relative celerity difference = '+str((c_f-c_dm)/(c_f+v_dm))+'\n')
    #f.close()

    final = consolidate('final.h5')
    exact_prof = RiemannProfile(
        Primitive(1.,1000.,0),
        Primitive(1.,0.1,0), 5./3.)
    ss_cord = (final['position']-0.5)/(final['time'][0])
    analytic = {}
    analytic['density'] = [exact_prof.calcPrimitive(ssc).d for ssc in ss_cord]
    analytic['pressure'] = [exact_prof.calcPrimitive(ssc).p for ssc in ss_cord]
    analytic['celerity'] = [exact_prof.calcPrimitive(ssc).w for ssc in ss_cord]

    l1_data = dict((key,l1_norm(analytic[key],final[key])) for key in analytic)

    #import pylab
    #for n,key in enumerate(['density','pressure','celerity']):
    #    pylab.subplot(3,1,n+1)
    #    pylab.plot(final['position'], final[key])
    #    pylab.plot(final['position'], analytic[key])
    #    pylab.ylabel(key)
    #pylab.xlabel('Position')
    #pylab.subplot(311)
    #pylab.title('time = '+str(final['time'][0]))
    #pylab.show()

    f = open('gradesheet.txt','w')
    for key in ['density','pressure','celerity']:
        f.write(str(l1_data[key])+'\n')
    f.close()

    return l1_data['pressure'] < 22 and l1_data['celerity'] < 0.37 and l1_data['density']<0.83

if __name__=='__main__':

    import os

    os.system('rm -rf test_passed.res test_failed.res')
    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')

