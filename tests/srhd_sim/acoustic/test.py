def L1_error_norm(a1, a2):

    abs_diff = [abs(x-y) for x,y in zip(a1,a2)]
    return sum(abs_diff)

import math

mean_density = 1
mean_pressure = 1e-3
pert_density = 1e-5
g = 4./3.
wavelength = 1
k = 2*math.pi/wavelength
mean_energy = mean_pressure/(g-1)+mean_density
mean_enthalpy = mean_energy+mean_pressure
sound_speed = math.sqrt(g*mean_pressure/mean_enthalpy)
pert_pressure = g*mean_pressure*pert_density/mean_density
pert_velocity = sound_speed*pert_density/mean_density

def initial_density(x):

    return mean_density+pert_density*math.sin(k*x)

def initial_pressure(x):

    return mean_pressure+pert_pressure*math.sin(k*x)

def initial_velocity(x):

    return pert_velocity*math.sin(k*x)

def load_snapshot(fname):

    import h5py
    import numpy

    with h5py.File(fname,'r') as f:
        return {field:numpy.array(f[field])
                for field in f}

def consolidate(pattern):

    from glob import glob
    import re
    import numpy

    file_list = sorted(glob(pattern),
                       key=lambda fname:int(re.search('(\d+)',fname)[0]))
    partitions = [load_snapshot(fname) for fname in file_list]
    return {field:numpy.concatenate([part[field] for part in partitions])
            for field in partitions[0]}

def main():

    import numpy
    from glob import glob

    if len(glob('initial_*.h5')):
        initial = consolidate('initial_*.h5')
        final = consolidate('final_*.h5')
    else:
        initial = load_snapshot('initial.h5')
        final = load_snapshot('final.h5')
    x_init = initial['position']
    d_init = initial['density']
    p_init = initial['pressure']
    v_init = initial['celerity']/numpy.sqrt(initial['celerity']**2+1)
    x_final = final['position']
    d_final = final['density']
    p_final = final['pressure']
    v_final = final['celerity']/numpy.sqrt(final['celerity']**2+1)
    time = final['time'][0]

    d_exact = [initial_density(x-sound_speed*time) for x in x_final]
    p_exact = [initial_pressure(x-sound_speed*time) for x in x_final]
    v_exact = [initial_velocity(x-sound_speed*time) for x in x_final]

    if False:

        import matplotlib
        matplotlib.use('Qt4Agg')
        import pylab

        pylab.subplot(311)
        pylab.plot(x_init,d_init)
        pylab.plot(x_init,[initial_density(x) for x in x_init])
        pylab.plot(x_final,d_final)
        pylab.plot(x_final,d_exact)

        pylab.subplot(312)
        pylab.plot(x_init,p_init)
        pylab.plot(x_init,[initial_pressure(x) for x in x_init])
        pylab.plot(x_final,p_final)
        pylab.plot(x_final,p_exact)

        pylab.subplot(313)
        pylab.plot(x_init,v_init)
        pylab.plot(x_init,[initial_velocity(x) for x in x_init])
        pylab.plot(x_final,v_final)
        pylab.plot(x_final,v_exact)

        pylab.legend(('initial numeric',
                      'initial exact',
                      'final numeric',
                      'final exact'))

        pylab.show()

    if False:

        import pylab
        import math

        entropy = [p/d**g for d,p in zip(d_final,p_final)]
        c_final = [math.sqrt(g*p/d) for d,p in zip(d_final,p_final)]
        riemann_right = [v+2*c/(g-1) for v,c in zip(v_final,c_final)]
        riemann_left = [v-2*c/(g-1) for v,c in zip(v_final,c_final)]
        pylab.subplot(211)
        pylab.plot(x_final,entropy)

        pylab.subplot(212)
        pylab.plot(x_final,riemann_right)
        pylab.plot(x_final,riemann_left)

        pylab.show()
    
    l1_vals = dict(density=L1_error_norm(d_exact,d_final),
                   pressure=L1_error_norm(p_exact,p_final),
                   velocity=L1_error_norm(v_exact,v_final))

    f = open('gradesheet.txt','w')
    for i in l1_vals:
        f.write(str(l1_vals[i])+'\n')
    f.close()

    res = True
    expected_l1 = {'density':1e-4,
                   'pressure':1e-6,
                   'velocity':1e-5};

    f = open('gradesheet.txt','w')
    for i in l1_vals:
        f.write(i+' '+str(l1_vals[i])+'\n')
    f.close()

    for i in l1_vals:
        res = res and l1_vals[i]<expected_l1[i]
    
    return res

if __name__=='__main__':

    import os

    os.system('rm -rf test_passed.res test_failed.res')
    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
