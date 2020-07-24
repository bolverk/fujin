def L1_error_norm(a1, a2):

    abs_diff = [abs(x-y) for x,y in zip(a1,a2)]
    return sum(abs_diff)

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
    import logging

    LOGLEVEL = os.environ.get('LOGLEVEL', 'WARNING').upper()
    logging.basicConfig(level=LOGLEVEL)
    
    if len(glob('initial_*.h5'))>1:
        initial = consolidate('initial_*.h5')
        final = consolidate('final_*.h5')
    else:
        initial = load_snapshot('initial.h5')
        final = load_snapshot('final.h5')
    x_init = initial['position']
    d_init = initial['density']
    p_init = initial['pressure']
    v_init = initial['celerity']
    x_final = initial['position']
    d_final = initial['density']
    p_final = initial['pressure']
    v_final = initial['celerity']

    l1_vals = dict(density=L1_error_norm(d_init,d_final),
                   pressure=L1_error_norm(p_init,p_final),
                   velocity=L1_error_norm(v_init,v_final))

    logging.debug(l1_vals)

    res = True
    for i in l1_vals:
        res = res and l1_vals[i]==0
    
    return res

if __name__=='__main__':

    import os

    os.system('rm -rf test_passed.res test_failed.res')
    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
