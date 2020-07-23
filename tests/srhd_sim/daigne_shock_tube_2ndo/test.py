def approx_equal(x1, x2, tol):
    """
    Checkes whether two values agree, within a given error margin
    """

    return abs((x1-x2)/(x1+x2))<tol

def calc_celerity(v):

    import math

    return v/math.sqrt(1-v**2)

def load_file(fname):

    import numpy
    import h5py

    with h5py.File(fname,'r') as f:
        return {field:numpy.array(f[field]) for field in f}

def consolidate_files():

    import logging
    LOGLEVEL = os.environ.get('LOGLEVEL', 'WARNING').upper()
    logging.basicConfig(level=LOGLEVEL)
    from glob import glob
    import re
    import numpy

    file_list = glob('./final_*.h5')
    file_list = sorted(file_list,
                       key=lambda fname:int(re.search('(\d+)',fname)[0]))
    logging.debug(file_list)
    partitions = [load_file(fname) for fname in file_list]
    return {field:numpy.concatenate([part[field] for part in partitions])
            for field in partitions[0]}

def main():

    from glob import glob
    from scipy import interpolate
    import numpy
    import logging

    LOGLEVEL = os.environ.get('LOGLEVEL', 'WARNING').upper()
    logging.basicConfig(level=LOGLEVEL)

    if len(glob('final_*.h5'))>0:
        data = consolidate_files()
    else:
        data = load_file('final.h5')

    index = int(len(data['position'])/2)+1
    p_f = data['pressure'][index]
    c_f = data['celerity'][index]
    p_dm = 20.7445
    v_dm = 0.956038
    c_dm = calc_celerity(v_dm)
    logging.debug([p_f, p_dm])
    logging.debug([c_f, c_dm])

    return approx_equal(p_f, p_dm, 0.02) and \
        approx_equal(c_f, c_dm, 0.01)

if __name__=='__main__':

    import os

    os.system('rm -rf test_passed.res test_failed.res')
    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
