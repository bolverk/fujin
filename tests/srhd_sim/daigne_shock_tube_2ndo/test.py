def approx_equal(x1, x2, tol):
    """
    Checkes whether two values agree, within a given error margin
    """

    return abs((x1-x2)/(x1+x2))<tol

def calc_celerity(v):

    import math

    return v/math.sqrt(1-v**2)

def main():

    import numpy
    rawd = numpy.loadtxt('res.txt')
    p_f = rawd[0]
    v_f = rawd[1]
    c_f = calc_celerity(v_f)
    p_dm = 20.7445
    v_dm = 0.956038
    c_dm = calc_celerity(v_dm)

    f = open('gradesheet.txt','w')
    f.write('Numerical pressure = '+str(p_f)+'\n')
    f.write('Numerical celerity = '+str(c_f)+'\n')
    f.write('Exact pressure = '+str(p_dm)+'\n')
    f.write('Exact celerity = '+str(c_dm)+'\n')
    f.write('Relative pressure difference = '+str((p_f-p_dm)/(p_f+p_dm))+'\n')
    f.write('Relative celerity difference = '+str((c_f-c_dm)/(c_f+c_dm))+'\n')
    f.close()

    return approx_equal(p_f, p_dm, 0.02) and \
        approx_equal(c_f, c_dm, 0.01)

if __name__=='__main__':

    import os

    os.system('rm -rf test_passed.res test_failed.res')
    if main():
        os.system('touch test_passed.res')
    else:
        os.system('touch test_failed.res')
