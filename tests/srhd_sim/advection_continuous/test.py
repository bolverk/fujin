def L1_error_norm(a1, a2):

    abs_diff = [abs(x-y) for x,y in zip(a1,a2)]
    return sum(abs_diff)

def main():

    import numpy

    x_init, d_init, p_init, v_init = numpy.loadtxt('init_cond.txt',
                                                   unpack=True)
    x_final, d_final, p_final, v_final = numpy.loadtxt('snapshot.txt',
                                                       unpack=True)
    
    l1_vals = dict(density=L1_error_norm(d_init,d_final),
                   pressure=L1_error_norm(p_init,p_final),
                   velocity=L1_error_norm(v_init,v_final))

    f = open('gradesheet.txt','w')
    for i in l1_vals:
        f.write(str(l1_vals[i])+'\n')
    f.close()

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
