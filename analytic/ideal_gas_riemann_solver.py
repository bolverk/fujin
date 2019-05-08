def calc_hugoniot_density(p,d0,p0,g):

    import math

    a = p/d0 + ((d0 + d0*g + g*p)*p0)/(d0**2*(-1 + g)) + \
        (g*p0**2)/(d0**2*(-1 + g)**2)
    b = -(((1 + g)*p)/(-1 + g)) - p0
    c = -((g*p**2)/(-1 + g)**2) - (g*p*p0)/(-1 + g)
    return (-b+math.sqrt(b**2-4*a*c))/(2*a)

def test_calc_hugoniot():

    import numpy
    import pylab

    p_list = numpy.exp(numpy.linspace(-12,12,1000))
    d0 = 1
    p0 = p_list[0]
    g = 4./3.
    d_list = [calc_hugoniot_density(p,d0,p0,g)
              for p in p_list]
    pylab.loglog(p_list, d_list)
    pylab.show()
              
def calc_hugoniot_density_dp(p, d0, p0, g):
    d = calc_hugoniot_density(p,d0,p0,g);
    return (d0*(d*(-1 + g)*(d + d0 - d*g + d0*g) + 2*d0*g*p) + \
                (-d**2 + d0**2)*(-1 + g)*g*p0)/\
                (-(d0**2*(-1 + g)*((1 + g)*p + (-1 + g)*p0)) + \
                      2*d*(g*p0*((-1 + g)*p + p0) + \
                               d0*(-1 + g)*((-1 + g)*p + p0 + g*p0)))

def test_calc_hugoniot_density_dp():

    import numpy
    import pylab

    p_list = numpy.exp(numpy.linspace(-12,12,1000))
    d0 = 1
    p0 = p_list[0]
    g = 4./3.
    d_list = [calc_hugoniot_density(p,d0,p0,g)
              for p in p_list]
    p_mid = [0.5*(p_list[i]+p_list[i+1])
             for i in range(len(p_list)-1)]
    dp_list = [calc_hugoniot_density_dp(p,d0,p0,g)
               for p in p_list]
    dp_numeric = numpy.diff(d_list)/numpy.diff(p_list)
    
    pylab.loglog(p_list, dp_list, '.')
    pylab.loglog(p_mid, dp_numeric, '.')
    pylab.show()

def calc_hugoniot_celerity(p, d0, p0, g):
    import math

    d = calc_hugoniot_density(p,d0,p0,g)
    numerator = (g-1)*(p-p0)*((g-1)*(d-d0)+p-p0)
    denominator = (g*p+d*(g-1))*(d0*(g-1)+g*p0)
    return math.sqrt(numerator/denominator)

def test_calc_hugoniot_celerity():

    import numpy
    import pylab

    p_list = numpy.exp(numpy.linspace(-12,12,1000))
    d0 = 1
    p0 = p_list[0]
    g = 4./3.
    w_list = [calc_hugoniot_celerity(p,d0,p0,g)
              for p in p_list]

    pylab.plot(p_list, w_list)
    pylab.show()

def calc_hugoniot_celerity_dp(p, d0, p0, g):
    d = calc_hugoniot_density(p,d0,p0,g)
    dddp = calc_hugoniot_density_dp(p,d0,p0,g)
    w = calc_hugoniot_celerity(p,d0,p0,g)
    dlnwdp = (d**2*(-1 + g)**2 - d*(-1 + g)*\
                  (d0*(-1 + g) - 2*p + 2*p0 - g*p0) + \
                  g*(p**2 - p0*(d0*(-1 + g) + p0)))/\
                  (2.*(d*(-1 + g) + g*p)*(p - p0)*\
                       (d0 + d*(-1 + g) - d0*g + p - p0))
    dlnwdd = ((-1 + g)*(d0*(-1 + g) + (-1 + g)*p + p0))/\
    (2.*(d*(-1 + g) + g*p)*(d0 + d*(-1 + g) - d0*g + p - p0))
    return w*(dlnwdp+dlnwdd*dddp)

def test_calc_hugoniot_celerity_dp():

    import numpy
    import pylab

    p_list = numpy.exp(numpy.linspace(-12,12,1000))
    d0 = 1
    p0 = p_list[0]
    g = 4./3.
    w_list = [calc_hugoniot_celerity(p,d0,p0,g)
              for p in p_list]
    p_mid = [0.5*(p_list[i]+p_list[i+1])
             for i in range(len(p_list)-1)]
    dw_list = [calc_hugoniot_celerity_dp(p,d0,p0,g)
               for p in p_list]
    dw_numeric = numpy.diff(w_list)/numpy.diff(p_list)
    
    pylab.loglog(p_list, dw_list, '.')
    pylab.loglog(p_mid, dw_numeric, '.')
    pylab.show()

def calc_isentrope_density(p,d0,p0,g):
    
    return d0*(float(p)/float(p0))**(1.0/g)

def test_calc_isentrope_density():

    import numpy
    import pylab

    d0 = 1
    p_list = numpy.exp(numpy.linspace(-12,12,1000))
    g = 4./3.
    p0 = numpy.max(p_list)
    d_list = [calc_isentrope_density(p,d0,p0,g)
              for p in p_list]
    
    pylab.loglog(p_list, d_list)
    pylab.show()

def calc_sound_speed(d,p,g):

    return (1.0/(g-1)+d/(g*p))**(-0.5)

def calc_isentrope_celerity(p, d0, p0, g):

    import math

    ba0 = calc_sound_speed(d0,p0,g)
    ba = calc_sound_speed(calc_isentrope_density(p,d0,p0,g),
                          p,g)
    return math.sinh((2.0/math.sqrt(g-1))*\
                         (math.atanh(ba/math.sqrt(g-1))-\
                              math.atanh(ba0/math.sqrt(g-1))))

def test_calc_isentrope_celerity():

    import numpy
    import pylab

    p_list = numpy.exp(numpy.linspace(-15,-3,1000))
    p0 = numpy.max(p_list)
    d0 = 1
    g = 4./3.
    w_list = [-calc_isentrope_celerity(p,d0,p0,g)
              for p in p_list]
    pylab.loglog(p_list, w_list)
    pylab.show()

def calc_isentrope_celerity_dp(p, d0, p0, g):

    import math

    ba0 = calc_sound_speed(d0,p0,g)
    ba = calc_sound_speed(calc_isentrope_density(p,d0,p0,g),
                          p, g)
    return (ba/(g*p))*\
        math.cosh((2.0/math.sqrt(g-1))*\
                      (math.atanh(ba/math.sqrt(g-1))-\
                           math.atanh(ba0/math.sqrt(g-1))))

def test_calc_isentrope_celerity_dp():

    import numpy
    import pylab

    p_list = numpy.exp(numpy.linspace(-5,-3,1000))
    d0 = 1
    g = 4./3.
    p0 = numpy.max(p_list)
    w_list = [calc_isentrope_celerity(p,d0,p0,g)
              for p in p_list]
    p_mid = [0.5*(p_list[i]+p_list[i+1])
             for i in range(len(p_list)-1)]
    dw_list = [calc_isentrope_celerity_dp(p,d0,p0,g)
               for p in p_list]
    dw_numeric = numpy.diff(w_list)/numpy.diff(p_list)
    
    pylab.plot(p_list, dw_list,'.')
    pylab.plot(p_mid, dw_numeric,'.')
    pylab.show()

def calc_hydrodynamic_celerity(p,d0,p0,g):

    if p>p0:
        return calc_hugoniot_celerity(p,d0,p0,g)
    else:
        return calc_isentrope_celerity(p,d0,p0,g)

def calc_hugoniot_celerity_dp(p,d0,p0,g):

    if p>p0:
        return calc_hugoniot_celerity_dp(p,d0,p0,g)
    else:
        return calc_isentrope_celerity_dp(p,d0,p0,g)

def celerity_addition(w1, w2):

    if(abs(w1)>1e10):
        raise NameError('blat')

    import math

    return w1*math.sqrt(1+w2**2)+w2*math.sqrt(1+w1**2)

def celerity_addition_diff(w1, w2, dw1, dw2):

    import math

    return (sqrt(w1**2+1)+w1*w2/sqrt(w2**2+1))*dw2+\
        (sqrt(w2**2+1)+w1*w2/sqrt(w1**2+1))*dw1;

def right_celerity(p,d0,p0,w0,g):

    return celerity_addition(w0,
                             calc_hydrodynamic_celerity(p,d0,p0,g))

def right_celerity_diff(p,d0,p0,w0,g):

    return celerity_addition_diff(w0,
                                  calc_hydrodynamic_celerity(p,d0,p0,g),
                                  0,
                                  calc_hydrodynamic_celerity_dp(p,d0,p0,g))

def left_celerity(p,d0,p0,w0,g):

    return celerity_addition(w0,
                             -calc_hydrodynamic_celerity(p,d0,p0,g))

def left_celerity_diff(p,d0,p0,w0,g):

    return celerity_addition_diff(w0,
                                  -calc_hydrodynamic_celerity(p,d0,p0,g),
                                  0,
                                  -calc_hydrodynamic_celerity_dp(p,d0,p0,g))

def eval_trans_eqn(p,dl,pl,wl,dr,pr,wr,g):

    return right_celerity(p,dr,pr,wr,g)-left_celerity(p,dl,pl,wl,g)

def eval_trans_eqn_diff(p,dl,pl,wl,dr,pr,wr,g):

    return right_celerity_diff(p,dr,pr,wr,g)-left_celerity_diff(p,dl,pl,wl,g)

def test_eval_trans_eqn():

    import numpy
    import pylab

    dl = 1
    pl = 1000
    wl = 0
    dr = 1
    pr = 1e-2
    wr = 0
    g = 5./3.
    p_list = numpy.linspace(10,20,1000)
    eqn_val = [eval_trans_eqn(p,dl,pl,wl,dr,pr,wr,g)
               for p in p_list]
    pylab.plot(p_list, eqn_val)
    pylab.plot(p_list,[0 for p in p_list],'k')
    pylab.show()

def overestimate_shock_pressure(w,d0,p0,g):

    import math

    ba0 = calc_sound_speed(d0,p0,g)
    return p0*(2*ba0**2+w*g*(w*g+math.sqrt(4*ba0**2+(w*g)**2)))/\
        (2*ba0**2)

def two_shock_bounds(dl,pl,wl,dr,pr,wr,g):

    pmin = min(pl,pr)
    dw = celerity_addition(wl,-wr)
    pmax = max(overestimate_shock_pressure(dw,dl,pl,g),
               overestimate_shock_pressure(dw,dr,pr,g))
    return pmin, pmax

def calc_isentrope_pressure(w,d0,p0,g):

    ba0 = calc_sound_speed(d0,p0,g)
    ba = math.sqrt(g-1)*\
        math.tanh(0.5*math.sqrt(g-1)*\
                      math.asinh(w)+math.atanh(ba0/sqrt(g-1)))
    return p0*((g*p0/d0)*(ba**-2-1.0/(g-1)))**(-g/(g-1))

def two_rarefaction_bounds(dl,pl,wl,dr,pr,wr,g):

    pmax = max(pl,pr)
    dw = celerity_addition(wl,-wr)
    pmin = min(calc_isentrope_pressure(dw,dl,pl,g),
               calc_isentrope_pressure(dw,dr,pr,g))
    return pmin, pmax

def solve_transcendental_equation(dl,pl,wl,dr,pr,wr,g):

    import scipy.optimize

    wrpl = right_celerity(pl,dr,pr,wr,g)
    wlpr = left_celerity(pr,dl,pl,wl,g)
    func = lambda p: eval_trans_eqn(p,dl,pl,wl,dr,pr,wr,g)
    fprime = lambda p: eval_trans_eqn_diff(p,dl,pl,wl,dr,pr,wr,g)
    if max(wrpl,wr)<min(wlpr,wl):
        pmin, pmax = two_shock_bounds(dl,pl,wl,dr,pr,wr,g)
    elif min(wrpl,wr)>max(wlpr,wl):
        pmin, pmax = two_rarefaction_bounds(dl,pl,wl,dr,pr,wr,g)
    else:
        pmin = min(pl,pr)
        pmax = max(pl,pr)
    func = lambda p: eval_trans_eqn(p,dl,pl,wl,dr,pr,wr,g)

    return scipy.optimize.brentq(func,pmin,pmax)

def riemann_solve(left, right, g):

    ps = solve_transcendental_equation(left.d,
                                       left.p,
                                       left.w,
                                       right.d,
                                       right.p,
                                       right.w,g)
    ws = 0.5*(right_celerity(ps,right.d,right.p,right.w,g)+\
                  left_celerity(ps,left.d,left.p,left.w,g))
    return ps, ws

def test_shock_rarefaction():

    dl = 1
    pl = 1000
    wl = 0
    dr = 1
    pr = 1e-2
    wr = 0
    g = 5./3.
    print(solve_transcendental_equation(dl,pl,wl,dr,pr,wr,g))

def calc_energy(d,p,g):

    return p/(g-1)+d

def celerity2velocity(w):

    import math

    return w/math.sqrt(pow(w,2)+1)

def shock_front_celerity(p,d0,p0,g):

    import math

    d = calc_hugoniot_density(p,d0,p0,g)
    e = calc_energy(d,p,g)
    e0 = calc_energy(d0,p0,g)
    numerator = (p-p0)*(e+p0)
    denominator = (e0+p0)*(e-e0-p+p0)
    return math.sqrt(numerator/denominator)

class Primitive:

    def __init__(self,d,p,w):

        self.d = d
        self.p = p
        self.w = w

class ShockProf:

    def __init__(self,p,d0,p0,g):

        self._ws = shock_front_celerity(p,d0,p0,g)
        self._vs = celerity2velocity(self._ws)
        self._upstream = Primitive(d0,p0,0)
        d = calc_hugoniot_density(p,d0,p0,g)
        w = calc_hugoniot_celerity(p,d0,p0,g)
        self._downstream = Primitive(d,p,w)

    def calcPrimitive(self, v):

        import copy

        if(v>self._vs):
            return copy.deepcopy(self._upstream)
        else:
            return copy.deepcopy(self._downstream)

def velocity_addition(v1, v2):

    return (v1+v2)/(1+v1*v2)

def calc_riemann_invariant(ba,g):

    import math

    return math.tanh(2*math.atanh(ba/math.sqrt(g-1))/math.sqrt(g-1))

def velocity2celerity(v):

    import math

    return v/math.sqrt(1-v**2)

class IsenProf:

    def __init__(self, p, d0, p0, g):

        self._g = g
        self._upstream = Primitive(d0,p0,0)
        self._csmax = calc_sound_speed(d0,p0,g)
        d = calc_isentrope_density(p,d0,p0,g)
        w = calc_isentrope_celerity(p,d0,p0,g)
        v = celerity2velocity(w)
        self._csmin = velocity_addition(v,calc_sound_speed(d,p,g))
        self._downstream = Primitive(d,p,w)

    def calcPrimitive(self, v):

        import scipy.optimize
        import copy
        import math

        g = self._g
        if v>self._csmax:
            return copy.deepcopy(self._upstream)
        elif v<self._csmin:
            return copy.deepcopy(self._downstream)
        else:

            trans_eqn = lambda ba: velocity_addition(v,-ba)-\
                math.tanh((2.0/math.sqrt(g-1))*
                     (math.atanh(ba/math.sqrt(g-1))-
                      math.atanh(self._csmax/math.sqrt(g-1))))
            ba = scipy.optimize.brentq(trans_eqn,
                                       calc_sound_speed(self._downstream.d,
                                                        self._downstream.p,
                                                        g),
                                       self._csmax)
            u = velocity_addition(v,-ba)
            w = velocity2celerity(u)
            d0 = self._upstream.d
            p0 = self._upstream.p
            d = d0*(g*p0*(ba**(-2)-1.0/(g-1))/d0)**(-1.0/(g-1))
            p = p0*(d/d0)**g
            return Primitive(d,p,w)

class HydroProf:

    def __init__(self,p,d0,p0,g):
        if p>p0:
            self._prof = ShockProf(p,d0,p0,g)
        else:
            self._prof = IsenProf(p,d0,p0,g)

    def calcPrimitive(self, v):

        return self._prof.calcPrimitive(v)

class RiemannProfile:

    def __init__(self, left, right, g):

        self._ps, self._ws = riemann_solve(left, right, g)
        self._vs = celerity2velocity(self._ws)
        self._left_prof = HydroProf(self._ps,
                                    left.d,
                                    left.p, g)
        self._right_prof = HydroProf(self._ps,
                                     right.d,
                                     right.p, g)
        self._left = left
        self._right = right

    def calcPrimitive(self, v):

        if v>self._vs:
            vr = celerity2velocity(self._right.w)
            res = self._right_prof.calcPrimitive(velocity_addition(v,-vr))
            res.w = celerity_addition(res.w,self._right.w)
        else:
            vl = celerity2velocity(self._left.w)
            res = self._left_prof.calcPrimitive(velocity_addition(vl,-v))
            res.w = celerity_addition(self._left.w,-res.w)
        return res

def test_riemann_profile():

    import numpy
    import pylab

    my_prof = RiemannProfile(Primitive(1,1,1),
                             Primitive(1,1,-1),5./3.)
    v_list = numpy.linspace(-0.99,0.99,1000)
    d_list = [my_prof.calcPrimitive(v).d for v in v_list]
    p_list = [my_prof.calcPrimitive(v).p for v in v_list]
    w_list = [my_prof.calcPrimitive(v).w for v in v_list]
    pylab.subplot(311)
    pylab.plot(v_list, d_list)
    pylab.subplot(312)
    pylab.plot(v_list, p_list)
    pylab.subplot(313)
    pylab.plot(v_list, w_list)
    pylab.show()

def test_riemann_profile_2():

    import numpy
    import pylab

    g = 5./3.
    my_prof = RiemannProfile(Primitive(1.,1000.,0.),
                             Primitive(1.,1e-2,0.),g)
    v_list = numpy.linspace(-0.99,0.99,10000)
    d_list = [my_prof.calcPrimitive(v).d for v in v_list]
    p_list = [my_prof.calcPrimitive(v).p for v in v_list]
    ba_list = [calc_sound_speed(d,p,g) for d,p in zip(d_list,p_list)]
    w_list = [my_prof.calcPrimitive(v).w for v in v_list]
    pylab.subplot(311)
    pylab.plot(v_list, d_list)
    pylab.subplot(312)
    pylab.plot(v_list, p_list)
    pylab.subplot(313)
    pylab.plot(v_list, w_list)
    #pylab.plot(v_list, ba_list)
    pylab.show()

if __name__ == '__main__':

    test_riemann_profile_2()
