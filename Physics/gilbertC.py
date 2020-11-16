import numpy as np
from itertools import product
from sympy.physics.quantum.cg import Wigner6j, CG, Wigner3j
from sympy import evalf


class GilbertCoefficients:


    def coefficient(self, q, F, F_prime, m, m_prime, I):
        cg = CG(F, m, 1, q, F_prime, m_prime)
        #cg = CG(F, m, F_prime, m_prime, F+F_prime, )
        cg = cg.doit()
        w6j = Wigner6j(0.5, I, F, F_prime, 1, 0.5)
        w6j = w6j.doit()
        expectiation_val = 2*pow(-1, F+I-0.5)*cg*w6j*pow(1.5*(2*F+1), 0.5)
        c_dum = pow(1/np.sqrt(2), abs(q))*expectiation_val
        self.Cval = c_dum.evalf()


if __name__ == '__main__':
    #test with cesium numbers and compare with Gilbert or Bennett's theses.
    q = -1
    Ff = 4
    mf = 3
    Fi = 4
    mi = mf-q
    gc= GilbertCoefficients()
    gc.coefficient(q=q, F=Fi, F_prime=Ff, m=mi, m_prime=mf, I=7/2)
    gc.Cval
    print(gc.Cval)