import numpy as np
from sympy import Symbol, im, I, Abs, sqrt
from sympy.physics.optics import *
from sympy.solvers import solve


lamb = 496e-7  # Wavelength (cm)
wb = 0.1  # Beam size
t = 0.635  # thickness of the mirror
f1 = -10.0  # focus of first lens
f2 = 10.0  # focus of second lens
n1 = 1.0  # index of refraction of air
n2 = 1.50  # index of refraction of glass
k = 7.35  # half the cavity length
rc = 100.00  # radius of curvature of the mirror
x1 = Symbol('x1', real=True)
x2 = 40.0
'''Calculate beam properties within the cavity'''
w0 = np.sqrt(lamb/(2*np.pi)*np.sqrt(2*k*(2*rc-2*k)))
div0 = lamb/(np.pi*w0)
R = k*(1+((np.pi*w0**2)/(lamb*k))**2)




M = FreeSpace(k)*CurvedRefraction(rc, n2, n1)*FreeSpace(t)*FlatRefraction(n1, n2)*FreeSpace(x2)*ThinLens(f2)\
    *FreeSpace(x1)*ThinLens(f1)
'''Calculate complex beam parameters'''
Qi = -(lamb*I/(np.pi*wb**2))  # initial beam parameters
q1 = [1/Qi, 1]
m11 = M[0, 0]
m12 = M[0, 1]
m21 = M[1, 0]
m22 = M[1, 1]
q2 = (m21*q1[0] + m22)/(m11*q1[0] + m12)
w = Abs(im(q2))
wf = sqrt(lamb/(np.pi*w))
print('Solving...')
print(solve(wf-w0, x1))
print()
