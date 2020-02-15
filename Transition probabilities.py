import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

def lorentzian(x, p, s0):
    # Here we use:(gamma/2)^2/((x-x0)^2+(gamma/2)^2)
    # Including and exponetial decay
    N = ((2/p[1])**2)*np.sqrt(1+s0)/(s0*np.pi)
    numerator = (s0/(s0+1))*(p[1] / 2)
    denominator = (1 + (2*(x - p[2])/(p[1]*np.sqrt(1+s0)))**2)
    L = N*(numerator / denominator)
    return L



laser = (266.65-15.0)
x = np.arange(-(laser+0.11), (laser+0.11), 0.01)
p = [1, 6.065, 0]
s0 = [0.1, 1, 10, 100]
section = np.arange(laser, laser+0.1, 0.001)


area0 = quad(lorentzian, x[0], x[len(x)-1], args=(p, s0[0]))
area1 = quad(lorentzian, x[0], x[len(x)-1], args=(p, s0[1]))
area2 = quad(lorentzian, x[0], x[len(x)-1], args=(p, s0[2]))
area3 = quad(lorentzian, x[0], x[len(x)-1], args=(p, s0[3]))


print(area0[0], area1[0], area2[0], area3[0])
'''a0 = np.trapz(norm_L0, dx=1)
a1 = np.trapz(norm_L1, dx=1)
a2 = np.trapz(norm_L2, dx=1)
a3 = np.trapz(norm_L3, dx=1)

print(a0*100, a1*100, a2*100, a3*100)'''

plt.style.use('ggplot')
plt.plot(x, lorentzian(x, p, s0[0]), label=r'$s0=%.2f$' %s0[0])
plt.fill_between(section, lorentzian(section, p, s0[0]), color='grey')
plt.legend()
plt.show()