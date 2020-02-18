import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import matplotlib.gridspec as gridspec

def lorentzian(x, p, s0):
    # Here we use:(gamma/2)^2/((x-x0)^2+(gamma/2)^2)
    # Including and exponetial decay
    N = ((2/p[1])**2)*np.sqrt(1+s0)/(s0*np.pi)
    numerator = (s0/(s0+1))*(p[1] / 2)
    denominator = (1 + (2*(x - p[2])/(p[1]*np.sqrt(1+s0)))**2)
    L = N*(numerator / denominator)
    return L



laser = (266.65-15.0)
x = np.arange(-(laser+0.5), (laser+1), 0.01)
p = [1, 6.065, 0]
s0 = [0.1, 1, 10, 100]
section = np.arange(laser, laser+0.105, 0.01)


area0 = quad(lorentzian, section[0], section[len(section)-1], args=(p, s0[0]))
area1 = quad(lorentzian, section[0], section[len(section)-1], args=(p, s0[1]))
area2 = quad(lorentzian, section[0], section[len(section)-1], args=(p, s0[2]))
area3 = quad(lorentzian, section[0], section[len(section)-1], args=(p, s0[3]))


print(area0[0]*100, area1[0]*100, area2[0]*100, area3[0]*100)
val = [area0[0]*100, area1[0]*100, area2[0]*100, area3[0]*100]
var = ord('%')
plt.style.use('ggplot')
gs = gridspec.GridSpec(2, 2)
fig = plt.figure(figsize=(16, 9))
ax1 = fig.add_subplot(gs[0, :])
ax1.set_title(r'Normalized Lorentzian with saturation parameter $s0$')
ax1.set_ylabel('Probability')
ax1.plot(x, lorentzian(x, p, s0[0]), label=r'$s0=%.2f$' % s0[0])
ax1.plot(x, lorentzian(x, p, s0[1]), label=r'$s0=%.2f$' % s0[1])
ax1.plot(x, lorentzian(x, p, s0[2]), label=r'$s0=%.2f$' % s0[2])
ax1.plot(x, lorentzian(x, p, s0[3]), label=r'$s0=%.2f$' % s0[3])
ax1.legend()
ax2 = fig.add_subplot(gs[1, :])
ax2.plot(x, lorentzian(x, p, s0[0]), label=r'$s0=%.2f$' % s0[0])
ax2.plot(x, lorentzian(x, p, s0[1]), label=r'$s0=%.2f$' % s0[1])
ax2.plot(x, lorentzian(x, p, s0[2]), label=r'$s0=%.2f$' % s0[2])
ax2.plot(x, lorentzian(x, p, s0[3]), label=r'$s0=%.2f$' % s0[3])
ax2.fill_between(section, lorentzian(section, p, s0[0]), alpha=0.5, label='A=%.6f %c' % (val[0], var))
ax2.fill_between(section, lorentzian(section, p, s0[1]), alpha=0.35, label='A=%.6f %c' % (val[1], var))
ax2.fill_between(section, lorentzian(section, p, s0[2]), alpha=0.25, label='A=%.6f %c' % (val[2], var))
ax2.fill_between(section, lorentzian(section, p, s0[3]), alpha=0.15, label='A=%.6f %c' % (val[3], var))
ax2.set_xlim(section[0]-0.01, section[len(section)-1]+0.01)
ax2.set_ylim(0, 0.00017)
ax2.legend()
plt.show()