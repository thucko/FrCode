import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import matplotlib.gridspec as gridspec

def lorentzian(x, p, s0):

    N = ((2/p[1])**2)*np.sqrt(1+s0)/(s0*np.pi)
    numerator = (s0/(s0+1))*(p[1] / 2)
    denominator = (1 + (2*(x - p[2])/(p[1]*np.sqrt(1+s0)))**2)
    L = N*(numerator / denominator)
    return L


h = 6.62607004E-34 # Plancks Constant
c = 3e8 # speed of light m/s
t = 26.24E-9 # lifetime
l = 780e-9

Is = (np.pi*h*c)/(3*t*l**3) #saturation intensity

Is = Is/1e-4

I = 156e-3/(np.pi*2.54**2) # light intensity W/cm^2
s1 = I/Is

laser = (266.65-15.0)
x = np.arange(-(laser+0.5), (laser+1), 0.1)
p = [1, 6.065, 0]
s0 = np.sort([0.1, 1, 10, 100, s1])
section = np.arange(laser, laser+0.105, 0.01)
plot_sec = np.arange(laser-.01, laser+0.2, 0.01)



area0 = quad(lorentzian, section[0], section[len(section)-1], args=(p, s0[0]))
area1 = quad(lorentzian, section[0], section[len(section)-1], args=(p, s0[1]))
area2 = quad(lorentzian, section[0], section[len(section)-1], args=(p, s0[2]))
area3 = quad(lorentzian, section[0], section[len(section)-1], args=(p, s0[3]))
area4 = quad(lorentzian, section[0], section[len(section)-1], args=(p, s0[4]))


print(area0[0]*100, area1[0]*100, area2[0]*100, area3[0]*100)
val = [area0[0]*100, area1[0]*100, area2[0]*100, area3[0]*100, area4[0]*100]
var = ord('%')
plt.style.use('ggplot')
gs = gridspec.GridSpec(2, 2)
fig = plt.figure(figsize=(16, 9))
ax1 = fig.add_subplot(gs[0, :])
ax1.set_title(r'Normalized Lorentzian with saturation parameter $s0$')
ax1.set_ylabel('Probability')
ax1.set_xlabel('Frequency Detuning (MHz)')
ax1.plot(x, lorentzian(x, p, s0[0]), label=r'$s0=$%.2E' % s0[0])
ax1.plot(x, lorentzian(x, p, s0[1]), label=r'$s0=%.2f$' % s0[1])
ax1.plot(x, lorentzian(x, p, s0[2]), label=r'$s0=%.2f$' % s0[2])
ax1.plot(x, lorentzian(x, p, s0[3]), label=r'$s0=%.2f$' % s0[3])
ax1.plot(x, lorentzian(x, p, s0[4]), label=r'$s0=%.2f$' % s0[4])
ax1.legend()
ax2 = fig.add_subplot(gs[1, :])
ax2.set_title('Probability of exciting F=2 transition @ %.2f MHz detuned' % laser)
ax2.plot(plot_sec, lorentzian(plot_sec, p, s0[0]))
ax2.plot(plot_sec, lorentzian(plot_sec, p, s0[1]))
ax2.plot(plot_sec, lorentzian(plot_sec, p, s0[2]))
ax2.plot(plot_sec, lorentzian(plot_sec, p, s0[3]))
ax2.plot(plot_sec, lorentzian(plot_sec, p, s0[4]))
ax2.fill_between(section, lorentzian(section, p, s0[0]), alpha=0.5, label='A=%.2E %c' % (val[0], var))
ax2.fill_between(section, lorentzian(section, p, s0[1]), alpha=0.35, label='A=%.2E %c' % (val[1], var))
ax2.fill_between(section, lorentzian(section, p, s0[2]), alpha=0.25, label='A=%.2E %c' % (val[2], var))
ax2.fill_between(section, lorentzian(section, p, s0[3]), alpha=0.15, label='A=%.2E %c' % (val[3], var))
ax2.fill_between(section, lorentzian(section, p, s0[4]), alpha=0.15, label='A=%.2E %c' % (val[4], var))
ax2.set_ylabel('Probability')
ax2.set_xlabel('Frequency Detuning (MHz)')
ax2.set_xlim(plot_sec[0], plot_sec[len(plot_sec)-1])

ax2.legend()
plt.show()