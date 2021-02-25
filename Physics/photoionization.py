import numpy as np
import matplotlib.pyplot as plt


h = 6.626e-34
c = 3e8
wavelength = 506e-9

power_buildup = 3882
input_power = 50e-3 # laser power in watts
area = np.pi*0.0207**2 # beam area in cm^2
extinction_ratio = [1, 1/10000, 1/400000, 1/1000000]

photo_cross_section = 2.08e-17 # taken for photoionization cross section paper
n_e = 0.1215
photon_energy = h*(c/wavelength)

rate = []

for ext in extinction_ratio:
    power_eff = input_power * power_buildup * ext
    intensity = power_eff / area
    photon_flux = intensity / photon_energy
    rate.append(0.2*photo_cross_section*photon_flux*n_e)

time = np.arange(0, 5, 0.00001)
N0 = 1e5

N = []
for R in rate:
    print(R,":",1/R)
    N.append(N0*np.exp(-R*time))

plt.style.use('../matplotlib_style/stylelib/cern_root.mplstyle')
plt.plot(time, N[0], color='r', label='No extinction')
plt.plot(time, N[1], color='g', label='10000:1')
plt.plot(time, N[2], color='b', label='400000:1')
plt.plot(time, N[3], color='m', label='1000000:1')
plt.xlabel('Time (s)', size=16)
plt.ylabel('Atom Number', size=16)
plt.title('Trap Lifetime for Different Extinction Ratios', size=24)
plt.legend()
plt.show()
