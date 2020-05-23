''' Fitter_test_tool
Author: Tim Hucko
Used for testing fitting routines using fake data
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from iminuit import Minuit
from scipy.special import wofz
import random as rd


'''def lorentzian(x, p):
    # Here we use:(gamma/2)^2/((x-x0)^2+(gamma/2)^2)
    # Including and exponetial decay
    numerator = (p[2] / 2)**2
    decay = (np.exp(-x / p[4]))
    denominator1 = ((x - p[3]) ** 2 + (p[2] / 2) ** 2)
    L = p[0] + decay * p[1] * (numerator / denominator1)
    return L'''


def chi2(p):
    func = voigt(x_step, p)
    err = y_e # Uncertainty in the data sqrt(N)
    # calculate the chi2
    delta = (y2 - func) / err
    chi2_v = sum(pow(delta, 2))
    return chi2_v


def voigt(x, p):
    sigma = p[3]/(2*np.sqrt(2*np.log(2)))
    gamma = p[2]/2
    #print(sigma)
    z = ((x-p[4]) + 1j*gamma)/(sigma*np.sqrt(2))
    decay = (np.exp(-x / p[5]))
    num = np.real(wofz(z))
    dem = sigma*np.sqrt(2*np.pi)
    v = p[0] + decay*p[1]*(num/dem)
    return v



x_step = np.arange(0, 80, 1)
y =[]
y2 = []
a = 6000
b_width = 2.8
n_of_scans = 119.0
pps = 6.0
for i in range(len(x_step)):
    k = rd.gauss(a, np.sqrt(a))
    k = np.sign(k)*k
    #y.append(lorentzian(x[i], [k, 100, 3.5, 0, 50]))
    y2.append(voigt(x_step[i], [k, 6500000.90, 3.877, 1.4899, 40.34, 342]))
    y2[i] = ((y2[i] + rd.gauss(0, np.sqrt(y2[i]))))
    if y2[i] < 0:
        y2[i] = np.sign(y2[i])*y2[i]
y_e = np.sqrt(y2)
p1 = [np.min(y2),  np.max(y2), 3.5, 1,  x_step[np.argmax(np.array(y2))], 350]
#p1 = [10, 40000, 3.5, 0.6, 0, 350]

m = Minuit.from_array_func(chi2, p1, error=(0.1, 100, 0.01, 0.01, 1, 10),
                                   limit=((0, None), (0, None), (0.01, None), (0.01, None), (0, None), (0, None)),
                                   fix=(False, False, False, False, False, False), errordef=1, pedantic=False)

m.migrad()  # This is minimization strategy
p_fit = [m.values["x0"], m.values["x1"], m.values["x2"], m.values["x3"], m.values["x4"], m.values["x5"]]
p_err = [m.errors["x0"], m.errors["x1"], m.errors["x2"], m.errors["x3"], m.errors["x4"], m.errors["x5"]]
print(m.fval / (len(y2) - len(p_fit)))
print(p_fit)
x_fit = np.arange(0, 80, 0.01)
y_fit = voigt(x_fit, p_fit)

residuals = y2 - voigt(x_step, p_fit)
normres = (y2 - voigt(x_step, p_fit))/y_e

plt.style.use('ggplot')
gs = gridspec.GridSpec(2, 2)
fig = plt.figure(figsize=(16, 9))
ax1 = fig.add_subplot(gs[0, :])
ax1.errorbar(x_step, np.float64(y2)/(b_width*n_of_scans*pps), yerr=y_e/(b_width*n_of_scans*pps) , fmt='ro', ecolor='black', capsize=5, zorder=1)
ax1.plot(x_fit, y_fit/(b_width*n_of_scans*pps),  color='#3392FF', zorder=2)

ax2 = fig.add_subplot(gs[1, 0])
ax2.set_title('Normalized Residuals', fontsize=16)
ax2.set_xlabel('Frequency (MHz)', fontsize=14)
ax2.plot(x_step, normres, '-o', color='#F06C00')
ax3 = fig.add_subplot(gs[1, 1])
ax3.set_title('Residuals', fontsize=16)
ax3.set_xlabel('Frequency (MHz)', fontsize=14)
ax3.plot(x_step, residuals, '-o', color='#D433FF')
#fig.subplots_adjust(left=0.06, bottom=0.08, right=0.95, top=0.95, wspace=0.07, hspace=0.25)

plt.show()
