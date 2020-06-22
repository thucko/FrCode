"""
Stark Shift Plotter
Author: Tim Hucko
Version: 1
- Using Matplotlib for Stark plots
- Plotting both Normalized residuals and regular residuals
- Using the stark data file produced by FrPNC_Spectra_fitter_v5.py, auto-detects if using forward or backward scans
"""


import numpy as np
from iminuit import Minuit
from tkinter.filedialog import askopenfilename
from tkinter.filedialog import asksaveasfile
from tkinter import *
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages


def quad(x, p):
    y = p[0]*x+p[1]
    return y


def file_save_plt():
    f = asksaveasfile(mode='wb', defaultextension=".pdf", title='Save Figure(s)',
                      initialfile='Stark_Shift_%s_scan_Plot' % sel2)
    if f is None:  # asksaveasfile return `None` if dialog closed with "cancel".
        return
    return f


def file_dialog():
    Tk().withdraw()
    file = askopenfilename()
    return file


def chi2(p):
    func = quad(x_data[0:len(x_data)-1], p)
    err = y_err[0:len(y_err)-1]
    # calculate the chi2
    delta = sum(((y_data[0:len(y_data)-1] - func)/err)**2)
    return delta



# Define some useful sci-notation
sci = {'m': 10**-3, 'c': 10**-2, 'k': 10**3, 'M': 10**6, 'G': 10**9}


f = file_dialog()
volt = np.genfromtxt(f, dtype='str', usecols=0)
pos = np.genfromtxt(f, dtype='float64', usecols=1)
pos_err = np.genfromtxt(f, dtype='float64', usecols=2)
data = pd.DataFrame()
global sel2
if '_Forward' in f:
    sel2 = 'Forward'
elif '_Backward' in f:
    sel2 = 'Backward'


E = np.array([])
for val in volt:
    for i, c in enumerate(val):
        if c.isalpha():
            num = np.float64(val[0:i])
            index = i
            V = num*sci.get(val[index])/1E3
            E = np.append(E, V / 2.858)
            break

for j in range(0, int(len(volt)/2)):
    data = data.append({
        'EField': (E[2*j]**2 - E[2*j+1]**2),
        'Frequency': -(pos[2*j] - pos[2*j+1]),
        'Uncertainty': np.sqrt(pos_err[2*j]**2+pos_err[2*j+1]**2)
    }, ignore_index=True)

y_err = data['Uncertainty']
y_data = data['Frequency']
x_data = data['EField']
p1 = [0.5, -0.03]

m = Minuit.from_array_func(chi2, p1, error=(0.05, 0.0), limit=(None, None), fix=(False, False)
                           , errordef=1)
m.migrad()


p_fit = [m.values['x0'], m.values['x1']]
p_err = [m.errors['x0'], m.errors['x1']]
Red_chi2 = chi2(p_fit) / (len(y_data) - len(p_fit))
print(m.fval/(len(y_data) - len(p_fit)))
print(Red_chi2)
print('delta m = %.4f' % p_err[0])
print('delta b = %.4f' % p_err[1])
x_fit = np.arange(np.min(x_data), np.max(x_data), 0.01)
y_fit = quad(x_fit, p_fit)
residual = y_data - quad(x_data, p_fit)

repack_data = pd.DataFrame({
    'EField': x_data,
    'Frequency': y_data,
    'Uncertainty': y_err,
    'Residuals': residual,
    'Normalized Residuals': residual/y_err
})



fit = pd.DataFrame({
    'y': y_fit,
    'x': x_fit
})


plt.style.use('ggplot')
gs = gridspec.GridSpec(2, 2)
fig = plt.figure(figsize=(16, 9))
ax1 = fig.add_subplot(gs[0, :])
#ax1.text(x_fit[0] + 1, np.sign(p_fit[0])*2.9, r'$\chi^2_{reduced}=%.3f$' % Red_chi2, fontsize=14)
ax1.text(x_fit[0] + 1, np.sign(p_fit[0])*2.6, r'$f(x) = %.3f x+(%.3f)$' % (p_fit[0], p_fit[1]), fontsize=14)
#ax1.text(x_fit[0] + 1, np.sign(p_fit[0])*2.3, r'$\delta m=%.3f$' % p_err[0], fontsize=14)
#ax1.text(x_fit[0] + 1, np.sign(p_fit[0])*2.0, r'$\delta b=%.3f$' % p_err[1], fontsize=14)
ax1.set_title('Frequency Shift vs Electric Field\n %s Scan' % sel2, fontsize=18)
ax1.set_xlabel(r'$E^2$(kV$^2$ cm$^{-2})$', fontsize=16)
ax1.set_ylabel(r'$\Delta\nu_{Stark}$ (MHz)', fontsize=16)
ax1.errorbar(repack_data['EField'], repack_data['Frequency'], yerr=repack_data['Uncertainty'], fmt='ro', ecolor='black',
                 capsize=5, zorder=1)
ax1.plot(fit['x'], fit['y'], color='#3392FF', zorder=2)
ax2 = fig.add_subplot(gs[1, 0])
ax2.set_title('Normalized Residuals', fontsize=16)
ax2.set_xlabel(r'$E^2$(kV$^2$ cm$^{-2})$', fontsize=14)
ax2.plot(repack_data['EField'], repack_data['Normalized Residuals'], '-o', color='#F06C00')
ax3 = fig.add_subplot(gs[1, 1])
ax3.set_title('Residuals', fontsize=16)
ax3.set_xlabel(r'$E^2$(kV$^2$ cm$^{-2})$', fontsize=14)
ax3.plot(repack_data['EField'], repack_data['Residuals'], '-o', color='#D433FF')
fig.subplots_adjust(left=0.06, bottom=0.08, right=0.95, top=0.93, wspace=0.12, hspace=0.25)
#plt.show()

fig_save = file_save_plt()
if fig_save is None:
    None
else:
    with PdfPages(fig_save) as pdf:
        pdf.savefig()

exit()