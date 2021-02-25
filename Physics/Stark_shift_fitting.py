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


def quad_sat(x, p):
    y = p[0]*(x**2/(1+(x/p[1])**2)) + p[2]
    return y

def quad(x, p):
    y = p[0]*x**2 + p[1]
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
    func = quad_sat(x_data[0:len(x_data)], p)
    err = y_err[0:len(y_err)]
    # calculate the chi2
    delta = (y_data[0:len(y_data)] - func)/err
    chi = sum(pow(delta, 2))
    return chi

def chi22(p):
    func = quad(x_data[0:len(x_data)], p)
    err = y_err[0:len(y_err)]
    # calculate the chi2
    delta = (y_data[0:len(y_data)] - func)/err
    chi = sum(pow(delta, 2))
    return chi



# Define some useful sci-notation
sci = {'m': 10**-3, 'c': 10**-2, 'k': 10**3, 'M': 10**6, 'G': 10**9}

f = file_dialog()

vals = pd.read_csv(f, sep='\t')

#data0 = pd.DataFrame(columns=['EField', 'Frequency', 'Uncertainty'])
data0 = pd.DataFrame(columns=['EField', 'Peak Rate', 'Uncertainty'])
global sel2
if 'forward' in f:
    sel2 = 'Forward'
elif 'backward' in f:
    sel2 = 'Backward'


i = 0
for var in vals['Voltage']:
    conversion = sci.get(var[(var.find('V')-1)])
    volt = (np.float(var[0:var.find(' ')])*conversion)/1E3
    E = volt / 2.858
    data0 = data0.append({
        'EField': E,
        #'Frequency': vals['Peak Position'].iloc[i],
        'Peak Rate': vals['Peak Rate'].iloc[i],
        #'Uncertainty': vals['Peak Position error'].iloc[i]
        'Uncertainty': vals['Error'].iloc[i]
    }, ignore_index=True)
    i = i+1

'''data = pd.DataFrame(columns=['EField', 'Frequency', 'Uncertainty'])
for j in range(0, len(vals)):
    dE = data0['EField'].iloc[j] - data0['EField'].iloc[-1]
    dF = (data0['Frequency'].iloc[-1]-data0['Frequency'].iloc[j])
    dU = np.sqrt(data0['Uncertainty'].iloc[-1]**2 + data0['Uncertainty'].iloc[j]**2)
    data = data.append({
        'EField': dE,
        'Frequency': dF,
        'Uncertainty': dU
    }, ignore_index=True)'''



'''for j in range(0, int(len(volt)/2)):
    data = data.append({
        'EField': (E[2*j]**2 - E[2*j+1]**2),
        'Frequency': -(pos[2*j] - pos[2*j+1]),
        'Uncertainty': np.sqrt(pos_err[2*j]**2+pos_err[2*j+1]**2)
    }, ignore_index=True)'''

y_err = data0['Uncertainty']
y_data = data0['Peak Rate']
x_data = data0['EField']
'''y_err = data['Uncertainty']
y_data = data['Frequency']
x_data = data['EField']'''
p1 = [1.2, 5, 0.0]

m1 = Minuit.from_array_func(chi2, p1, error=(0.1, 0.1, 0.1), limit=(None, None, None), fix=(False, False, False), errordef=1)
m1.migrad()

p2 = [1.2, 0.0]
m2 = Minuit.from_array_func(chi22, p2, error=(0.1, 0.1), limit=(None, None), fix=(False, False), errordef=1)
m2.migrad()


p1_fit = [m1.values['x0'], m1.values['x1'], m1.values['x2']]
p1_err = [m1.errors['x0'], m1.errors['x1'], m1.errors['x2']]
p2_fit = [m2.values['x0'], m2.values['x1']]
p2_err = [m2.errors['x0'], m2.errors['x1']]
Red1 = m1.fval/(len(y_data) - len(p1_fit))
Red2 = m2.fval/(len(y_data) - len(p2_fit))
print(Red1)
print(Red2)
print('delta x0 = %.10f' % p1_fit[0])
print('delta x1 = %.10f' % p1_fit[1])
print('delta x2 = %.10f' % p1_fit[2])
#print('delta x3 = %.10f' % p_fit[3])

x_fit = np.arange(np.min(x_data), np.max(x_data), 0.01)
y1_fit = quad_sat(x_fit, p1_fit)
y2_fit = quad(x_fit, p2_fit)
residual = y_data - quad_sat(x_data, p1_fit)

repack_data = pd.DataFrame({
    'EField': x_data,
    'Frequency': y_data,
    'Uncertainty': y_err,
    'Residuals': residual,
    'Normalized Residuals': residual/y_err
})



fit = pd.DataFrame({
    'y1': y1_fit,
    'y2': y2_fit,
    'x': x_fit
})


plt.style.use('../matplotlib_style/stylelib/cern_root.mplstyle')
plt.errorbar(repack_data['EField'], repack_data['Frequency'], yerr=repack_data['Uncertainty'], fmt='ro', ecolor='black',
                 capsize=5, zorder=1)
plt.plot(fit['x'], fit['y1'], color='b', label=r'$\chi^2_{red}$=%.2f' %Red1)
plt.plot(fit['x'], fit['y2'], color='g', linestyle='--', label=r'$\chi^2_{red}$=%.2f' %Red2)
plt.title('Peak Rates vs Electric Field\n %s Scan' % sel2, fontsize=18)
plt.xlabel(r'$E$(kV/cm)', fontsize=16)
plt.ylabel(r'Peak Rates (kHz)', fontsize=16)
#plt.ylabel(r'$\Delta\nu_{Stark}$ (MHz)', fontsize=16)
'''gs = gridspec.GridSpec(2, 2)
fig = plt.figure(figsize=(16, 9))
ax1 = fig.add_subplot(gs[0, :])
#ax1.text(x_fit[0] + 1, np.sign(p_fit[0])*2.9, r'$\chi^2_{reduced}=%.3f$' % Red_chi2, fontsize=14)
ax1.text(x_fit[0] + 1, np.sign(p_fit[0])*2.6, r'$f(x) = %.3f x+(%.3f)$' % (p_fit[0], p_fit[1]), fontsize=14)
#ax1.text(x_fit[0] + 1, np.sign(p_fit[0])*2.3, r'$\delta m=%.3f$' % p_err[0], fontsize=14)
#ax1.text(x_fit[0] + 1, np.sign(p_fit[0])*2.0, r'$\delta b=%.3f$' % p_err[1], fontsize=14)
ax1.set_title('Frequency Shift vs Electric Field\n %s Scan' % sel2, fontsize=18)
ax1.set_xlabel(r'$E$(V/cm)', fontsize=16)
ax1.set_ylabel(r'$\Delta\nu_{Stark}$ (MHz)', fontsize=16)
ax1.errorbar(repack_data['EField'], repack_data['Frequency'], yerr=repack_data['Uncertainty'], fmt='ro', ecolor='black',
                 capsize=5, zorder=1)
ax1.plot(fit['x'], fit['y'], color='#3392FF', zorder=2)
ax2 = fig.add_subplot(gs[1, 0])
ax2.set_title('Normalized Residuals', fontsize=16)
ax2.set_xlabel(r'$E$(V/cm)', fontsize=14)
ax2.plot(repack_data['EField'], repack_data['Normalized Residuals'], '-o', color='#F06C00')
ax3 = fig.add_subplot(gs[1, 1])
ax3.set_title('Residuals', fontsize=16)
ax3.set_xlabel(r'$E$(V/cm)', fontsize=14)
ax3.plot(repack_data['EField'], repack_data['Residuals'], '-o', color='#D433FF')
fig.subplots_adjust(left=0.06, bottom=0.08, right=0.95, top=0.93, wspace=0.12, hspace=0.25)'''
plt.legend()
plt.show()

'''fig_save = file_save_plt()
if fig_save is None:
    None
else:
    with PdfPages(fig_save) as pdf:
        pdf.savefig()'''

exit()
