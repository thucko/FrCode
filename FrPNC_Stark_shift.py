"""
Stark Shift Plotter
Author: Tim Hucko
Version: 0
"""


import numpy as np
from plotnine import *
from iminuit import Minuit
from tkinter.filedialog import askopenfilename
from tkinter import *
import pandas as pd


def quad (x, p):
    y = p[0]*x+p[1]
    return y


def file_dialog():
    Tk().withdraw()
    file = askopenfilename()
    return file


def chi2(p):
    func = quad(x_data, p)
    err = y_err
    # calculate the chi2
    delta = (y_data - func) / err
    chi2_v = sum(pow(delta, 2))
    return chi2_v

# Define some useful sci-notation
sci = {'m': 10**-3, 'c': 10**-2, 'k': 10**3, 'M': 10**6, 'G': 10**9}


f = file_dialog()
volt = np.genfromtxt(f, dtype='str', usecols=0)
pos = np.genfromtxt(f, dtype='float64', usecols=1)
pos_err = np.genfromtxt(f, dtype='float64', usecols=2)
data = pd.DataFrame()
E = np.array([])
for val in volt:
    for i, c in enumerate(val):
        if c.isalpha():
            num = np.float64(val[0:i])
            index = i
            V = num*sci.get(val[index])
            E = np.append(E, num / 2.858)
            break

for j in range(0, int(len(volt)/2)):
    data = data.append({
        'EField': (E[2*j]**2 - E[2*j+1]**2),
        'Frequency': -(pos[2*j] - pos[2*j+1]),
        'Uncertainty': np.sqrt(pos_err[2*j]**2+pos_err[2*j+1]**2)
    }, ignore_index=True)

y_err = data['Uncertainty']*np.sqrt(38.660597)
y_data = data['Frequency']
x_data = data['EField']
p1 = [-0.5,0]

m = Minuit.from_array_func(chi2, p1, error=(0.01, 0.01), limit=(None, None), fix=(False, False)
                           , errordef=1, pedantic=True)
m.migrad()


p_fit = [m.values['x0'], m.values['x1']]
p_err = [m.errors['x0'], m.errors['x1']]
Red_chi2 = chi2(p_fit) / (len(y_data) - len(p_fit))
print(Red_chi2)
print(p_fit)

x_fit = np.arange(np.min(x_data), np.max(x_data), 0.01)
y_fit = quad(x_fit, p_fit)


fit = pd.DataFrame({
    'y': y_fit,
    'x': x_fit
})

g1 = (ggplot(data, aes(x='EField', y='Frequency'))
            + ggtitle('Frequency Shift vs Electric Field')
            + geom_point(color='red')
            + geom_errorbar(aes(x='EField', ymin='Frequency-Uncertainty'
                                , ymax='Frequency+Uncertainty'))
            + geom_line(fit, aes(x='x', y='y'), color='blue')
            +ylab(r'$\Delta\nu_{Stark}$(MHZ)')
            +xlab('E(kV cm$^{-1}$)')
      )
print(g1)


exit()