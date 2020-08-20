"""
To read and plot scope waveforms from CSV files

Author: Tim Hucko
Version 0

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit
import pandas as pd
from tkinter.filedialog import askopenfilenames
from tkinter.filedialog import askdirectory
from tkinter import *
import os
from itertools import product

is_fft = False
is_norm = True
is_log = False
do_fit = True

def file_dialog():
    Tk().withdraw()
    file = askopenfilenames(initialdir='~/Documents')
    return file


def dir_dialog():
    Tk().withdraw()
    pwd = askdirectory()
    return pwd

def powermW(x):
    P = np.float64(10**(x/10))
    return P

def fun(x, a, b, c):
    f = a * np.exp(-x/b) + c
    # f = c*(a-np.exp(-x/b))
    return f


class ScopeWaveform:
    def __init__(self):
        self.channels = {}

    def get_data(self, files):

        for x in files:
            file_name = os.path.basename(x)
            d = pd.read_csv(x, dtype='a')
            col_name = d.columns
            if d[col_name[0]][0] != 0:
               d[col_name[0]] = np.float64(d[col_name[0]])-min(np.float64(d[col_name[0]]))
            if is_fft == True:
                d[col_name[0]] = (np.float64(d[col_name[0]]))*(np.float64(d[col_name[-1]][0]))
            else:
                d[col_name[0]] = (np.float64(d[col_name[0]]))/(np.float64(d[col_name[-1]][0]))
            for j in col_name[1:-1]:
                d[j] =(np.float64(d[j]))
                if is_norm == True:
                    mx = np.max(d[j])
                    d[j] = (np.float64(d[j]/mx))
                if is_log == True:
                    d[j] = np.log(d[j])





            d.drop('Increment', axis=1, inplace=True)
            self.channels.update({file_name: d})


if __name__ == '__main__':
    files = file_dialog()
    data = ScopeWaveform()
    data.get_data(files)
    plt.style.use('ggplot')
    for f in files:
        fn = os.path.basename(f)
        df = data.channels[fn]
        col_name = df.columns
        pair = list(product(col_name[0], col_name[1:]))
        x = np.array(np.float64(df[pair[0][0]]))
        if do_fit == True:
            y = np.array(np.float64(df[pair[0][1]]))
            popt, pcov = curve_fit(fun, x, y, maxfev=1000)
            print(pcov)
        for i in range(0, len(pair)):
            plt.plot(x, df[pair[i][1]], label=pair[i][1])
    if do_fit == True:
        plt.plot(x, fun(x, *popt), 'k', label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
        # plt.text(1, 0.8, r'Fit Function: $f(x) = c(a-e^{-bx})$', size= 16)
        plt.text(1, 0.8, r'Fit Function: $f(x) = ae^{-x/b}+c$', size=16)


    plt.title(r'AOM decay')
    if is_fft == True:
        plt.xlabel('Frequency (kHz)')
        plt.ylabel('Power (dBm)')
    else:
        if is_norm == True:
            plt.ylabel('Normalized Voltage')
        else:
            plt.ylabel('Voltage')
        plt.xlabel(r'Time ($\mu$s)')
    plt.legend()
    plt.show()


