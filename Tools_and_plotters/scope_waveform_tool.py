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
from re import search

is_fft = False
is_norm = False
x1 = 0
x2 = 5
is_log = False
is_chop = False
do_fit = False
do_avg = False
laser_norm = False
ttl_time = 5
t_avg = 2
t_int = 2


def file_dialog():
    Tk().withdraw()
    file = askopenfilenames(initialdir='~/Documents')
    return file


def dir_dialog():
    Tk().withdraw()
    pwd = askdirectory()
    return pwd


def powermW(x):
    P = np.float64(10 ** (x / 10))
    return P

if is_log is True:
    def fun(x, a, b):
        f = a+b*x
        return f
else:
    def fun(x, a, b, c):
        f = a * np.exp(-x / b) + c
        #f = c * (a - np.exp(-x / b)
        return f

def lorentzian(x, a, b, c):
    numerator = (b / 2) ** 2
    denominator1 = ((x - (c)) ** 2 + (b / 2) ** 2)
    L = a * (numerator / denominator1)
    return L


'''def lorentzian(x, a, b, c):
    s0 =1e-10
    N = ((2/a)**2)*np.sqrt(1+s0)/(s0*np.pi)
    numerator = (s0/(s0+1))*(a / 2)
    denominator = (1 + (2*(x - b)/(a*np.sqrt(1+s0)))**2)
    L = c*N*(numerator / denominator)
    return L'''


class ScopeWaveform:
    def __init__(self):
        self.channels = {}
        self.x_vals = []
        self.y_avg = []
        self.x_avg = []
        self.y_chop = []
        self.x_chop = []

    def get_data(self, files):
        for x in files:
            file_name = os.path.basename(x)
            d = pd.read_csv(x, dtype='a')
            col_name = d.columns
            if d[col_name[0]][0] != 0:
                d[col_name[0]] = np.float64(d[col_name[0]]) - min(np.float64(d[col_name[0]]))
            if is_fft is True:
                d[col_name[0]] = (np.float64(d[col_name[0]])) / (np.float64(d[col_name[-1]][0]))
            else:
                d[col_name[0]] = (np.float64(d[col_name[0]])) / (np.float64(d[col_name[-1]][0]))
            if is_chop is True:
                chop_test = np.array_split(d[col_name[0]], 20)
                comb = chop_test[1::2]
                self.x_chop = np.concatenate(comb)
                if do_avg is True:
                    n = int(t_avg / (np.float64(d[col_name[-1]][0])))
                    mod_x = np.mod(len(self.x_chop), n)
                    if mod_x != 0:
                        x_vals = self.x_chop[:-mod_x]
                    else:
                        x_vals = self.x_chop
                    # ratio = (len(x_vals) - 1) / x_vals.iloc[-1]
                    # get for reshape, number in front is average over
                    # time
                    # len_range = d.shape[0] / (n * 10)
                    self.x_avg = x_vals[0::n]
            else:
                if do_avg is True:
                    n = 2000 #int(t_avg / (np.float64(d[col_name[-1]][0])))
                    #n = int(t_avg / 1e-6)
                    mod_x = np.mod(len(d[col_name[0]]), n)
                    if mod_x != 0:
                        x_vals = d[col_name[0]].to_numpy()[:-mod_x]
                    else:
                        x_vals = d[col_name[0]].to_numpy()
                    # ratio = (len(x_vals) - 1) / x_vals.iloc[-1]
                    # get for reshape, number in front is average over
                    # time
                    #len_range = d.shape[0] / (n * 10)
                    self.x_avg = x_vals[0::n]+t_avg
            i = 0
            for j in col_name[1:-1]:
                d[j] = (np.float64(d[j]))
                d_avg = 0.011
                if is_log is True:
                    d[j] = -np.log(d[j])
                if laser_norm is True and j != 'Laser':
                    laser_val = np.float64(d['Laser'])
                    max_laser = np.max(laser_val)
                    norm_laser = laser_val / max_laser
                    d[j] = (np.float64(d[j] / norm_laser))
                if is_norm is True:
                    pos1 = d[d[col_name[0]] == x1].index.values
                    pos2 = d[d[col_name[0]] == x2].index.values
                    mx = np.max(d[j][pos1[0]:pos2[0]])
                    d[j] = (np.float64(d[j]) / mx)

                if is_chop is True and j != 'TTL':
                    chop_test = np.array_split(d[j], 20)
                    comb = chop_test[1::2]
                    self.y_chop = np.concatenate(comb)


                if do_avg is True and j != 'Laser' and j != 'TTL':

                    if is_chop is True:
                        mod_val = np.mod(len(self.y_chop), n)
                        mod_val_2 = np.mod(len(self.y_chop) - mod_val, n)
                        yvals = self.y_chop[:-mod_val]
                        self.y_chop = self.y_chop[:-mod_val]
                        self.x_chop = self.x_chop[:-mod_val]
                        self.y_avg.append(np.mean(yvals.reshape(-1, n), axis=1))
                        r = x_vals[-1] / len(self.y_avg[i])
                        #self.x_avg = np.arange(x_vals[0], x_vals[-1], r)
                        print(j, np.mean(self.y_avg[i]), np.std(self.y_avg[i]))
                        min_val = min(yvals)
                        max_val = max(yvals)
                        precent = 100 * (max_val - min_val) / max_val
                        print(min_val, max_val, 'Percentage dip: %.3f' % precent)
                        i = i + 1

                    else:
                        mod_val = np.mod(len(d[j]), n)
                        if mod_val != 0:
                            yvals = d[j].to_numpy()[:-mod_val]
                        else:
                           yvals = d[j].to_numpy()
                        self.y_avg.append(np.mean(yvals.reshape(-1, n), axis=1))
                        print(j, np.mean(self.y_avg[i]), np.std(self.y_avg[i]))
                        min_val = min(yvals)
                        max_val = max(yvals)
                        precent = 100 * (max_val - min_val) / max_val
                        print(min_val, max_val, 'Percentage dip: %.3f' % precent)
                        i = i + 1

                    print('done')

            d.drop('Increment', axis=1, inplace=True)
            self.channels.update({file_name: d})


if __name__ == '__main__':
    files = file_dialog()
    data = ScopeWaveform()
    data.get_data(files)
    plt.style.use('../matplotlib_style/stylelib/cern_root.mplstyle')
    for f in files:
        fn = os.path.basename(f)
        df = data.channels[fn]
        col_name = df.columns
        pair = list(product(col_name[0], col_name[1:]))
        x = np.array(np.float64(df[pair[0][0]]))
        if do_fit is True:
            y = np.array(np.float64(df[pair[0][1]]))
            popt, pcov = curve_fit(fun, x, y)
            perr = np.sqrt(np.diag(pcov))
            print(popt, perr)
        n = 0
        for i in range(0, len(pair)):

            if laser_norm is True:
                if pair[i][1] != 'Laser':
                    plt.plot(x, df[pair[i][1]], '-', label=(pair[i][1]))
                    if do_avg is True and (pair[i][1] != 'Laser'):
                        plt.plot(data.x_avg, data.y_avg[n], '-*', label='Average of %i ms' % t_int)
                        n = n + 1

            else:
                plt.plot(x, df[pair[i][1]], '-.', label=pair[i][1])
                '''if is_chop is True and (pair[i][1] != 'TTL'):
                    #plt.plot(data.x_chop, data.y_chop, '-*', label=' Light ON ')'''
                if do_avg is True and (pair[i][1] != 'Laser') and (pair[i][1] != 'TTL'):
                    plt.plot(data.x_avg, data.y_avg[n], '-*', label='Average of %i ms' % t_int)
                    n = n + 1

    if do_fit is True:
        if is_log is True:
            plt.plot(x, fun(x, *popt), 'k', label='fit: a=%5.2f, b=%5.2f' % (popt[0], popt[1]))
            plt.text(2.5, 6, r'Fit Function: $f(x) = a+bx$', size=16)
        else:
            plt.plot(x, fun(x, *popt), 'k', label='fit: a=%5.2f, b=%5.2f, c=%5.2f' % tuple(popt))
            # plt.plot(x, lorentzian(x,1, 2.12715, 26.21020155))
            # plt.text(1, 0.8, r'Fit Function: $f(x) = c(a-e^{-bx})$', size= 16)
            plt.text(2, np.mean(y), r'Fit Function: $f(x) = ae^{-x/b} + c$', size=16)

    plt.title(r'PBC (locked) Transmission', size=30)
    if is_fft is True:
        plt.xlabel('Frequency (kHz)', size = 20)
        plt.ylabel('Power (dBm)')
    else:
        if is_norm is True:
            plt.ylabel('Normalized Voltage')
        elif is_log is True:
            plt.ylabel('log(V)')
        else:
            plt.ylabel('Voltage (V)', size = 20)
        plt.xlabel(r'Time (ms)', size = 20)
    #plt.ylim(-0.25,1.1)
    plt.ylim(0)
  #  plt.xlim(x1, x2)
    plt.legend()
    plt.show()
