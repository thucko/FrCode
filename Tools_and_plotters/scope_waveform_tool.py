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
is_norm = False
is_log = False
do_fit = False
do_avg = True
laser_norm = True
t_avg = 2


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


def fun(x, a, b, c):
    f = a * np.exp(-x / b) + c
    # f = c*(a-np.exp(-x/b))
    return f


class ScopeWaveform:
    def __init__(self):
        self.channels = {}
        self.x_vals = []
        self.y_avg = []
        self.x_avg = []

    def get_data(self, files):

        for x in files:
            file_name = os.path.basename(x)
            d = pd.read_csv(x, dtype='a')
            col_name = d.columns
            if d[col_name[0]][0] != 0:
                d[col_name[0]] = np.float64(d[col_name[0]]) - min(np.float64(d[col_name[0]]))
            if is_fft is True:
                d[col_name[0]] = (np.float64(d[col_name[0]])) * (np.float64(d[col_name[-1]][0]))
            else:
                d[col_name[0]] = (np.float64(d[col_name[0]])) * (np.float64(d[col_name[-1]][0]))
            if do_avg is True:
                x_vals = d[col_name[0]]
                #ratio = (len(x_vals) - 1) / x_vals.iloc[-1]
                n = int(t_avg/(np.float64(d[col_name[-1]][0])))  # get for reshape, number in front is average over
                # time
                len_range = d.shape[0]/(n*10)
                self.x_avg = np.arange(2, (x_vals.iloc[-1]+(np.float64(d[col_name[-1]][0]))), 2)/len_range
            i = 0
            for j in col_name[1:-1]:
                d[j] = (np.float64(d[j]))
                if is_norm is True:
                    mx = np.max(d[j])
                    d[j] = (np.float64(d[j] / mx))
                if is_log is True:
                    d[j] = np.log(d[j])
                if laser_norm is True and j != 'Laser':
                    laser_val = np.float64(d['Laser'])
                    max_laser = np.max(laser_val)
                    norm_laser = laser_val/max_laser
                    d[j] = (np.float64(d[j] / norm_laser))

                if do_avg is True and j != 'Laser':
                    yvals = d[j].to_numpy()
                    self.y_avg.append(np.mean(yvals.reshape(-1, n), axis=1))
                    print(j, np.mean(self.y_avg[i]), np.std(self.y_avg[i]))
                    min_val = min(yvals)
                    max_val = max(yvals)
                    precent = 100*(max_val-min_val)/max_val
                    print(min_val, max_val, 'Percentage dip: %.3f' %precent)
                    i = i+1

                    print('done')

            d.drop('Increment', axis=1, inplace=True)
            self.channels.update({file_name: d})


if __name__ == '__main__':
    files = file_dialog()
    data = ScopeWaveform()
    data.get_data(files)
    plt.style.use('../matplotlib_style/stylelib/cern_root.mplstyle')
    colours = ['b', 'r', 'g']
    for f in files:
        fn = os.path.basename(f)
        df = data.channels[fn]
        col_name = df.columns
        pair = list(product(col_name[0], col_name[1:]))
        x = np.array(np.float64(df[pair[0][0]]))
        if do_fit is True:
            y = np.array(np.float64(df[pair[0][1]]))
            popt, pcov = curve_fit(fun, x, y, maxfev=1000)
            perr = np.sqrt(np.diag(pcov))
            print(popt, perr)
        n = 0
        for i in range(0, len(pair)):

            if laser_norm is True:
                if pair[i][1] != 'Laser':
                    #plt.plot(x, df[pair[i][1]], label=(pair[i][1]+' (Normalized to laser power)'))
                    if do_avg is True and (pair[i][1] != 'Laser'):
                        plt.plot(data.x_avg, data.y_avg[n], '-*',label='Average of %i ms' % t_avg, linewidth=3.0, ms=10.0)
                        n = n+1

            else:
                plt.plot(x, df[pair[i][1]], 'o', label=pair[i][1])
                if do_avg is True and (pair[i][1] != 'Laser'):
                    plt.plot(data.x_avg, data.y_avg[n], '-*',label='Average of %i ms' % t_avg, linewidth=3.0, ms=10.0)
                    n = n + 1

    if do_fit is True:
        plt.plot(x, fun(x, *popt), 'k', label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
        # plt.text(1, 0.8, r'Fit Function: $f(x) = c(a-e^{-bx})$', size= 16)
        plt.text(1, 0.8, r'Fit Function: $f(x) = ae^{-x/b}+c$', size=16)


    plt.title(r'PBC Transmission (decay)')
    if is_fft is True:
        plt.xlabel('Frequency (kHz)')
        plt.ylabel('Power (dBm)')
    else:
        if is_norm is True:
            plt.ylabel('Normalized Voltage')
        else:
            plt.ylabel('Voltage')
        plt.xlabel(r'Time (minutes)')
    plt.ylim(0)
    plt.legend()
    plt.show()
