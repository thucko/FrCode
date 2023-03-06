#!/usr/bin/env python3

'''
- Using Matplotlib for subplots of data and residuals
- plots now contain normalized residuals and regular residuals
-
Author: Tim Hucko

'''

import numpy as np
import pandas as pd
from tkinter.filedialog import askdirectory
from tkinter.filedialog import askopenfilename
from tkinter.filedialog import asksaveasfile
from tkinter import messagebox
import glob
from iminuit import Minuit
from iminuit.util import propagate
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
from scipy.special import wofz, voigt_profile
import tkinter as tk
import matplotlib
import matplotlib.pyplot as plt
import warnings
import time
import pyautogui
import sys
import os
from scipy.integrate import quad
from iminuit.cost import LeastSquares

warnings.simplefilter(action='ignore', category=FutureWarning)

transition_type = r'Burst signal at 110 MHz, Power = 165 $\mu$W'

def lorentzian(x, p):
    # Here we use:(gamma/2)^2/((x-x0)^2+(gamma/2)^2) Normalized to hieght
    # (gamma/2pi) is normalized to area
    # Including an exponetial decay
    numerator = (p[2]/(2*np.pi))
    decay = np.exp(-x/p[4])
    denominator1 = ((x - (p[3])) ** 2 + (p[2] / 2) ** 2)
    background = p[0]
    L = background + p[1] * decay*(numerator / denominator1)
    return L


# Voigt profile
# see https://en.wikipedia.org/wiki/Voigt_profile for info on the Voigt profile
def voigt(x, p):
    sigma = p[3] / (2 * np.sqrt(2 * np.log(2)))
    gamma = (p[2] / 2)
    z = (x - p[4] + 1j * gamma) / (sigma * np.sqrt(2))
    decay = (np.exp(-x / p[5]))
    num = np.real(wofz(z))
    dem = sigma * np.sqrt(2 * np.pi)
    background = p[0]
    V = background + p[1]*voigt_profile((x-p[4]), sigma, gamma)*decay
    return V


# Parameters for minimization includes Voigt and Lorentzian parameters
def Parameters(data, x, f):
    # Lorentzain: p = [Offset, Height, FWHM, Center, lifetime]
    # Voigt: p = [Offset, Height, FWHM_L, FWHM_G, Center, Liftime]
    if f == 'lorentzian':
        p = np.array((np.min(data), np.max(data), 3.5, x[np.argmax(data)],400))
        return p
    elif f == 'voigt':

        p = np.array((np.min(data), np.max(data), 3.6, 1, x[np.argmax(data)], 450))
        # p = [np.min(data), np.max(data), 3.6, 1, device_put(x)[np.argmax(data)], 350]

        return p


# Chi squared function used in the minimization
def chi2(p):
    if user_states.fit_select == 'lorentzian':
        func = lorentzian(x_step, p)  # Our function we need
    elif user_states.fit_select == 'voigt':
        func = voigt(x_step, p)
    err = y_err  # Uncertainty in the data sqrt(N)
    # calculate the chi2
    delta = (y_data - func) / err
    chi2_v = sum(pow(delta, 2))
    return chi2_v


# This class handles the data for forward and backward scans
class DataManager:
    def __init__(self, name):
        self.name = name
        self.y_val = []
        self.x_val = []
        self.y_err = []

    def data_select(self, name, data, x, b_width, pps):

        if name == 'forward':
            a = 10
            b = 70
        elif name == 'backward':
            a = 80
            b = 160
        elif name == 'full':  # New option in case we sit at one frequency
            a = 5
            b = 155

        self.y_val.extend(data[a:b]/(b_width * pps))  # converts counts into rate (Hz)
        self.x_val.extend(x[a:b]*120/160)
        self.y_err.extend(np.sqrt(data[a:b]) / (b_width * pps)) # transforms the uncertainty into rate


# minimizaion class using iminuit
class MinimizationRoutine:

    def __init__(self, name):
        self.name = name
        self.fit = []
        self.err = []
        self.covar = np.float64
        self.mvalues = np.float64

    def minimzation(self, name, par):
        m = Minuit(chi2, par)
        if name == 'voigt':
            errors = (1, 1, 0.01, 0.001, .1, 10)
            limits = (None, None, (0.1, 9), (1, 6), None, (1, 600))
            fixed = (False, False, False, False, False, False)

        elif name == 'lorentzian':
            m.errors = (0.1,0.1,0.1,0.1,0.1)
            m.errordef = Minuit.LEAST_SQUARES
            m.fixed[0] = False
            m.fixed[1] = False
            m.fixed[2] = False 
            m.fixed[3] = False
            m.fixed[4] = False
            m.limits = ((0,None), (0,None), (0,None), (0,None), (0,None))
            

           

        m.migrad()

        if name == 'voigt':
            #               [Offset        , Height        , FWHM_L        , FWHM_G       , Center         , liftime       ]
            self.fit.extend([m.values["x0"], m.values["x1"], m.values["x2"], m.values["x3"], m.values["x4"], m.values["x5"]])
            self.err.extend([m.errors["x0"], m.errors["x1"], m.errors["x2"], m.errors["x3"], m.errors["x4"], m.errors["x5"]])

        elif name == 'lorentzian':
            #               [Offset        , Height        , FWHM_L        , Center        , lifetime      ]
            self.fit.extend([m.values["x0"], m.values["x1"], m.values["x2"], m.values["x3"], m.values["x4"]])
            self.err.extend([m.errors["x0"], m.errors["x1"], m.errors["x2"], m.errors["x3"], m.errors["x4"]])
            self.covar = m.covariance
            self.mvalues = m.values
            # self.fit.extend([m.values["x0"], m.values["x1"], m.values["x2"], m.values["x3"], 0])
            # self.err.extend([m.errors["x0"], m.errors["x1"], m.errors["x2"], m.errors["x3"], 0])
            # self.fit.extend([m.values["x0"], m.values["x1"], m.values["x2"], m.values["x3"], m.values["x4"], m.values["x5"]])
            #self.err.extend([m.errors["x0"], m.errors["x1"], m.errors["x2"], m.errors["x3"], m.errors["x4"], m.errors["x5"]])


# get time stamps from raw data file
class GetTimes:
    def __init__(self):
        self.list_files = pd.DataFrame(columns=['Time', 'File path'])

    def get_time(self, file_list):
        for file_path in file_list:
            f = open(file_path, 'r')
            cont = f.readlines()
            start_time = cont[6]
            time_ = start_time[13:28]
            time_ = time_.replace(' ', '')
            # pm_am = time.find('PM')
            if (time_.endswith('PM')):

                if time_[0:2] == '12':
                    time_ = time_.replace('PM', '')
                else:
                    time_ = time_.replace('PM', '')
                    time_ = str(0) + time_
                    hour = int(time_[0:2])
                    hour = hour + 12
                    time_ = str(hour) + time_[2::]
            elif (time_.endswith('AM')):
                time_ = time_.replace('AM', '')
                if time_[0:2] == '12':
                    time_ = time_.replace(time_[0:2], '00')
                else:
                    if (len(time_) != 12):
                        time_ = str(0) + time_
            f.close()
            self.list_files = self.list_files.append({
                'Time': time_,
                'File path': file_path
            }, ignore_index=True)
            self.list_files = self.list_files.sort_values(by='Time')

class Stark_data:
    def __init__(self, direction, sOi, func):
        self.direction = direction
        self.sum_or_indvidual = sOi
        self.user_function = func

    def get_stark_data(self, volt, peakpos, peakpos_err):
        stark_data = pd.DataFrame(columns=['Voltage', 'Peak Position', 'Peak Position error'])
        if self.sum_or_indvidual == 'individual':
            peakpos = np.array(peakpos)
            pp = np.mean(peakpos)
            pperr = np.std(peakpos)
        elif self.sum_or_indvidual == 'sum':
            peakpos = np.float(peakpos)
            peakpos_err = np.float(peakpos_err)
            pp = peakpos
            pperr = peakpos_err
        stark_data = stark_data.append({
            'Voltage': volt,
            'Peak Position': pp,
            'Peak Position error': pperr
        }, ignore_index=True)
        stark_data.to_csv('Stark_%s_%s_data' %(self.direction, self.sum_or_indvidual), sep='\t', index=False, mode='a',
                          header=False)
    def calc_peakrate(self, volt, fit_vals):
        rate_data = pd.DataFrame(columns=['Voltage', 'Peak Rate', 'Error'])
        pars = np.array(fit_vals.to_numpy())
        if self.user_function == 'voigt':
            par = pars[0][1:7]
            err = pars[0][8:14]
            rate = voigt(fit_vals['Peak Position'], par)
            drate = rate*np.sqrt((sum(pow(err/par,2))))

        rate_data = rate_data.append({
            'Voltage': volt,
            'Peak Rate': rate[0],
            'Error': drate[0]
        }, ignore_index = True)

        rate_data.to_csv('Stark_%s_%s_rates' % (self.direction, self.sum_or_indvidual), sep='\t', index=False, mode='a',
                          header=False)






class GUI:
    def __init__(self):
        self.fit_select = ''
        self.volt_select = ''
        self.scan_direc_select = ''
        self.summing_scan = ''
        self.data_dir = ''
        self.fit_path = ''
        self.err_path = ''
        self.figr_path = ''
        self.fit_name = ''  # file name for fit
        self.err_name = ''  # file name for err
        self.figr_name = '' # file name for figure
        self.save_figq = ''
        self.root = tk.Tk()
        self.var1 = tk.StringVar(self.root)  # variable for fit
        self.var2 = tk.StringVar(self.root)  # variable for scan direction
        self.var3 = tk.StringVar(self.root)  # variable for summing scans
        self.var4 = tk.StringVar(self.root)  # variable for voltage
        self.vardata = tk.StringVar(self.root)  # string variable for dir selection
        self.varfit = tk.StringVar(self.root)  # variable for saving fit data
        self.varerr = tk.StringVar(self.root)  # variable for saving fit error data
        self.varfig = tk.StringVar(self.root)  # variable for saving figure data
        self.varsavefig = tk.StringVar(self.root)


    def gui_loop(self):
        self.root.title('Spectra Fitter')
        screen_width, screen_height = pyautogui.size()
        popup_width, popup_height = 350, 500
        self.root.geometry('%dx%d+%d+%d' % (
            popup_width, popup_height, (screen_width - popup_width) / 2, (screen_height - popup_height) / 2))

        # Get dir for data
        #self.vardata.set=('~/Documents/')
        tk.Label(self.root, text='Select Data Directory:').place(x=0, y=5)
        dir_path = tk.Entry(self.root, textvariable=self.vardata, width=44)
        dir_path.place(x=0, y=25)
        img = tk.PhotoImage(file='../ui_icons/folder.png')
        img = img.subsample(100)
        dir_button = tk.Button(self.root, image=img, command=self.dir_dialog, width=25, height=20)
        dir_button.place(x=316, y=20)


        # State select buttons
        funcs = [('Lorentzian', 'lorentzian'), ('Voigt', 'voigt'), ('None', 'none')]
        scan_direction = [('Forward Scans', 'forward'), ('Backward Scans', 'backward'), ('Full', 'full')]
        scan_sum_state = [('Sum', 'sum'), ('Individual', 'individual')]
        y_state_sel = 75
        y1_init = 20
        y2_init = 20
        y3_init = 20
        tk.Label(self.root, text='Fit Function:').place(x=0, y=y_state_sel)
        for text, val in funcs:
            rb_fun = tk.Radiobutton(self.root, text=text, value=val, variable=self.var1, height=2, width=15, justify=tk.CENTER,
                                    indicatoron=0, selectcolor='green', command=self.get_vals)
            rb_fun.place(x=0, y=y_state_sel + y1_init)
            # rb_fun.pack(side=tk.LEFT)
            y1_init = y1_init + 40
        tk.Label(self.root, text='Scan Direction:').place(x=120, y=y_state_sel)
        for text, val in scan_direction:
            rb_scan = tk.Radiobutton(self.root, text=text, value=val, variable=self.var2, height=2, width=15, justify=tk.CENTER,
                                     indicatoron=0, selectcolor='green', command=self.get_vals)
            rb_scan.place(x=120, y=y_state_sel + y2_init)
            y2_init = y2_init + 40

        tk.Label(self.root, text='Sum/Indvidual:').place(x=240, y=y_state_sel)
        for text, val in scan_sum_state:
            rb_sum = tk.Radiobutton(self.root, text=text, value=val, variable=self.var3, height=2, width=15, justify=tk.CENTER,
                                    indicatoron=0, selectcolor='green', command=self.get_vals)
            rb_sum.place(x=240, y=y_state_sel + y3_init)
            y3_init = y3_init + 40

        tk.Label(self.root, text='Voltage Used:').place(x=0, y=y_state_sel + 130)
        voltage_ent = tk.Entry(self.root, textvariable=self.var4, width=30)
        voltage_ent.place(x=0, y=y_state_sel + 150)

        button_state_val = tk.Button(self.root, text="Okay", command=self.get_vals, height=1, width=15)
        button_state_val.place(x=218, y=y_state_sel + 145)


        # Save the fit data

        self.data_dir = self.vardata.get()
        self.varfit.set(self.data_dir + '/analysis/')
        self.varerr.set(self.data_dir + '/analysis/')
        tk.Label(self.root, text='Saving Fit Values to:').place(x=0, y=255)
        fit_save_ent = tk.Entry(self.root, textvariable=self.varfit, width=44)
        fit_save_ent.place(x=0, y=275)


        save_figs_q = [('Yes', 'yes'), ('No', 'no')]
        x_init = 50
        tk.Label(self.root, text='Save Figures?:').place(x=0, y=320)
        for text, val in save_figs_q:
            rb_sfig = tk.Radiobutton(self.root, text=text, value=val, variable=self.varsavefig, height=2, width=15, justify=tk.CENTER,
                                    indicatoron=0, selectcolor='green', command=self.disable_figs)
            rb_sfig.place(x=x_init, y=340)
            x_init = x_init + 120

        self.data_dir = self.vardata.get()
        self.varfig.set(self.data_dir + '/analysis/figures/')
        tk.Label(self.root, text='Saving Figures to:').place(x=0, y=380)
        fig_save_ent = tk.Entry(self.root, textvariable=self.varfig, width=44)
        fig_save_ent.place(x=0, y=400)


        button_state_val = tk.Button(self.root, text="Enter", command=self.validate_all, height=2, width=10, justify=tk.CENTER)
        button_state_val.place(x=60, y=440)
        button_cancel = tk.Button(self.root, text="Cancel", command=self.end_program, height=2, width=10, justify=tk.CENTER)
        button_cancel.place(x=180, y=440)

        self.root.mainloop()

    # Definition for selecting directory
    def dir_dialog(self):
        self.vardata.set(askdirectory(initialdir='~/Documents/Francium'))
        self.data_dir = self.vardata.get()
        self.varfit.set(self.data_dir + '/analysis/')
        self.varerr.set(self.data_dir + '/analysis/')
        if self.save_figq == 'no':
            self.varfig.set(self.data_dir + '/analysis/figure/')
        else:
            self.varfig.set(self.data_dir + '/analysis/figure/' + self.figr_name + '.pdf')

        if self.fit_name != '':
            self.varfit.set(self.data_dir + '/analysis/' + self.fit_name+'.txt')
            self.varerr.set(self.data_dir + '/analysis/' + self.err_name + '.txt')
            if self.save_figq == 'no':
                self.varfig.set(self.data_dir + '/analysis/figure/')
            elif self.save_figq == 'yes':
                self.varfig.set(self.data_dir + '/analysis/figure/' + self.figr_name + '.pdf')


    def get_vals(self):
        self.fit_select = self.var1.get()
        self.scan_direc_select = self.var2.get()
        self.summing_scan = self.var3.get()
        self.volt_select = self.var4.get()

        if self.volt_select == '' and self.scan_direc_select == 'none':
            self.fit_name = '%s_%s_fit' % (self.fit_select, self.summing_scan)
            self.err_name = '%s_%s_err' % (self.fit_select, self.summing_scan)
            self.figr_name = '%s_%s_figure' % (self.fit_select, self.summing_scan)
        elif self.volt_select == '':
            self.fit_name = '%s_%s_%s_fit' % (self.fit_select, self.scan_direc_select, self.summing_scan)
            self.err_name = '%s_%s_%s_err' % (self.fit_select, self.scan_direc_select, self.summing_scan)
            self.figr_name = '%s_%s_%s_figure' % (self.fit_select, self.scan_direc_select, self.summing_scan)
        elif self.scan_direc_select == 'None':
            volt_str_reformat = self.volt_select.replace(' ', '_')
            self.fit_name = '%s_%s_%s_fit' % (self.fit_select, volt_str_reformat, self.summing_scan)
            self.err_name = '%s_%s_%s_err' % (self.fit_select, volt_str_reformat, self.summing_scan)
            self.figr_name = '%s_%s_%s_figure' % (self.fit_select, volt_str_reformat, self.summing_scan)
        else:
            volt_str_reformat = self.volt_select.replace(' ', '_')
            volt_str_reformat = volt_str_reformat.replace('.', '_0')
            self.fit_name = '%s_%s_%s_%s_fit' % (self.fit_select, volt_str_reformat, self.scan_direc_select, self.summing_scan)
            self.err_name = '%s_%s_%s_%s_err' % (self.fit_select, volt_str_reformat, self.scan_direc_select, self.summing_scan)
            self.figr_name = '%s_%s_%s_%s_figure' % (self.fit_select, volt_str_reformat, self.scan_direc_select,
                                                     self.summing_scan)

        self.varfit.set(self.data_dir + '/analysis/' + self.fit_name + '.txt')
        self.varerr.set(self.data_dir + '/analysis/' + self.err_name + '.txt')
        if self.save_figq == 'no':
            self.varfig.set(self.data_dir + '/analysis/figure/')
        elif self.save_figq == 'yes':
            self.varfig.set(self.data_dir + '/analysis/figure/' + self.figr_name + '.pdf')

        print(self.fit_select, self.scan_direc_select, self.summing_scan, self.volt_select)

    def validate_all(self):
        self.fit_select = self.var1.get()
        self.scan_direc_select = self.var2.get()
        self.summing_scan = self.var3.get()
        self.volt_select = self.var4.get()
        self.data_dir = self.vardata.get()
        self.fit_path = self.varfit.get()
        self.err_path = self.varerr.get()
        self.figr_path = self.varfig.get()
        self.root.quit()


    def disable_figs(self):
        self.save_figq = self.varsavefig.get()
        if self.save_figq == 'no':
            fig_save_ent = tk.Entry(self.root, textvariable=self.varfig, width=44, state=tk.DISABLED)
            fig_save_ent.place(x=0, y=400)
        elif self.save_figq == 'yes':
            fig_save_ent = tk.Entry(self.root, textvariable=self.varfig, width=44)
            fig_save_ent.place(x=0, y=400)
            self.varfig.set(self.data_dir + '/analysis/figure/' + self.figr_name + '.pdf')


    def end_program(self):
        print('Terminating program. Goodbye!')
        run_program = False
        exit()
run_program = True
user_states = GUI()

while run_program is True:
    # GUI section
    user_states.gui_loop()

    # Grabs the directory for the data
    path = user_states.data_dir
    file_list = glob.glob(path + '/*.txt')
    file_list = sorted(file_list)[0:-1]

    bin_num = 79           
    # binwidth in ms
    b_width = 0.2
    first_p = 0
    last_p = 79
    pps = last_p - first_p

    # Dataframe for storing fit values

    fit_vals = pd.DataFrame(columns=['Time', 'Offset', 'Peak Height', 'FWHM_L', 'FWHM_G', 'Peak Position', 'Lifetime',
                                    'Reduced Chi Squared', 'Offset err', 'Peak Height err', 'FWHM_L err', 'FWHM_G err',
                                    'Peak Position err', 'Lifetime err'])


    sum_data = pd.DataFrame(columns=['Counts'])
    data_raw = np.array([])
    data_time = []
    n_of_scans = len(file_list)

    # get time when file was created, this is important for looking at indvidual scans
    file_ls = GetTimes()
    file_ls.get_time(file_list)
    acq_time = file_ls.list_files

    # Prep for saving fit data
    if not os.path.exists(os.path.dirname(user_states.fit_path)):
        os.makedirs(os.path.dirname(user_states.fit_path))
        f1 = open((user_states.fit_name+'.txt'), 'a')
        f1.close()

    if not os.path.exists(os.path.dirname(user_states.err_path)):
        os.makedirs(os.path.dirname(user_states.err_path))
        f2 = open((user_states.err_name+'.txt'), 'a')
        f2.close()
    if user_states.save_figq == 'yes':
        if not os.path.exists(os.path.dirname(user_states.figr_path)):
            os.makedirs(os.path.dirname(user_states.figr_path))
            f3 = open((user_states.figr_name+'.txt'), 'a')
            f3.close()



    # begin minimization
    print("Starting Minimization...")
    data_for_plotting = []

    # if user selects to look at each scan individually
    if user_states.summing_scan == 'individual':
        # Use file name to get data from file, the files should be chronologically ordered at this point
        for index, row in acq_time.iterrows():
            times = row['Time']
            file_path = row['File path']
            data_raw = (np.genfromtxt(file_path))
            data_binned = np.array([])
            a = []

            # Gather data points an bin them
            for i in range(0, 160):
                y_bin = (data_raw[i * bin_num * 1:(i * bin_num* 1) + bin_num * 1])
                a = sum(y_bin[first_p:last_p])
                data_binned = np.append(data_binned, a)

            x = np.arange(0, len(data_binned), 1)

            ''' Now to minimize the scans'''
            data = DataManager(user_states.scan_direc_select)
            data.data_select(data.name, data_binned, x, b_width, pps)
            y_data = np.array(data.y_val)
            x_step = np.array(data.x_val)
            y_err = np.array(data.y_err)
            time_stamp = times
            if user_states.fit_select != 'none':
                print("Starting Minimization...")
                p1 = Parameters(y_data, x_step, user_states.fit_select)  # Generate initial parameter guesses for scans
                val = MinimizationRoutine(user_states.fit_select)
                val.minimzation(val.name, p1)
                p_fit = val.fit
                p_err = val.err
                x_fit = np.arange(min(x_step), max(x_step) + 1, 0.1)
                if val.name == 'lorentzian':
                    exp_val = lorentzian(x_step, p_fit)
                    y_fit = lorentzian(x_fit, p_fit)
                    Red_chi2 = chi2(p_fit) / (len(y_data) - len(p_fit))
                    Counting_Rate = lorentzian(p_fit[3], p_fit)
                    p_fit.insert(3, 'N/A')
                    p_err.insert(3, 'N/A')

                elif val.name == 'voigt':
                    exp_val = voigt(x_step, p_fit)
                    y_fit = voigt(x_fit, p_fit)
                    Red_chi2 = chi2(p_fit) / (len(y_data) - len(p_fit))
                    p_fit.insert(6, 'N/A')
                    p_err.insert(6, 'N/A')

                repack_data = pd.DataFrame({
                    'rate': y_data,
                    'err': y_err,
                    'Normalized Residuals': (y_data - exp_val) / y_err,
                    'Residuals': (y_data - exp_val),
                    'Freq': x_step,
                })

                fit_vals = fit_vals.append({
                    'Time': time_stamp,
                    'Offset': p_fit[0],
                    'Peak Height': p_fit[1],
                    'FWHM_L': p_fit[2],
                    'FWHM_G': p_fit[3],
                    'Peak Position': p_fit[4],
                    'Lifetime': p_fit[5],
                    'Reduced Chi Squared': Red_chi2,
                    'Offset err': p_err[0],
                    'Peak Height err': p_err[1],
                    'FWHM_L err': p_err[2],
                    'FWHM_G err': p_err[3],
                    'Peak Position err': p_err[4],
                    'Lifetime err': p_err[5]
                }, ignore_index=True)

                fit = pd.DataFrame({
                    'x_fit': x_fit,
                    'y_fit': y_fit
                })

                data_for_plotting.append([repack_data, fit])

                print("Minimization complete")
                # save the data and errors to csv files
                fit_vals.to_csv(user_states.fit_path, sep='\t', index=False)
            elif user_states.fit_select == 'none':
                repack_data = pd.DataFrame({
                    'rate': y_data,
                    'err': y_err,
                    'Freq': x_step,
                })
                data_for_plotting.append(repack_data)


    # if user selects to sum all the scans together
    elif user_states.summing_scan == 'sum':
        start_end_time = '%s-%s' % (acq_time['Time'].iloc[0], acq_time['Time'].iloc[-1])
        data0 = []
        # Use file name to get data from file, the files are chronologically ordered
        for index, row in acq_time.iterrows():
            file_path = row['File path']
            data0.append(np.genfromtxt(file_path))
        shape = list(data0[0].shape)
        shape[:0] = [len(data0)]
        data1 = np.concatenate(data0).reshape(shape)
        y = data1.sum(axis=0)
        data_binned = np.array([])
        a = []
        # Gather data points an bin them
        for i in range(0, 160):
            y_bin = (y[i * bin_num * 1:(i * bin_num * 1) + bin_num * 1])
            a = sum(y_bin[first_p:last_p])
            data_binned = np.append(data_binned, a)

        x = np.arange(0, len(data_binned), 1)
        ''' Now to minimize the scans'''
        data = DataManager(user_states.scan_direc_select)
        data.data_select(data.name, data_binned, x, b_width, pps)
        y_data = np.array(data.y_val) / n_of_scans
        x_step = np.array(data.x_val)
        y_err = np.array(data.y_err) / n_of_scans
        if user_states.fit_select != 'none':
            print("Starting Minimization...")
            p1 = Parameters(y_data, x_step, user_states.fit_select)  # Generate initial parameter guesses for scans
            val = MinimizationRoutine(user_states.fit_select)
            val.minimzation(val.name, p1)
            p_fit = val.fit
            p_err = val.err
            x_fit = np.arange(min(x_step), max(x_step) + 1, 0.1)

            if val.name == 'lorentzian':
                exp_val = lorentzian(x_step, p_fit)
                y_fit = lorentzian(x_fit, p_fit)
                Red_chi2 = chi2(p_fit) / (len(y_data) - len(p_fit))
                Counting_Rate = lorentzian(p_fit[3], p_fit)
                crate_err = np.sqrt(Counting_Rate/(b_width*pps))
                p_fit.insert(3, 'N/A')
                p_err.insert(3, 'N/A')
                p_fit.insert(5, 'N/A')
                p_err.insert(5, 'N/A')

            elif val.name == 'voigt':
                exp_val = voigt(x_step, p_fit)
                y_fit = voigt(x_fit, p_fit)
                Red_chi2 = chi2(p_fit) / (len(y_data) - len(p_fit))
                p_fit.insert(6, 'N/A')
                p_err.insert(6, 'N/A')

            repack_data = pd.DataFrame({
                'rate': y_data,
                'err': y_err,
                'Normalized Residuals': (y_data - exp_val) / y_err,
                'Residuals': (y_data - exp_val),
                'Freq': x_step,

            })

            fit_vals = fit_vals.append({
                'Time': start_end_time,
                'Offset': p_fit[0],
                'Peak Height': p_fit[1],
                'FWHM_L': p_fit[2],
                'FWHM_G': p_fit[3],
                'Peak Position': p_fit[4],
                'Lifetime': p_fit[5],
                'Reduced Chi Squared': Red_chi2,
                'Offset err': p_err[0],
                'Peak Height err': p_err[1],
                'FWHM_L err': p_err[2],
                'FWHM_G err': p_err[3],
                'Peak Position err': p_err[4],
                'Lifetime err': p_err[5],
                'Counting Rate (kHz)': Counting_Rate,
                'rate err': crate_err
            }, ignore_index=True)

            fit = pd.DataFrame({
                'x_fit': x_fit,
                'y_fit': y_fit
            })
            y_fit, ycov = propagate(lambda p: lorentzian(x_fit, p), val.mvalues, val.covar)
            yerr_prop = np.diag(ycov) ** 0.5
            data_for_plotting.append([repack_data, fit])
            print("Minimization complete")
            fit_vals.to_csv(user_states.fit_path, sep='\t', index=False)
        elif user_states.fit_select == 'none':
            repack_data = pd.DataFrame({
                'rate': y_data,
                'err': y_err,
                'Freq': x_step,

            })
            data_for_plotting.append(repack_data)





    # Plotting routine
    #matplotlib.use('Cairo')
    plt.style.use('../matplotlib_style/stylelib/cern_root.mplstyle')

    if user_states.save_figq == 'no':
        pdf = None
    elif user_states.save_figq == 'yes':
        pdf = PdfPages(user_states.figr_path)
    q = 1
    l = 0
    for data in data_for_plotting:
        if user_states.fit_select != 'none':
            fig = plt.figure(figsize=(8,6))
            gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
        
            fig.subplots_adjust(hspace=0)
            ax1 = fig.add_subplot(gs[0, 0])

            
            ax1.set_ylabel(r'\textbf{Counting Rate (kHz)}', fontsize=16)
            ax2 = fig.add_subplot(gs[1, 0])
            ax2.set_xlabel(r'\textbf{Laser detuning (MHz)}', fontsize=16)
            ax2.set_ylabel(r'\textbf{Normalized}'+'\n'+r'\textbf{Residuals}', fontsize=16)
            ax1.sharex(ax2)
            ax1.set_xlabel('', fontsize=16)
            ax2.set_ylim(-3,3)

            ax1.fill_between(x_fit, y_fit - yerr_prop,  y_fit + yerr_prop, facecolor="b", alpha=0.35)
            ax1.errorbar(x_step, y_data, yerr=y_err, fmt='ro', ecolor='black', capsize=5, zorder=1)

            ax1.plot(x_fit, y_fit, color='b', zorder=2, label='Lorentzain')

            tick_div = (ax1.get_yticks()[1] - ax1.get_yticks()[0]) / 5
            ax2.plot(x_step, data[0]['Normalized Residuals'], 'o', color='b')
            ax2.axhline(y=0, color='black')
            plt.show()

        elif user_states.fit_select == 'none':
            fig = plt.figure()
            plt.xlabel('Time (s)', fontsize=16)
            plt.ylabel('Rate (kHz)', fontsize=16)
            plt.errorbar(data['Freq'], data['rate'], yerr=data['err'], fmt='r-o', ecolor='black', capsize=5, zorder=1)
            if user_states.summing_scan == 'individual':
                if user_states.volt_select =='':
                    plt.title(r'%s %s Scan \#%i' % (transition_type, user_states.scan_direc_select, q), fontsize=20)
                else:
                    plt.title('%s %s Scan \#%i @ %s' % (transition_type, user_states.scan_direc_select, q,
                                                        user_states.volt_select), fontsize=20)
            elif user_states.summing_scan == 'sum':
                if user_states.volt_select == '':
                    plt.title(r'%s %i %s Scan' % (transition_type, n_of_scans,user_states.scan_direc_select), fontsize=20)
                else:
                    plt.title('%s %i %s Scan @ %s' % (transition_type, n_of_scans, user_states.scan_direc_select, user_states.volt_select),
                              fontsize=20)
                plt.title("")
                plt.show()


        if pdf is not None:
            pdf.savefig(fig)
        plt.close(fig)

        q = q + 1
        l = l + 1
    if pdf is not None:
        pdf.close()
    print("Graph(s) completed")


