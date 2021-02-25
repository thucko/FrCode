'''
- Using Matplotlib for subplots of data and residuals
- plots now contain normalized residuals and reqular residuals
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
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
from scipy.special import wofz
import tkinter as tk
import matplotlib
import matplotlib.pyplot as plt
import warnings
import time
import pyautogui
import os
warnings.simplefilter(action='ignore', category=FutureWarning)

# Lorentzian function
def lorentzian(x, p_lorentzian):
    # Here we use:(gamma/2)^2/((x-x0)^2+(gamma/2)^2)
    # Including and exponetial decay
    numerator = (p_lorentzian[2] / 2) ** 2
    decay = (np.exp(-x / p_lorentzian[4]))
    denominator1 = ((x - (p_lorentzian[3])) ** 2 + (p_lorentzian[2] / 2) ** 2)
    L = p_lorentzian[0] + decay * p_lorentzian[1] * (numerator / denominator1)
    return L


# Voigt profile
# see https://en.wikipedia.org/wiki/Voigt_profile for info on the Voigt profile
def voigt(x, p):
    sigma = p[3] / (2 * np.sqrt(2 * np.log(2)))
    gamma = p[2] / 2
    z = (x - p[4] + 1j * gamma) / (sigma * np.sqrt(2))
    decay = (np.exp(-x / p[5]))
    num = np.real(wofz(z))
    dem = sigma * np.sqrt(2 * np.pi)
    V = p[0] + decay * p[1] * (num / dem)
    return V


# Parameters for minimization includes Voigt and Lorentzian parameters
def Parameters(data, x, f):
    # Lorentzain: p = [Offset, Height, FWHM, Center, lifetime]
    # Voigt: p = [Offset, Height, FWHM_L, FWHM_G, Center, Liftime]
    if f == 'lorentzian':
        p = np.array((np.min(data), np.max(data), 3.5, x[np.argmax(data)], 350))
        return p
    elif f == 'voigt':

        p = np.array((np.min(data), np.max(data), 3.6, 1, x[np.argmax(data)], 350))
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

    def data_select(self, name, data, x):

        if name == 'forward':
            a = 0
            b = 80
        elif name == 'backward':
            a = 80
            b = 160
        elif name == 'none':  # New option in case we sit at one frequency
            a = 0
            b = 160

        self.y_val.extend(data[a:b] / (b_width * pps))  # converts counts into rate (Hz)
        self.x_val.extend(x[a:b])
        self.y_err.extend(np.sqrt(data[a:b]) / (b_width * pps))  # transforms the uncertainty into rate


# minimizaion class using iminuit
class MinimizationRoutine:

    def __init__(self, name):
        self.name = name
        self.fit = []
        self.err = []

    def minimzation(self, name, par):
        if name == 'voigt':
            errors = (10, 10, 0.01, 0.001, .1, 1)
            limits = (None, None, (3, 6), (0, 0.9), None, (17, 500))
            fixed = (False, False, False, False, False, False)

        elif name == 'lorentzian':
            errors = (10, 10, 0.01, .1, 1)
            limits = (None, None, (2, 6), None, None)
            fixed = (False, False, False, False, False)

        m = Minuit.from_array_func(chi2, par, error=errors, limit=limits, fix=fixed, errordef=1, pedantic=False)
        m.migrad()

        if name == 'voigt':
            #               [Offset        , Height        , FWHM_L        , FWHM_G       , Center         , liftime       ]
            self.fit.extend([m.values["x0"], m.values["x1"], m.values["x2"], m.values["x3"], m.values["x4"], m.values["x5"]])
            self.err.extend([m.errors["x0"], m.errors["x1"], m.errors["x2"], m.errors["x3"], m.errors["x4"], m.errors["x5"]])

        elif name == 'lorentzian':
            #               [Offset        , Height        , FWHM_L        , Center        , lifetime      ]
            self.fit.extend([m.values["x0"], m.values["x1"], m.values["x2"], m.values["x3"], m.values["x4"]])
            self.err.extend([m.errors["x0"], m.errors["x1"], m.errors["x2"], m.errors["x3"], m.errors["x4"]])


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
        self.vardata.set('~/Documents')
        tk.Label(self.root, text='Select Data Directory:').place(x=0, y=5)
        dir_path = tk.Entry(self.root, textvariable=self.vardata, width=44)
        dir_path.place(x=0, y=25)
        img = tk.PhotoImage(file='../ui_icons/folder.png')
        img = img.subsample(100)
        dir_button = tk.Button(self.root, image=img, command=self.dir_dialog, width=25, height=20)
        dir_button.place(x=316, y=20)

        # State select buttons
        funcs = [('Lorentzian', 'lorentzian'), ('Voigt', 'voigt')]
        scan_direction = [('Forward Scans', 'forward'), ('Backward Scans', 'backward'), ('None', 'none')]
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
        # button_cancel = tk.Button(self.root, text="Cancel", command=self.end_program, height=1, width=10)
        # button_cancel.place(x=252, y=145)

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
        self.vardata.set(askdirectory(initialdir='~/Documents'))
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
            else:
                print('Do you want to save fiures?')

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
        else:
            print('Do you want to save fiures?')

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
        else:
            print('Do you want to save fiures?')

    def end_program(self):
        print('Terminating program. Goodbye!')
        exit()

# GUI section
user_states = GUI()
user_states.gui_loop()

# Grabs the directory for the data
path = user_states.data_dir
file_list = glob.glob(path + '/*.txt')

# binwidth in ms
b_width = 2.8
first_p = 2
last_p = 8
pps = last_p - first_p

# Dataframe for storing fit values
'''fit_vals = pd.DataFrame(columns=['Time', 'Offset', 'Peak Height', 'FWHM_L', 'FWHM_G', 'Peak Position', 'Lifetime',
                                 'Reduced Chi Squared'])
error_vals = pd.DataFrame(columns=['Time', 'Offset err', 'Peak Height err', 'FWHM_L err', 'FWHM_G err',
                                   'Peak Position err', 'Lifetime err'])'''
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
            y_bin = (data_raw[i * 8 * 1:(i * 8 * 1) + 8 * 1])
            a = sum(y_bin[first_p:last_p])
            data_binned = np.append(data_binned, a)

        x = np.arange(0, len(data_binned), 1)

        ''' Now to minimize the scans'''
        data = DataManager(user_states.scan_direc_select)
        data.data_select(data.name, data_binned, x)
        y_data = np.array(data.y_val)
        x_step = np.array(data.x_val)
        y_err = np.array(data.y_err)
        time_stamp = times
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
            p_fit.insert(3, 'N/A')
            p_err.insert(3, 'N/A')

        elif val.name == 'voigt':
            exp_val = voigt(x_step, p_fit)
            y_fit = voigt(x_fit, p_fit)
            Red_chi2 = chi2(p_fit) / (len(y_data) - len(p_fit))

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
        y_bin = (y[i * 8 * 1:(i * 8 * 1) + 8 * 1])
        a = sum(y_bin[first_p:last_p])
        data_binned = np.append(data_binned, a)

    x = np.arange(0, len(data_binned), 1)
    ''' Now to minimize the scans'''
    data = DataManager(user_states.scan_direc_select)
    data.data_select(data.name, data_binned, x)
    y_data = np.array(data.y_val) / len(data0)
    x_step = np.array(data.x_val)
    y_err = np.array(data.y_err) / len(data0)
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
        p_fit.insert(3, 'N/A')
        p_err.insert(3, 'N/A')

    elif val.name == 'voigt':
        exp_val = voigt(x_step, p_fit)
        y_fit = voigt(x_fit, p_fit)
        Red_chi2 = chi2(p_fit) / (len(y_data) - len(p_fit))

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
        'Lifetime err': p_err[5]
    }, ignore_index=True)

    fit = pd.DataFrame({
        'x_fit': x_fit,
        'y_fit': y_fit
    })
    data_for_plotting.append([repack_data, fit])

    fit_vals.to_csv(user_states.fit_path, sep='\t', index=False)


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

    gs = gridspec.GridSpec(2, 2)
    fig = plt.figure(q - 1, figsize=(16, 9))
    fig.subplots_adjust(left=0.06, bottom=0.08, right=0.95, top=0.95, wspace=0.07, hspace=0.25)
    ax1 = fig.add_subplot(gs[0, :])
    ax1.set_xlabel('Frequency (MHz)', fontsize=16)
    ax1.set_ylabel('Rate (kHz)', fontsize=16)
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.set_title('Normalized Residuals', fontsize=16)
    ax2.set_xlabel('Frequency (MHz)', fontsize=14)
    ax3 = fig.add_subplot(gs[1, 1])
    ax3.set_title('Residuals', fontsize=16)
    ax3.set_xlabel('Frequency (MHz)', fontsize=14)

    ax1.text(data[0]['Freq'].iloc[0] + 20, np.max(data[0]['rate']) - 0.1 * np.max(data[0]['rate']),
             r'$\chi^2_{reduced}=%.3f$' % fit_vals['Reduced Chi Squared'].iloc[l], fontsize=14)
    ax1.text(data[0]['Freq'].iloc[0] + 20, np.max(data[0]['rate']) - 0.17 * np.max(data[0]['rate']),
             r'$\nu_{peak}=%.3f$ MHz' % fit_vals['Peak Position'].iloc[l], fontsize=14)
    ax1.text(data[0]['Freq'].iloc[0] + 20, np.max(data[0]['rate']) - 0.24 * np.max(data[0]['rate']),
             r'$\delta\nu_{peak}=%.3f$ MHz' % fit_vals['Peak Position err'].iloc[l], fontsize=14)

    if user_states.summing_scan == 'individual':
        if user_states.volt_select =='':
            ax1.set_title(r'%s Scan \#%i (%s fit) ' % (user_states.scan_direc_select, q, user_states.fit_select),
                          fontsize=20)
        else:
            ax1.set_title('%s Scan \#%i @ %s (%s fit) ' % (user_states.scan_direc_select, q, user_states.volt_select,
                                                           user_states.fit_select), fontsize=20)
    elif user_states.summing_scan == 'sum':
        if user_states.volt_select is None:
            ax1.set_title(r'%s Scan (%s fit) ' % (user_states.scan_direc_select, user_states.fit_select), fontsize=20)
        else:
            ax1.set_title('%s Scan @ %s (%s fit) ' % (user_states.scan_direc_select, user_states.volt_select,
                                                      user_states.fit_select), fontsize=20)

    ax1.errorbar(data[0]['Freq'], data[0]['rate'], yerr=data[0]['err'], fmt='ro', ecolor='black',
                 capsize=5, zorder=1)
    ax1.plot(data[1]['x_fit'], data[1]['y_fit'], color='b', zorder=2)
    ax2.plot(data[0]['Freq'], data[0]['Normalized Residuals'], '-o', color='#F06C00')
    ax3.plot(data[0]['Freq'], data[0]['Residuals'], '-o', color='#D433FF')
    if user_states.summing_scan == 'sum':
        plt.show()

    if pdf is not None:
        pdf.savefig(fig)
    plt.close(fig)

    q = q + 1
    l = l + 1
if pdf is not None:
    pdf.close()

