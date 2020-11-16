'''
- Using Matplotlib for subplots of data and residuals
- plots now contain normalized residuals and reqular residuals
-
Author: Tim Hucko

'''

import numpy as np
import pandas as pd
from tkinter import simpledialog
from tkinter.filedialog import askdirectory
from tkinter.filedialog import askopenfilename
from tkinter.filedialog import asksaveasfile
from tkinter import messagebox
import glob
from iminuit import Minuit
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
from scipy.special import wofz
from tkinter import *
import matplotlib.pyplot as plt
import warnings
import time


warnings.simplefilter(action='ignore', category=FutureWarning)


# Definition for selecting directory
def dir_dialog():
    Tk().withdraw()
    pwd = askdirectory(initialdir='~/Documents')
    return pwd


# Definition for saving files
def file_save():
    f = asksaveasfile(mode='w', defaultextension=".txt", initialfile="%s_%s_%s_%s_data" % (volt, sel2, sel, sel3))
    if f is None:  # asksaveasfile return `None` if dialog closed with "cancel".
        return
    return f


# Definition for saving uncertainties
def file_save_err():
    f = asksaveasfile(mode='w', defaultextension=".txt", initialfile="%s_%s_%s_%s_err" % (volt, sel2, sel, sel3))
    if f is None:  # asksaveasfile return `None` if dialog closed with "cancel".
        return
    return f


# Definition for saving figures
def file_save_plt():
    f = asksaveasfile(mode='wb', defaultextension=".pdf", title='Save Figure(s)',
                      initialfile='%s_scan_%s_fit_%s_%s' % (sel2, sel, sel3, volt))
    if f is None:  # asksaveasfile return `None` if dialog closed with "cancel".
        return
    return f



# Definition for choosing a file
def file_dialog():
    Tk().withdraw()
    file = askopenfilename()
    return file


# Definition for choosing to sum the scans or look at individual scans
def sum_scan():
    global sel3
    MsgBox = messagebox.askquestion('Sum Spectra', 'Would you like to sum all the spectra together? (Selecting "No" will'
                                                   ' look at each file indvidually)')
    if MsgBox == 'yes':
        messagebox.showinfo('Sum Spectra', 'Summing the Spectra!')
        sel3 = 'summed'
        return True
    else:
        messagebox.showinfo('Sum Spectra', 'Looking at scans indvidually!')
        sel3 = 'Indviduals'
        return False


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
        z = (x - p[4] +1j*gamma) / (sigma * np.sqrt(2))
        decay = (np.exp(-x / p[5]))
        num = np.real(wofz(z))
        dem = sigma * np.sqrt(2 * np.pi)
        V = p[0] + decay * p[1] * (num / dem)
        return V


# Parameters for minimization includes Voigt and Lorentzian parameters
def Parameters(data, x, f):
    # Lorentzain: p = [Offset, Height, FWHM, Center, lifetime]
    # Voigt: p = [Offset, Height, FWHM_L, FWHM_G, Center, Liftime]
    if f == 'Lorentzian':
        p = (np.min(data), np.max(data), 3.5, x[np.argmax(data)], 350)
        return p
    elif f == 'Voigt':

         p = (np.min(data), np.max(data), 3.6, 1, x[np.argmax(data)], 350)
        # p = [np.min(data), np.max(data), 3.6, 1, device_put(x)[np.argmax(data)], 350]

         return p


# Chi squared function used in the minimization

def chi2(p):
    if sel == 'Lorentzian':
        func = lorentzian(x_step, p)  # Our function we need
    elif sel == 'Voigt':
        func = voigt(x_step, p)
    err = y_err  # Uncertainty in the data sqrt(N)
    # calculate the chi2
    delta = (y_data - func) / err
    chi2_v = sum(pow(delta, 2))
    return chi2_v



# Used for selecting Lorentzian or Voigt fitting fucntions
def select_fun():
    global sel
    sel = var.get()
    roo.quit()
    roo.destroy()


def select_FB():
    global sel2
    sel2 = var2.get()
    roo2.quit()
    roo2.destroy()


class DataManager:
    def __init__(self, name):
        self.name = name
        self.y_val = []
        self.x_val = []
        self.y_err = []
    def data_select(self, name, data, x):

        if name == 'Forward':
            a = 0
            b = 80
        elif name == 'Backward':
            a = 80
            b = 160

        self.y_val.extend(data[a:b]/(b_width*pps))
        self.x_val.extend(x[a:b])
        self.y_err.extend(np.sqrt(data[a:b])/(b_width*pps))


class MinimizationRoutine:

    def __init__(self, name):
        self.name = name
        self.fit = []
        self.err = []
    def minimzation(self, name, par):
        if name == 'Voigt':
            errors = (10, 10, 0.01, 0.001, .1, 1)
            limits = (None, None, (3, 6), (0, 0.9), None, (17, 500))
            fixed = (False, False, False, False, False, False)

        elif name == 'Lorentzian':
            errors = (10, 10, 0.01, .1, 1)
            limits = (None, None, (2, 6), None, None)
            fixed = (False, False, False, False, False)

        m = Minuit.from_array_func(chi2, par, error=errors, limit=limits, fix=fixed, errordef=1, pedantic=False)
        m.migrad()

        if name == 'Voigt':
            #               [Offset        , Height        , FWHM_L        , FWHM_G       , Center         , liftime       ]
            self.fit.extend([m.values["x0"], m.values["x1"], m.values["x2"], m.values["x3"], m.values["x4"], m.values["x5"]])
            self.err.extend([m.errors["x0"], m.errors["x1"], m.errors["x2"], m.errors["x3"], m.errors["x4"], m.errors["x5"]])

        elif name == 'Lorentzian':
            #               [Offset        , Height        , FWHM_L        , Center        , lifetime      ]
            self.fit.extend([m.values["x0"], m.values["x1"], m.values["x2"], m.values["x3"], m.values["x4"]])
            self.err.extend([m.errors["x0"], m.errors["x1"], m.errors["x2"],  m.errors["x3"], m.errors["x4"]])


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




# Grabs the directory for the data
path = dir_dialog()
file_list = glob.glob(path + '/*.txt')
#list_file = pd.DataFrame(columns=['Time', 'File path'])

# GUI for selecting fitting functions
funcs = ["Lorentzian", "Voigt"]
roo = Tk()
roo.geometry('260x75')
roo.resizable(0, 0)
roo.title('Fitting Function')
var = StringVar(roo)
var.set(funcs[1])
w = OptionMenu(roo, var, *funcs)
w.pack()
button = Button(roo, text="Select", command=select_fun)
button.pack()
mainloop()

# GUI for selecting forward or backwards scans
direction = ['Forward', 'Backward']
roo2 = Tk()
roo2.geometry('260x75')
roo2.resizable(0, 0)
roo2.title('Scan Selection')
var2 = StringVar(roo2)
var2.set(direction[0])
w2 = OptionMenu(roo2, var2, *direction)
w2.pack()
button2 = Button(roo2, text="Select", command=select_FB)
button2.pack()
mainloop()


# binwidth in ms
b_width = 2.8
first_p = 2
last_p = 8
pps = last_p - first_p

fit_vals = pd.DataFrame(columns=['Time', 'Offset', 'Peak Height', 'FWHM_L', 'FWHM_G', 'Peak Position', 'Lifetime',
                                     'Reduced Chi Squared'])
error_vals = pd.DataFrame(columns=['Time', 'Offset err', 'Peak Height err', 'FWHM_L err', 'FWHM_G err',
                               'Peak Position err', 'Lifetime err'])

sum_data = pd.DataFrame(columns=['Counts'])
data_raw = np.array([])
data_time = []
n_of_scans = len(file_list)

# Ask user if they would like to save the individual spectra
'''fig_save = save_fig()
if fig_save is True:
    print('Choose Directory...')
    fig_pwd = dir_dialog()
if fig_save is False:
    print('Not saving spectra...')'''


# get time when file was created, this is important for looking at indvidual scans

file_ls = GetTimes()
file_ls.get_time(file_list)
acq_time = file_ls.list_files
prompt1 = sum_scan()
root = Tk()
root.withdraw()
# Ask user for Voltage
volt = simpledialog.askstring(title='Voltage', prompt='Voltage used:')


# begin minimization
print("Starting Minimization...")
begin = time.time()
# if user selects to look at each scan individually
if prompt1 is False:
    q = 1
    fig_save = file_save_plt()
    if fig_save is None:
        None
        pdf = None
    else:
        pdf = PdfPages(fig_save)

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
        data = DataManager(sel2)
        data.data_select(data.name, data_binned, x)
        y_data = np.array(data.y_val)
        x_step = np.array(data.x_val)
        y_err = np.array(data.y_err)
        time_stamp = times
        p1 = Parameters(y_data, x_step, sel)  # Generate initial parameter guesses for scans
        val = MinimizationRoutine(sel)
        val.minimzation(val.name, p1)
        p_fit = val.fit
        p_err = val.err
        x_fit = np.arange(min(x_step), max(x_step) + 1, 0.1)
        if val.name == 'Lorentzian':
            exp_val = lorentzian(x_step, p_fit)
            y_fit = lorentzian(x_fit, p_fit)
            Red_chi2 = chi2(p_fit) / (len(y_data) - len(p_fit))
            p_fit.insert(3, 'N/A')
            p_err.insert(3, 'N/A')

        elif val.name == 'Voigt':
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
            'Reduced Chi Squared': Red_chi2
        }, ignore_index=True)

        error_vals = error_vals.append({
            'Time': time_stamp,
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
        plt.style.use('ggplot')
        gs = gridspec.GridSpec(2, 2)
        fig = plt.figure(q-1, figsize=(16, 9))

        ax1 = fig.add_subplot(gs[0, :])
        ax1.text(x_fit[0] + 20, np.max(y_data) - 0.1 * np.max(y_data), r'$\chi^2_{reduced}=%.3f$' % Red_chi2,
                 fontsize=14)
        ax1.text(x_fit[0] + 20, np.max(y_data) - 0.17 * np.max(y_data), r'$\nu_{peak}=%.3f$ MHz' % p_fit[4],
                 fontsize=14)
        ax1.text(x_fit[0] + 20, np.max(y_data) - 0.24 * np.max(y_data), r'$\delta\nu_{peak}=%.3f$ MHz' % p_err[4],
                 fontsize=14)
        ax1.set_title('%s Scan #%d @ %s (%s fit) ' % (sel2, q, volt, sel), fontsize=20)
        ax1.set_xlabel('Frequency (MHz)', fontsize=16)
        ax1.set_ylabel('Rate (kHz)', fontsize=16)
        ax1.errorbar(repack_data['Freq'], repack_data['rate'], yerr=repack_data['err'], fmt='ro', ecolor='black',
                     capsize=5, zorder=1)
        ax1.plot(fit['x_fit'], fit['y_fit'], color='#3392FF', zorder=2)
        ax2 = fig.add_subplot(gs[1, 0])
        ax2.set_title('Normalized Residuals', fontsize=16)
        ax2.set_xlabel('Frequency (MHz)', fontsize=14)
        ax2.plot(repack_data['Freq'], repack_data['Normalized Residuals'], '-o', color='#F06C00')
        ax3 = fig.add_subplot(gs[1, 1])
        ax3.set_title('Residuals', fontsize=16)
        ax3.set_xlabel('Frequency (MHz)', fontsize=14)
        ax3.plot(repack_data['Freq'], repack_data['Residuals'], '-o', color='#D433FF')
        fig.subplots_adjust(left=0.06, bottom=0.08, right=0.95, top=0.95, wspace=0.07, hspace=0.25)
        # Saving figures
        if pdf is not None:
            pdf.savefig(fig)
        plt.close()
        q = q + 1
    if pdf is not None:
        pdf.close()
    fin = time.time()
    print("Minimization complete")
    print("Time: ", fin - begin, "s")
    # save the data and errors to csv files
    fit_vals.to_csv(file_save(), sep='\t', index=False)
    error_vals.to_csv(file_save_err(), sep='\t', index=False)

# if user selects to sum all the scans together
elif prompt1 is True:
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
    data = DataManager(sel2)
    data.data_select(data.name, data_binned, x)
    y_data = np.array(data.y_val)/len(data0)
    x_step = np.array(data.x_val)
    y_err = np.array(data.y_err)/len(data0)
    p1 = Parameters(y_data, x_step, sel)  # Generate initial parameter guesses for scans
    val = MinimizationRoutine(sel)
    val.minimzation(val.name, p1)
    p_fit = val.fit
    p_err = val.err
    x_fit = np.arange(min(x_step), max(x_step) + 1, 0.1)
    if val.name == 'Lorentzian':
        exp_val = lorentzian(x_step, p_fit)
        y_fit = lorentzian(x_fit, p_fit)
        Red_chi2 = chi2(p_fit) / (len(y_data) - len(p_fit))
        p_fit.insert(3, 'N/A')
        p_err.insert(3, 'N/A')

    elif val.name == 'Voigt':
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
        'Offset': p_fit[0],
        'Peak Height': p_fit[1],
        'FWHM_L': p_fit[2],
        'FWHM_G': p_fit[3],
        'Peak Position': p_fit[4],
        'Lifetime': p_fit[5],
        'Reduced Chi Squared': Red_chi2
    }, ignore_index=True)

    error_vals = error_vals.append({
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
    #df = pd.melt(repack_data, id_vars=['Freq'], value_vars=['Normalized Residuals', 'Residuals'], var_name='Type')
    plt.style.use('cern_root')
    '''gs = gridspec.GridSpec(2, 2)
    fig = plt.figure(figsize=(16, 9))

    ax1 = fig.add_subplot(gs[0, :])
    ax1.text(x_fit[0]+20, np.max(y_data)-0.1*np.max(y_data), r'$\chi^2_{reduced}=%.3f$' % Red_chi2, fontsize=14)
    ax1.text(x_fit[0] + 20, np.max(y_data) - 0.17 * np.max(y_data), r'$\nu_{peak}=%.3f$ MHz' % p_fit[4], fontsize=14)
    ax1.text(x_fit[0] + 20, np.max(y_data) - 0.24 * np.max(y_data), r'$\delta\nu_{peak}=%.3f$ MHz' % p_err[4], fontsize=14)
    ax1.set_title('%s Scan @ %s (%s fit) Summed %d Files' % (sel2, volt, sel, n_of_scans), fontsize=20)
    ax1.set_xlabel('Frequency (MHz)', fontsize=16)
    ax1.set_ylabel('Rate (kHz)', fontsize=16)
    ax1.errorbar(repack_data['Freq'], repack_data['rate'], yerr=repack_data['err'], fmt='ro', ecolor='black',
                 capsize=5, zorder=1)
    ax1.plot(fit['x_fit'], fit['y_fit'], color='b', zorder=2)
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.set_title('Normalized Residuals', fontsize=16)
    ax2.set_xlabel('Frequency (MHz)', fontsize=14)
    ax2.plot(repack_data['Freq'], repack_data['Normalized Residuals'], '-o', color='#F06C00')
    ax3 = fig.add_subplot(gs[1, 1])
    ax3.set_title('Residuals', fontsize=16)
    ax3.set_xlabel('Frequency (MHz)', fontsize=14)
    ax3.plot(repack_data['Freq'], repack_data['Residuals'], '-o', color='#D433FF')
    fig.subplots_adjust(left=0.06, bottom=0.08, right=0.95, top=0.95, wspace=0.07, hspace=0.25)'''

    #quick and dirty plotting
    plt.errorbar(repack_data['Freq'], repack_data['rate'], yerr=repack_data['err'], fmt='ro', ecolor='black',
                 capsize=5, zorder=1)
    plt.plot(fit['x_fit'], fit['y_fit'], color='b', zorder=2)
    #plt.text(x_fit[0]+20, np.max(y_data)-0.1*np.max(y_data), r'$\chi^2_{reduced}=%.3f$' % Red_chi2, fontsize=14)
    #plt.text(x_fit[0] + 20, np.max(y_data) - 0.17 * np.max(y_data), r'$\nu_{peak}=%.3f$ MHz' % p_fit[4], fontsize=14)
    #plt.text(x_fit[0] + 20, np.max(y_data) - 0.24 * np.max(y_data), r'$\delta\nu_{peak}=%.3f$ MHz' % p_err[4],
    #         fontsize=14)
    plt.title(r'$^{211}$Fr $\beta$ Transition', size=24)
    plt.xlabel('Frequency (MHz)', fontsize=16)
    plt.ylabel('Rate (kHz)', fontsize=16)
    plt.grid()
    plt.show()

    #plt.savefig('test.pdf')


    #print(g2)
    # save figures
    fig_save = file_save_plt()
    if fig_save is None:
        None
    else:
        with PdfPages(fig_save) as pdf:
            pdf.savefig()

    fit_vals.to_csv(file_save(), sep='\t', index=False)
    error_vals.to_csv(file_save_err(), sep='\t', index=False)
    # Saving Stark Data\
    P = np.float64(fit_vals['Peak Position'])
    E = np.float64(error_vals['Peak Position err'])
    stark_data = '%s\t%.6f\t%.6f' % (volt, P, E)
    file_name = "Stark_%s_%s_%s.txt" % (sel2, sel, sel3)
    fn = open(file_name, 'a+')
    fn.write(stark_data + '\n')
    fn.close()


exit()
