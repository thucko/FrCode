'''
Version 5
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

warnings.simplefilter(action='ignore', category=FutureWarning)


# Definition for selecting directory
def dir_dialog():
    Tk().withdraw()
    pwd = askdirectory()
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
    numerator = (p_lorentzian[2] / 2)**2
    decay = (np.exp(-x / p_lorentzian[4]))
    denominator1 = ((x - (p_lorentzian[3])) ** 2 + (p_lorentzian[2] / 2) ** 2)
    L = p_lorentzian[0] + decay * p_lorentzian[1] * (numerator / denominator1)
    return L


# Voigt profile
# see https://en.wikipedia.org/wiki/Voigt_profile for info on the Voigt profile
def voigt(x, p):
    sigma = p[5]/(2*np.sqrt(2*np.log(2)))
    gamma = p[2]/2
    z = ((x-p[3]) + 1j*gamma)/(sigma*np.sqrt(2))
    decay = (np.exp(-x / p[4]))
    num = np.real(wofz(z))
    dem = sigma*np.sqrt(2*np.pi)
    V = p[0]+decay*p[1]*(num/dem)
    return V


# Parameters for minimization includes Voigt and Lorentzian parameters
def Parameters(data, x, f):
    # Lorentzain: p = [Offset, Height, FWHM, Center, lifetime]
    # Voigt: p = [Offset, Height, FWHM_L, FWHM_G, Center, Liftime]
    if f == 'Lorentzian':
        p = [np.min(data), np.max(data)-np.min(data), 3.5, np.argmax(np.array(data)), 50]
        return p
    elif f == 'Voigt':

        p = [np.min(data), np.max(data), 3.6, x[np.argmax(np.array(data))], 50, 1]
        return p


# Chi squared function used in the minimization
def chi2(p):

    if sel == 'Lorentzian':
        func = lorentzian(x_step, p)  # Our function we need
    elif sel == 'Voigt':
        func = voigt(x_step, p)
    err = y_err   # Uncertainty in the data sqrt(N)
    # calculate the chi2
    delta = (y_data - func) / err
    chi2_v = sum(pow(delta, 2))
    return chi2_v


# Used for selecting Lorentzian or Voigt fitting fucntions
def ok():
    global sel
    sel = var.get()
    roo.quit()
    roo.destroy()


def select():
    global sel2
    sel2 = var2.get()
    roo2.quit()
    roo2.destroy()




# Grabs the directory for the data
path = dir_dialog()
file_list = glob.glob(path + '/*.txt')
list_file = pd.DataFrame(columns=['Time', 'File path'])

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
button = Button(roo, text="Select", command=ok)
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
button2 = Button(roo2, text="Select", command=select)
button2.pack()
mainloop()


# binwidth in ms
b_width = 2.8


if sel == 'Lorentzian':
    the_fits = pd.DataFrame(columns=['Time', 'Offset', 'Peak Height', 'FWHM', 'Peak Position', 'Lifetime',
                                     'Reduced Chi Squared'])
    the_error = pd.DataFrame(columns=['Time', 'Offset err', 'Peak Height err', 'FWHM err', 'Peak Position err',
                                      'Lifetime err'])

elif sel == 'Voigt':
    the_fits = pd.DataFrame(columns=['Time', 'Offset', 'Peak Height', 'FWHM_L', 'FWHM_G', 'Peak Position', 'Lifetime',
                                     'Reduced Chi Squared'])
    the_error = pd.DataFrame(columns=['Time', 'Offset err', 'Peak Height err', 'FWHM_L err', 'FWHM_G err'
                             , 'Peak Position err', 'Lifetime err'])

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
for file_path in file_list:
    f = open(file_path, 'r')
    cont = f.readlines()
    start_time = cont[6]
    time = start_time[13:28]
    time = time.replace(' ', '')
    # pm_am = time.find('PM')
    if (time.endswith('PM')):

        if time[0:2] == '12':
            time = time.replace('PM', '')
        else:
            time = time.replace('PM', '')
            time = str(0) + time
            hour = int(time[0:2])
            hour = hour + 12
            time = str(hour) + time[2::]
    elif (time.endswith('AM')):
        time = time.replace('AM', '')
        if time[0:2] == '12':
            time = time.replace(time[0:2], '00')
        else:
            if (len(time) != 12):
                time = str(0) + time
    f.close()
    list_file = list_file.append({
        'Time': time,
        'File path': file_path
    }, ignore_index=True)

list_file = list_file.sort_values(by='Time')
prompt1 = sum_scan()
root = Tk()
root.withdraw()
# Ask user for Voltage
volt = simpledialog.askstring(title='Voltage', prompt='Voltage used:')


# begin minimization
print("Starting Minimization...")
# if user selects to look at each scan individually
if prompt1 is False:
    q = 1
    fig_save = file_save_plt()
    pdf = PdfPages(fig_save)
    # Use file name to get data from file, the files should be chronologically ordered at this point
    for index, row in list_file.iterrows():
        times = row['Time']
        file_path = row['File path']
        data_raw = (np.genfromtxt(file_path))
        first_p = 2
        last_p = 8
        pps = last_p - first_p
        data_binned = []
        a = []

        # Gather data points an bin them
        for i in range(0, 160):
            y_bin = (data_raw[i * 8 * 1:(i * 8 * 1) + 8 * 1])
            a = sum(y_bin[first_p:last_p])
            data_binned = np.append(data_binned, a)

        x = np.arange(0, len(data_binned), 1)

        df_Forward = pd.DataFrame({
            'binned PMT data': data_binned[0:80],
            'steps': x[0:80],
            'error': data_binned[0:80] ** 0.5,
        })
        df_Backward = pd.DataFrame({
            'binned PMT data': data_binned[80:159],
            'steps': x[80:159],
            'error': data_binned[80:159] ** 0.5,
        })
        ''' Now to minimize the scans'''

        if sel2 == 'Forward':
            x_step = df_Forward['steps']
            y_data = df_Forward['binned PMT data']/(b_width*pps)
            y_err = df_Forward['error']/(b_width*pps)
            time_stamp = times

        elif sel2 == 'Backward':
            x_step = df_Backward['steps']
            y_data = df_Backward['binned PMT data']/(b_width*pps)
            y_err = df_Backward['error']/(b_width*pps)
            time_stamp = times

        p1 = Parameters(y_data, x_step, sel)  # Generate initial parameter guesses for scans

        # If Lorentzain function is selected for fitting
        if sel == 'Lorentzian':
            #Use for lorentzian fit
            m = Minuit.from_array_func(chi2, p1, error=(10, 1, 0.001, 0.001, .1),
                                       limit=(None, None, (2, 6), None, (0, 500)),
                                       fix=(False, False, False, False, False), errordef=1, pedantic=False)
            m.migrad() # This is minimization strategy
            #       [Offset        , Height        , FWHM          , Center        , lifetime      ]
            p_fit = [m.values["x0"], m.values["x1"], m.values["x2"], m.values["x3"], m.values["x4"]]
            p_err = [m.errors["x0"], m.errors["x1"], m.errors["x2"], m.errors["x3"], m.errors["x4"]]

            Red_chi2 = chi2(p_fit) / (len(y_data) - len(p_fit))
            repack_data = pd.DataFrame({
                'rate': y_data,
                'err': y_err,
                'Normalized Residuals': (y_data - voigt(x_step, p_fit)) / y_err,
                'Residuals': (y_data - voigt(x_step, p_fit)),
                'Freq': x_step,
            })

            the_fits = the_fits.append({
                'Time': time_stamp,
                'Offset': p_fit[0],
                'Peak Height': p_fit[1],
                'FWHM': p_fit[2],
                'Peak Position': p_fit[3],
                'Lifetime': p_fit[4],
                'Reduced Chi Squared': Red_chi2
            }, ignore_index=True)

            the_error = the_error.append({
                'Time': time_stamp,
                'Offset err': p_err[0],
                'Peak Height err': p_err[1],
                'FWHM err': p_err[2],
                'Peak Position err': p_err[3],
                'Lifetime err': p_err[4]
            }, ignore_index=True)
            x_fit = np.arange(min(x_step), max(x_step) + 1, 0.1)
            y_fit = lorentzian(x_fit, p_fit)

        # If Voigt function is selected for fitting
        elif sel == 'Voigt':

            # use for Voigt fit
            m = Minuit.from_array_func(chi2, p1, error=(10, 1, 0.001, 0.001, .1, 0.001),
                                       limit=(None, None, (3, 6), None, (17, 500), (0.6, 4)),
                                       fix=(False, False, False, False, False, False), errordef=1, pedantic=False)
            m.migrad()  # This is minimization strategy
            p_fit = [m.values["x0"], m.values["x1"], m.values["x2"], m.values["x3"], m.values["x4"], m.values["x5"]]
            p_err = [m.errors["x0"], m.errors["x1"], m.errors["x2"], m.errors["x3"], m.errors["x4"], m.errors["x5"]]

            Red_chi2 = chi2(p_fit) / (len(y_data) - len(p_fit))
            repack_data = pd.DataFrame({
                'rate': y_data,
                'err': y_err,
                'Normalized Residuals': (y_data - voigt(x_step, p_fit)) / y_err,
                'Residuals': (y_data - voigt(x_step, p_fit)),
                'Freq': x_step,
            })
            the_fits = the_fits.append({
                'Time': time_stamp,
                'Offset': p_fit[0],
                'Peak Height': p_fit[1],
                'FWHM_L': p_fit[2],
                'FWHM_G': p_fit[5],
                'Peak Position': p_fit[3],
                'Lifetime': p_fit[4],
                'Reduced Chi Squared': Red_chi2
            }, ignore_index=True)

            the_error = the_error.append({
                'Time': time_stamp,
                'Offset err': p_err[0],
                'Peak Height err': p_err[1],
                'FWHM_L err': p_err[2],
                'FWHM_G err': p_err[5],
                'Peak Position err': p_err[3],
                'Lifetime err': p_err[4]
            }, ignore_index=True)
            x_fit = np.arange(min(x_step), max(x_step) + 1, 0.1)
            y_fit = voigt(x_fit, p_fit)

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
        ax1.text(x_fit[0] + 20, np.max(y_data) - 0.17 * np.max(y_data), r'$\nu_{peak}=%.3f$ MHz' % p_fit[3],
                 fontsize=14)
        ax1.text(x_fit[0] + 20, np.max(y_data) - 0.24 * np.max(y_data), r'$\delta\nu_{peak}=%.3f$ MHz' % p_err[3],
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
        pdf.savefig(fig)
        plt.close()
        q = q + 1
    pdf.close()
    print("Minimization complete")

    # save the data and errors to csv files
    the_fits.to_csv(file_save(), sep='\t', index=False)
    the_error.to_csv(file_save_err(), sep='\t', index=False)

# if user selects to sum all the scans together
elif prompt1 is True:
    data = []
    # Use file name to get data from file, the files are chronologically ordered
    for index, row in list_file.iterrows():
        file_path = row['File path']
        data.append(np.genfromtxt(file_path))
    shape = list(data[0].shape)
    shape[:0] = [len(data)]
    data1 = np.concatenate(data).reshape(shape)
    y = data1.sum(axis=0)
    first_p = 2
    last_p = 8
    pps = last_p - first_p
    data_binned = []
    a = []
    # Gather data points an bin them
    for i in range(0, 160):
        y_bin = (y[i * 8 * 1:(i * 8 * 1) + 8 * 1])
        a = sum(y_bin[first_p:last_p])
        data_binned = np.append(data_binned, a)

    x = np.arange(0, len(data_binned), 1)
    df_Forward = pd.DataFrame({
        'binned PMT data': data_binned[0:80],
        'steps': x[0:80],
        'error': data_binned[0:80] ** 0.5,
    })
    df_Backward = pd.DataFrame({
        'binned PMT data': data_binned[80:159],
        'steps': x[80:159],
        'error': data_binned[80:159] ** 0.5,
    })
    ''' Now to minimize the scans'''

    if sel2 == 'Forward':
        x_step = df_Forward['steps']
        y_data = df_Forward['binned PMT data']/(b_width*n_of_scans*pps)
        y_err = df_Forward['error']/(b_width*n_of_scans*pps)
        time_stamp = 'N/A'

    elif sel2 == 'Backward':
        x_step = df_Backward['steps']
        y_data = df_Backward['binned PMT data']/(b_width*n_of_scans*pps)
        y_err = df_Backward['error']/(b_width*n_of_scans*pps)
        time_stamp = 'N/A'


    p1 = Parameters(y_data, x_step ,sel)  # Generate initial parameter guesses for scans

    # If Lorentzain function is selected for fitting
    if sel == 'Lorentzian':
        # Use for lorentzian fit
        m = Minuit.from_array_func(chi2, p1, error=(10, 1, 0.001, 0.001, .1),
                                   limit=(None, None, (2, 6), None, (0, 500)),
                                   fix=(False, False, False, False, False), errordef=1, pedantic=False)
        m.migrad()  # This is minimization strategy
        #       [Offset        , Height        , FWHM          , Center        , lifetime      ]
        p_fit = [m.values["x0"], m.values["x1"], m.values["x2"], m.values["x3"], m.values["x4"]]
        p_err = [m.errors["x0"], m.errors["x1"], m.errors["x2"], m.errors["x3"], m.errors["x4"]]
        Red_chi2 = chi2(p_fit) / (len(y_data) - len(p_fit))
        print(Red_chi2)
        repack_data = pd.DataFrame({
            'rate': y_data,
            'err': y_err,
            'Normalized Residuals': (y_data - voigt(x_step, p_fit)) / y_err,
            'Residuals': (y_data - voigt(x_step, p_fit)),
            'Freq': x_step,
        })

        the_fits = the_fits.append({
            'Time': time_stamp,
            'Offset': p_fit[0],
            'Peak Height': p_fit[1],
            'FWHM': p_fit[2],
            'Peak Position': p_fit[3],
            'Lifetime': p_fit[4],
            'Reduced Chi Squared': Red_chi2
        }, ignore_index=True)

        the_error = the_error.append({
            'Time': time_stamp,
            'Offset err': p_err[0],
            'Peak Height err': p_err[1],
            'FWHM err': p_err[2],
            'Peak Position err': p_err[3],
            'Lifetime err': p_err[4]
        }, ignore_index=True)
        x_fit = np.arange(min(x_step), max(x_step) + 1, 0.1)
        y_fit = lorentzian(x_fit, p_fit)

    # If Voigt function is selected for fitting
    elif sel == 'Voigt':

        m = Minuit.from_array_func(chi2, p1, error=(10, 1, 0.001, 0.001, .1, 0.001),
                                   limit=(None, None, (3, 6), None, (17, 500), (0.6, 6)),
                                   fix=(False, False, False, False, False, False), errordef=1, pedantic=False)
        m.migrad()  # This is minimization strategy
        #       [Offset        , Height        , FWHM_L         , Center        , lifetime ,     FWHM_G        ]
        p_fit = [m.values["x0"], m.values["x1"], m.values["x2"], m.values["x3"], m.values["x4"], m.values["x5"]]
        p_err = [m.errors["x0"], m.errors["x1"], m.errors["x2"], m.errors["x3"], m.errors["x4"], m.errors["x5"]]
        Red_chi2 = chi2(p_fit) / (len(y_data) - len(p_fit))
        print(Red_chi2)
        # repackage the data
        repack_data = pd.DataFrame({
            'rate': y_data,
            'err':  y_err,
            'Normalized Residuals': (y_data - voigt(x_step, p_fit)) / y_err,
            'Residuals': (y_data - voigt(x_step, p_fit)),
            'Freq': x_step,
        })

        the_fits = the_fits.append({
            'Time': time_stamp,
            'Offset': p_fit[0],
            'Peak Height': p_fit[1],
            'FWHM_L': p_fit[2],
            'FWHM_G': p_fit[5],
            'Peak Position': p_fit[3],
            'Lifetime': p_fit[4],
            'Reduced Chi Squared': Red_chi2
        }, ignore_index=True)

        the_error = the_error.append({
            'Time': time_stamp,
            'Offset err': p_err[0],
            'Peak Height err': p_err[1],
            'FWHM_L err': p_err[2],
            'FWHM_G err': p_err[5],
            'Peak Position err': p_err[3],
            'Lifetime err': p_err[4]
        }, ignore_index=True)
        x_fit = np.arange(min(x_step), max(x_step) + 1, 0.1)
        y_fit = voigt(x_fit, p_fit)

    fit = pd.DataFrame({
        'x_fit': x_fit,
        'y_fit': y_fit
    })

    #df = pd.melt(repack_data, id_vars=['Freq'], value_vars=['Normalized Residuals', 'Residuals'], var_name='Type')
    plt.style.use('ggplot')
    gs = gridspec.GridSpec(2, 2)
    fig = plt.figure(figsize=(16, 9))

    ax1 = fig.add_subplot(gs[0, :])
    ax1.text(x_fit[0]+20, np.max(y_data)-0.1*np.max(y_data), r'$\chi^2_{reduced}=%.3f$' % Red_chi2, fontsize=14)
    ax1.text(x_fit[0] + 20, np.max(y_data) - 0.17 * np.max(y_data), r'$\nu_{peak}=%.3f$ MHz' % p_fit[3], fontsize=14)
    ax1.text(x_fit[0] + 20, np.max(y_data) - 0.24 * np.max(y_data), r'$\delta\nu_{peak}=%.3f$ MHz' % p_err[3], fontsize=14)
    ax1.set_title('%s Scan @ %s (%s fit) Summed %d Files' % (sel2, volt, sel, n_of_scans), fontsize=20)
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
    #plt.show()
    #plt.savefig('test.pdf')


    #print(g2)
    # save figures
    fig_save = file_save_plt()
    if fig_save is None:
        None
    else:
        with PdfPages(fig_save) as pdf:
            pdf.savefig()

    the_fits.to_csv(file_save(), sep='\t', index=False)
    the_error.to_csv(file_save_err(), sep='\t', index=False)
    # Saving Stark Data\
    P = np.float64(the_fits['Peak Position'])
    E = np.float64(the_error['Peak Position err'])
    stark_data = '%s\t%.6f\t%.6f' % (volt, P, E)
    file_name = "Stark_%s_%s_%s.txt" % (sel2, sel, sel3)
    fn = open(file_name, 'a+')
    fn.write(stark_data + '\n')
    fn.close()


exit()
