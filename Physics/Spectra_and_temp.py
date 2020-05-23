'''
Version 2 of data analysis program Inculding Voigt profile
Author: Tim Hucko


'''

import numpy as np
import pandas as pd
from tkinter import Tk
from tkinter import simpledialog
from tkinter.filedialog import askdirectory
from tkinter.filedialog import askopenfilename
from tkinter.filedialog import asksaveasfile
from tkinter import messagebox
import glob
from iminuit import Minuit
from plotnine import *
from astropy.modeling.models import Voigt1D
from scipy.special import wofz
from tkinter import *


import matplotlib.pyplot as plt


def dir_dialog():
    Tk().withdraw()
    pwd = askdirectory()
    return pwd


def file_save():
    f = asksaveasfile(mode='w', defaultextension=".txt")
    if f is None: # asksaveasfile return `None` if dialog closed with "cancel".
        return
    return f


def file_dialog():
    Tk().withdraw()
    file = askopenfilename()
    return file


def file_dialog_temp():
    Tk().withdraw()
    file = askopenfilename()
    x = np.genfromtxt(file, dtype='unicode_', usecols=1)#Times
    y_ULE = np.genfromtxt(file, dtype=None, usecols=3)#ULE temp
    y_box = np.genfromtxt(file, dtype=None, usecols=7)#ULE Box temp
    y_outside = np.genfromtxt(file, dtype=None, usecols=8)#Temo of room near ULE
    return [x, y_ULE, y_box, y_outside]


def save_fig():
    MsgBox = messagebox.askquestion('Saving Spectra', 'Would you like to save the spectra(s)?')
    if MsgBox == 'yes':
        messagebox.showinfo('Saving Spectra', 'Saving the Spectra!')
        return True
    else:
        messagebox.showinfo('Saving Spectra', 'Not saving the Spectra!')
        return False

def sum_scan():
    MsgBox = messagebox.askquestion('Sum Spectra', 'Would you like to sum all the spectra together? (Selecting "No" will'
                                                   ' look at each file indvidually)')
    if MsgBox == 'yes':
        messagebox.showinfo('Sum Spectra', 'Summing the Spectra!')
        return True
    else:
        messagebox.showinfo('Sum Spectra', 'Looking at scans indvidually!')
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

# voigt profile p = [Offset, Height, FWHM, Center, Liftime, Standard deviation]
# see https://en.wikipedia.org/wiki/Voigt_profile for info on the voigt profile
def voigt2(x, p):
    decay = (np.exp(-x/p[4]))
    v1 = Voigt1D(x_0=p[3], amplitude_L=0.18, fwhm_L=p[2], fwhm_G=p[5])
    V = p[0] + decay * p[1] * v1(x)  #voigt function
    return V


def voigt(x, p):
    sigma = p[5]/(2*np.sqrt(2*np.log(2)))
    gamma = p[2]/2
    z = ((x-p[3]) + 1j*gamma)/(sigma*np.sqrt(2))
    decay = (np.exp(-x / p[4]))
    num = wofz(z).real
    dem = sigma*np.sqrt(2*np.pi)
    V = p[0]+decay*p[1]*(num/dem)
    return V


def poly_2 (x,p):
    y = p[0]+p[1]*x+p[2]*x**2
    return y


def Parameters(data, f):
    # p = [Offset, Height, FWHM, Center, lifetime]
    if f == 'Lorentzian':
        p = [np.min(data), np.max(data)-np.min(data), 3.5, np.argmax(np.array(data)), 50]
        return p
    elif f == 'Voigt':
        p = [np.min(data), np.max(data), 3.6, np.argmax(np.array(data)), 50, 1] # use for voigt
        return p



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


def chi2_poly(p):
    func = poly_2(x_poly, p)
    err = the_error['Peak Position err']
    dat = the_fits['Peak Position']
    chi2_p = sum(pow((dat-func) / err, 2))
    return chi2_p


def grab(x):
    start_time = input("Start Time: ")
    end_time = input("End Time: ")
    start_time = str(start_time + ":")  # always take the hh:mm:
    end_time = str(end_time + ":")
    st_index = [i for i, s in enumerate(x[0]) if start_time in s]
    et_index = [i for i, s in enumerate(x[0]) if end_time in s]
    st_index = np.asanyarray(st_index)
    et_index = np.asanyarray(et_index)
    Time = x[0][st_index[0]:et_index[0]]
    ULE_Temp = x[1][st_index[0]:et_index[0]]
    box_Temp = x[2][st_index[0]:et_index[0]]
    outside_temp = x[3][st_index[0]:et_index[0]]
    return[Time, ULE_Temp, box_Temp, outside_temp]


def append_array(x, y):
    x[0] = np.append(x[0], y[0])
    x[1] = np.append(x[1], y[1])
    x[2] = np.append(x[2], y[2])
    x[3] = np.append(x[3], y[3])
    return [x[0], x[1], x[2], x[3]]


def cal(x):
    avg = np.average(x)
    std = np.std(x)
    return [avg, std]

def ok():
    global sel
    sel = var.get()
    roo.quit()


# Grabs the directory for the data
path = dir_dialog()
file_list = glob.glob(path + '/*.txt')
list_file = pd.DataFrame(columns=['Time', 'File path'])
funcs = ["Lorentzian", "Voigt"]
roo = Tk()
roo.geometry('260x75')
roo.resizable(0, 0)
roo.title('Fitting Function')

var = StringVar(roo)
var.set(funcs[0])
w = OptionMenu(roo, var, *funcs)
w.pack()

button = Button(roo, text="Select", command=ok)
button.pack()
mainloop()



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
fig_save = save_fig()
if fig_save is True:
    print('Choose Directory...')
    fig_pwd = dir_dialog()
if fig_save is False:
    print('Not saving spectra...')
# Ask user for Voltage


# begin minimization
print("Starting Minimization...")
q = 1
# get time when file was created, this is important for looking at indvidual scans
for file_path in file_list:
    f = open(file_path, 'r')
    cont = f.readlines()
    start_time = cont[6]
    time = start_time[13:28]
    time = time.replace(' ', '')
    # pm_am = time.find('PM')
    if (time.endswith('PM')):
        time = time.replace('PM', '')
        time = str(0) + time
        hour = int(time[0:2])
        hour = hour + 12
        time = str(hour) + time[2::]
    elif (time.endswith('AM')):
        time = time.replace('AM', '')
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
volt = simpledialog.askstring(title='Voltage', prompt='Voltage used:')

# if user selects to look at each scan indvidually
if prompt1 is False:
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
        df_Backwards = pd.DataFrame({
            'binned PMT data': data_binned[80:159],
            'steps': x[80:159],
            'error': data_binned[80:159] ** 0.5,
        })
        ''' Now to minimize the scans'''
        x_step = df_Forward['steps']
        y_data = df_Forward['binned PMT data']
        y_err = df_Forward['error']
        time_stamp = times

        p1 = Parameters(y_data, sel)  # Generate initial parameter guesses for scans

        # If Lorentzain function is selected for fitting
        if sel == 'Lorentzian':
            #Use for lorentzian fit
            m = Minuit.from_array_func(chi2, p1, error=(100, 1, 0.001, 0.001, .1),
                                       limit=(None, None, (2, 6), None, (0, 500)),
                                       fix=(False, False, False, False, False), errordef=1, pedantic=False)
            m.migrad() # This is minimization strategy
            #       [Offset        , Height        , FWHM          , Center        , lifetime      ]
            p_fit = [m.values["x0"], m.values["x1"], m.values["x2"], m.values["x3"], m.values["x4"]]
            p_err = [m.errors["x0"], m.errors["x1"], m.errors["x2"], m.errors["x3"], m.errors["x4"]]

            Red_chi2 = chi2(p_fit) / (len(y_data) - len(p_fit))
            res = pd.DataFrame({
                '(Observed-Fit)/Error': (y_data - lorentzian(x_step, p_fit)) / df_Forward['error'],
                'Frequency (MHz)': x_step
            })
            # repackage the data
            repack_data = pd.DataFrame({
                'Events': y_data,
                'err_min': y_data - y_err,
                'err_max': y_data + y_err,
                'Frequency (MHz)': x_step,
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
            x_fit = np.arange(min(df_Forward['steps']), max(df_Forward['steps']) + 1, 0.1)
            y_fit = lorentzian(x_fit, p_fit)

        # If Voigt function is selected for fitting
        elif sel == 'Voigt':

            # use for Voigt fit
            m = Minuit.from_array_func(chi2, p1, error=(100, 1, 0.001, 0.001, .1, 0.001),
                                       limit=(None, None, (2, 6), None, (17, 500), (0.5, 4)),
                                       fix=(False, False, False, False, False, False), errordef=1, pedantic=False)
            m.migrad()  # This is minimization strategy
            p_fit = [m.values["x0"], m.values["x1"], m.values["x2"], m.values["x3"], m.values["x4"], m.values["x5"]]
            p_err = [m.errors["x0"], m.errors["x1"], m.errors["x2"], m.errors["x3"], m.errors["x4"], m.errors["x5"]]

            Red_chi2 = chi2(p_fit) / (len(y_data) - len(p_fit))
            res = pd.DataFrame({
                '(Observed-Fit)/Error': (y_data - voigt(x_step, p_fit)) / df_Forward['error'],
                'Frequency (MHz)': x_step
            })
            # repackage the data
            repack_data = pd.DataFrame({
                'Events': y_data,
                'err_min': y_data - y_err,
                'err_max': y_data + y_err,
                'Frequency (MHz)': x_step,
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
            x_fit = np.arange(min(df_Forward['steps']), max(df_Forward['steps']) + 1, 0.1)
            y_fit = voigt(x_fit, p_fit)


        fit = pd.DataFrame({
            'x_fit': x_fit,
            'y_fit': y_fit
        })
        g1 = (ggplot()
              + ggtitle('Forward Scan #%d @ %s (%s fit)' % (q, volt, sel))
              + geom_point(repack_data, aes(x='Frequency (MHz)', y='Events'), color='red')
              + geom_errorbar(repack_data, aes(x='Frequency (MHz)', ymin='err_min', ymax='err_max'))
              + geom_line(fit, aes(x='x_fit', y='y_fit'), color='blue')
              )
        g2 = (ggplot()
              + ggtitle('Residual #%d @ %s (%s fit)' % (q, volt, sel))
              + geom_point(res, aes(x='Frequency (MHz)', y='(Observed-Fit)/Error'), color='red')
             )
        #print(g1, g2)


        # Saving figures
        if fig_save is True:
            ggplot.save(g1, filename='forward_scan_%d_%s_fit_%s.pdf' % (q, sel, volt),
                    path=fig_pwd)
            ggplot.save(g2, filename='residuals_%d_%s_fit_%s.pdf' % (q, sel, volt),
                    path=fig_pwd)
        q = q + 1


    print("Minimization complete")
    theme_set(theme_void())
    # save the data and errors to csv files
    the_fits.to_csv(file_save(), sep='\t', index=False)
    #the_error.to_csv(file_save(), sep='\t', index=False)
    error_bar = pd.DataFrame({
        'x': the_error['Time'],
        'err_min': the_fits['Peak Position'] - the_error['Peak Position err'],
        'err_max': the_fits['Peak Position'] + the_error['Peak Position err']
    })

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
    df_Backwards = pd.DataFrame({
        'binned PMT data': data_binned[80:159],
        'steps': x[80:159],
        'error': data_binned[80:159] ** 0.5,
    })
    ''' Now to minimize the scans'''
    x_step = df_Forward['steps']
    y_data = df_Forward['binned PMT data']
    y_err = df_Forward['error']
    time_stamp = "N/A"

    p1 = Parameters(y_data, sel)  # Generate initial parameter guesses for scans

    # If Lorentzain function is selected for fitting
    if sel == 'Lorentzian':
        # Use for lorentzian fit
        m = Minuit.from_array_func(chi2, p1, error=(100, 1, 0.001, 0.001, .1),
                                   limit=(None, None, (2, 6), None, (0, 500)),
                                   fix=(False, False, False, False, False), errordef=1, pedantic=False)
        m.migrad()  # This is minimization strategy
        #       [Offset        , Height        , FWHM          , Center        , lifetime      ]
        p_fit = [m.values["x0"], m.values["x1"], m.values["x2"], m.values["x3"], m.values["x4"]]
        p_err = [m.errors["x0"], m.errors["x1"], m.errors["x2"], m.errors["x3"], m.errors["x4"]]
        Red_chi2 = chi2(p_fit) / (len(y_data) - len(p_fit))
        print(Red_chi2)
        res = pd.DataFrame({
            '(Observed-Fit)/Error': (y_data - lorentzian(x_step, p_fit)) / df_Forward['error'],
            'Frequency (MHz)': x_step
        })
        # repackage the data
        repack_data = pd.DataFrame({
            'Events': y_data,
            'err_min': y_data - y_err,
            'err_max': y_data + y_err,
            'Frequency (MHz)': x_step,
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
        x_fit = np.arange(min(df_Forward['steps']), max(df_Forward['steps']) + 1, 0.1)
        y_fit = lorentzian(x_fit, p_fit)

    # If Voigt function is selected for fitting
    elif sel == 'Voigt':

        m = Minuit.from_array_func(chi2, p1, error=(100, 1, 0.001, 0.001, .1, 0.001),
                                   limit=(None, None, (2, 6), None, (17, 500), (0.3, 6)),
                                   fix=(False, False, False, False, False, False), errordef=1, pedantic=False)
        m.migrad()  # This is minimization strategy
        #       [Offset        , Height        , FWHM_L         , Center        , lifetime ,     FWHM_G        ]
        p_fit = [m.values["x0"], m.values["x1"], m.values["x2"], m.values["x3"], m.values["x4"], m.values["x5"]]
        p_err = [m.errors["x0"], m.errors["x1"], m.errors["x2"], m.errors["x3"], m.errors["x4"], m.errors["x5"]]
        Red_chi2 = chi2(p_fit) / (len(y_data) - len(p_fit))
        print(Red_chi2)
        res = pd.DataFrame({
            '(Observed-Fit)/Error': (y_data - voigt(x_step, p_fit)) / df_Forward['error'],
            'Frequency (MHz)': x_step
        })
        # repackage the data
        repack_data = pd.DataFrame({
            'Events': y_data,
            'err_min': y_data - y_err,
            'err_max': y_data + y_err,
            'Frequency (MHz)': x_step,
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
        x_fit = np.arange(min(df_Forward['steps']), max(df_Forward['steps']) + 1, 0.1)
        y_fit = voigt(x_fit, p_fit)

    fit = pd.DataFrame({
        'x_fit': x_fit,
        'y_fit': y_fit
    })
    g1 = (ggplot()
          + ggtitle('Forward Scan #%d @ %s (%s fit)' % (q, volt, sel))
          + geom_point(repack_data, aes(x='Frequency (MHz)', y='Events'), color='red')
          + geom_errorbar(repack_data, aes(x='Frequency (MHz)', ymin='err_min', ymax='err_max'))
          + geom_line(fit, aes(x='x_fit', y='y_fit'), color='blue')
          )
    g2 = (ggplot()
          + ggtitle('Residual #%d @ %s (%s fit)' % (q, volt, sel))
          + geom_point(res, aes(x='Frequency (MHz)', y='(Observed-Fit)/Error'), color='red')
          )

    print(g1, g2)
    # save figures
    if fig_save is True:
        ggplot.save(g1, filename='forward_scan_summed_%s_fit_%s.pdf' % (sel, volt),
                    path=fig_pwd)
        ggplot.save(g2, filename='residuals_summed_%s_fit_%s.pdf' % (sel, volt),
                    path=fig_pwd)
    #save minuit output to files
    the_fits.to_csv(file_save(), sep='\t', index=False)
    the_error.to_csv(file_save(), sep='\t', index=False)
# grab the temperature file un-comment if wanting to use
'''file_data = file_dialog_temp()
data = grab(file_data)


user1 = input("Add another file? (y/n): ")
# loop to add more temperature data
while(user1 == "y"):
    file_data2 = file_dialog_temp()
    new_data = grab(file_data2)
    data = append_array(data, new_data)
    user1 = input("Add another file? (y/n): ")'''

# plot peak positions

'''print('Load peak positions')
pos = np.genfromtxt(file_dialog(), dtype='float64', usecols=4,
                    skip_header=True)
print('Load peak positions uncertainties')
err = np.genfromtxt(file_dialog(), dtype='float64', usecols=4,
                    skip_header=True)
print('Load reduced chi square values')
chi = np.genfromtxt(file_dialog(), dtype='float64', usecols=6,
                    skip_header=True)
x = np.arange(0, len(pos), 1)*20/60

df = pd.DataFrame({
        'Peak Position (MHz)': pos,
        'err_m': pos - err,
        'err_mx': pos + err,
        'Time (mins)': x,
        'Reduced Chi Squared': chi}
)

g3 = (ggplot(df, aes(x='Time (mins)', y='Peak Position (MHz)'))
     + ggtitle('Peak Positions')
     + geom_point(color='red')
     + geom_errorbar(aes(x='Time (mins)', ymin='err_m', ymax='err_mx'))
)
g4 = (ggplot(df, aes(x='Time (mins)', y='Reduced Chi Squared'))
     + ggtitle('Reduced Chi Squareds')
     + geom_point(color='green')
            )
g5 = (ggplot(df, aes('Peak Position (MHz)'))
     + ggtitle('Histogram of Peak Positions')
     + geom_histogram(color='red')
             )
g6 = (ggplot(df, aes('Reduced Chi Squared'))
     + ggtitle('Histogram of Reduced Chi Squared')
     + geom_histogram(color='green')
             )
print(g3, g4, g5, g6)
print('Saving figures...')
save_as_pdf_pages([g3, g4, g5, g6], filename=file_save())'''


# if using temperature un-comment this section
#xtic = np.arange(0, len(the_fits))*20/60
"""fig = plt.figure()

#ax1 = fig.add_subplot(111, label="ULE Temperature")
ax2 = fig.add_subplot(111, label="Peak positions", frame_on=False)

#plt.title("ULE Temperatures and Peak Positions", fontsize="20")
plt.title("Peak Positions @ %s" % volt, fontsize="20")



#ax1.get_xaxis().set_visible(False)
'''ax1.plot(data[0], data[1], "bo", label="ULE Temperature", ms=5)
ax1.yaxis.tick_right()
ax1.yaxis.set_label_position('right')
ax1.set_xlabel('Hours', fontsize='15')
ax1.set_ylabel('Temperatures (degree C)', fontsize="15")
'''

#ax2.get_xaxis().set_visible(False)
ax2.set_ylabel('Peak position (MHz)', fontsize="15")
#ax2.xaxis.set_ticks(xtic)
ax2.set_xlabel('Hours', fontsize='15')
ax2.set_ylim([35, 50])
ax2.errorbar(the_fits['Time'], the_fits['Peak Position'], the_error['Peak Position err'],  label="Peak positions",
             fmt='ro')
fig.tight_layout()
plt.legend()
plt.grid()
plt.show()
plt.close('all')"""
exit()
