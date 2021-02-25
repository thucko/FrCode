"""
Author: Tim Hucko

Version 3:
- Renamed to peak_plotter
- Using Matlibplot for subplot plotting

"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
from tkinter.filedialog import askopenfilename
from tkinter.filedialog import asksaveasfile
from tkinter import simpledialog
from tkinter.filedialog import askdirectory
from tkinter import *
from scipy.stats import norm
from scipy.stats import iqr
from matplotlib.backends.backend_pdf import PdfPages



# Definition for saving figures
def file_save_plt():
    f = asksaveasfile(mode='wb', defaultextension=".pdf", title='Save Figure(s)',
                      initialfile='%s_scan_%s_Plot' % (sel2, volt))
    if f is None:  # asksaveasfile return `None` if dialog closed with "cancel".
        return
    return f


def file_dialog():
    Tk().withdraw()
    file = askopenfilename()
    return file


def dir_dialog():
    Tk().withdraw()
    pwd = askdirectory()
    return pwd


def parameters(d):
    p = [np.average(d), np.std(d)]
    return p

def gauss(x, p):
    denom = p[1]*np.sqrt(2*np.pi)
    num = np.exp(-0.5*((x-p[0])/p[1])**2)
    g = num/denom
    return g




def select():
    global sel2
    sel2 = var2.get()
    roo2.quit()
    roo2.destroy()

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
root = Tk()
root.withdraw()
volt = simpledialog.askstring(title='Voltage', prompt='Voltage used:')

print('Load peak positions')
pos = np.genfromtxt(file_dialog(), dtype='float64', usecols=5,
                    skip_header=True)
print('Load peak positions uncertainties')
err = np.genfromtxt(file_dialog(), dtype='float64', usecols=5,
                    skip_header=True)
print('Load reduced chi square values')
chi = np.genfromtxt(file_dialog(), dtype='float64', usecols=7,
                    skip_header=True)
x = np.arange(0, len(pos), 1)*20/60

if sel2 == 'Forward':
    df = pd.DataFrame({
            'Peak': pos,
            'err': err,
            'Time': x,
            'Chi': chi}
    )
elif sel2 == 'Backward':
    df = pd.DataFrame({
        'Peak': 158 - pos,
        'err': err,
        'Time': x,
        'Chi': chi}
    )



mu, std = norm.fit(df['Peak'])
bin_num = int((np.max(df['Peak'])-np.min(df['Peak']))/(2*iqr(df['Peak']*len(df['Peak'])**(-0.3333))))
bin_num2 = int((np.max(df['Chi'])-np.min(df['Chi']))/(2*iqr(df['Chi']*len(df['Chi'])**(-0.3333))))
p = [mu, std]

x_fit = np.arange(np.min(df['Peak'])-0.15, np.max(df['Peak'])+0.15, 0.01)
y_fit = norm.pdf(x_fit, p[0], p[1])

y1 = max(df['Peak'])+0.1
y2 = min(df['Peak'])-0.1

plt.style.use('../matplotlib_style/stylelib/cern_root.mplstyle')
gs = gridspec.GridSpec(2, 2)
fig = plt.figure(figsize=(16, 9))
ax1 = fig.add_subplot(gs[0, 0])
ax1.set_title('%s Scans Peak Positions @ %s' % (sel2, volt))
ax1.set_ylabel('Peak Position (MHz)')
ax1.set_xlabel('Time (mins)')
ax1.set_ylim(y2, y1)
ax1.errorbar(df['Time'], df['Peak'], yerr=df['err'], fmt='ro', ecolor='black', capsize=5, zorder=1)

ax2 = fig.add_subplot(gs[1, 0])
ax2.set_title('%s Scans Histogram of Peak Positions @ %s' % (sel2, volt))
ax2.set_ylabel('Density')
ax2.set_xlabel('Peak Position (MHz)')
ax2.text(x_fit[0] + 0.1, np.max(y_fit)*1.1, r'$\mu=%.3f$' % p[0], fontsize=14)
ax2.text(x_fit[0] + 0.1, np.max(y_fit), r'$\sigma=%.3f$' % p[1], fontsize=14)
ax2.hist(df['Peak'], bins=bin_num, density=True, edgecolor='red', facecolor='grey')
ax2.plot(x_fit, y_fit, color='blue')

ax3 = fig.add_subplot(gs[0, 1])
ax3.set_title(r'%s Scans $\chi^2_{Reduced}$ @ %s' % (sel2, volt))
ax3.set_ylabel(r'$\chi^2_{Reduced}$ values')
ax3.set_xlabel('Time (mins)')
ax3.scatter(df['Time'], df['Chi'], color='green')

ax4 = fig.add_subplot(gs[1, 1])
ax4.set_title(r'%s Scans Histogram of $\chi^2_{Reduced}$ @ %s' % (sel2, volt))
ax4.set_ylabel('Density')
ax4.set_xlabel(r'$\chi^2_{Reduced}$ values')
ax4.hist(df['Chi'], bins=bin_num2, density=True, edgecolor='green', facecolor='grey')
fig.subplots_adjust(left=0.06, bottom=0.08, right=0.95, top=0.95, wspace=0.12, hspace=0.25)
#plt.show()

fig_save = file_save_plt()
if fig_save is None:
    None
else:
    with PdfPages(fig_save) as pdf:
        pdf.savefig()

exit()