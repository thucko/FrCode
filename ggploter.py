"""
Author: Tim Hucko

Version 2

"""




import numpy as np
from plotnine import *
import pandas as pd
from tkinter.filedialog import askopenfilename
from tkinter.filedialog import asksaveasfile
from tkinter import simpledialog
from tkinter.filedialog import askdirectory
from tkinter import *
from scipy.stats import norm



def file_save():
    f = asksaveasfile(mode='w', defaultextension=".txt")
    if f is None: # asksaveasfile return `None` if dialog closed with "cancel".
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
            'Peak Position (MHz)': pos,
            'err_m': pos - err,
            'err_mx': pos + err,
            'Time (mins)': x,
            'Reduced Chi Squared': chi}
    )
elif sel2 == 'Backward':
    df = pd.DataFrame({
        'Peak Position (MHz)': 158 - pos,
        'err_m': (158 - pos) - err,
        'err_mx': (158 - pos) + err,
        'Time (mins)': x,
        'Reduced Chi Squared': chi}
    )



mu, std = norm.fit(df['Peak Position (MHz)'])

p = [mu, std]

x = np.arange(np.min(df['Peak Position (MHz)'])-0.15, np.max(df['Peak Position (MHz)'])+0.15, 0.01)

f = norm.pdf(x, p[0], p[1])

g = pd.DataFrame({
    'x': x,
    'G': f
})


l = np.sqrt(len(df['Peak Position (MHz)']))

width = (np.max(df['Peak Position (MHz)'])-np.min(df['Peak Position (MHz)']))/l

width2 = (np.max(df['Reduced Chi Squared'])-np.min(df['Reduced Chi Squared']))/(l*5)

y1 = max(df['Peak Position (MHz)'])+0.1
y2 = min(df['Peak Position (MHz)'])-0.1

g1 = (ggplot(df, aes(x='Time (mins)', y='Peak Position (MHz)'))
     + ggtitle('%s Scans Peak Positions @ %s' % (sel2, volt))
     + geom_point(color='red')
     + geom_errorbar(aes(x='Time (mins)', ymin='err_m', ymax='err_mx'))
     + ylim(y2, y1)

     )
g2 = (ggplot(df, aes(x='Time (mins)', y='Reduced Chi Squared'))
     + ggtitle('%s Scans Reduced Chi Squared @ %s' % (sel2, volt))
     + geom_point(color='green')
            )
g3 = (ggplot(df, aes(x='Peak Position (MHz)', y='stat(density)'))
     + ggtitle('%s Scans Histogram of Peak Positions @ %s \n mu = %.3f, sigma=%.3f' % (sel2, volt, p[0], p[1]))
     + geom_histogram(color='red', binwidth=width)
     + geom_line(g, aes(x='x', y='G'), color='blue')
     + ylab('Density')
              )

g4 = (ggplot(df, aes('Reduced Chi Squared'))
     + ggtitle('%s Scans Histogram of Reduced Chi Squared @ %s' % (sel2, volt))
     + geom_histogram(color='green', binwidth=width2)
      )

#print(g3)
save_as_pdf_pages([g1, g2, g3, g4], filename='%s_scan_peak_positions_%s.pdf' %(sel2, volt), path=dir_dialog())