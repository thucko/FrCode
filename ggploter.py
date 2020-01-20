import numpy as np
from plotnine import *
import pandas as pd
from tkinter import Tk
from tkinter.filedialog import askopenfilename
from tkinter.filedialog import asksaveasfile
from tkinter import simpledialog
from tkinter.filedialog import askdirectory



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

root = Tk()
root.withdraw()
volt = simpledialog.askstring(title='Voltage', prompt='Voltage used:')

print('Load peak positions')
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

y1 = max(df['Peak Position (MHz)'])+0.1
y2 = min(df['Peak Position (MHz)'])-0.1

g1 = (ggplot(df, aes(x='Time (mins)', y='Peak Position (MHz)'))
     + ggtitle('Peak Positions @ %s' % volt)
     + geom_point(color='red')
     + geom_errorbar(aes(x='Time (mins)', ymin='err_m', ymax='err_mx'))
     + ylim(y2, y1)

     )
g2 = (ggplot(df, aes(x='Time (mins)', y='Reduced Chi Squared'))
     + ggtitle('Reduced Chi Squared @ %s' % volt)
     + geom_point(color='green')
            )
g3 = (ggplot(df, aes('Peak Position (MHz)'))
     + ggtitle('Histogram of Peak Positions @ %s' % volt)
     + geom_histogram(color='red')
             )
g4 = (ggplot(df, aes('Reduced Chi Squared'))
     + ggtitle('Histogram of Reduced Chi Squared @ %s' % volt)
     + geom_histogram(color='green')
      )

print(g1)
save_as_pdf_pages([g1, g2, g3, g4], filename='peak_pos_%s' % volt, path=dir_dialog())