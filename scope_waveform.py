"""
To read and plot scope waveforms from CSV files

Author: Tim Hucko
Version 0

"""

import numpy as np
from plotnine import *
import pandas as pd
from tkinter.filedialog import askopenfilenames
from tkinter.filedialog import askdirectory
from tkinter import *



def file_dialog():
    Tk().withdraw()
    file = askopenfilenames()
    return file


def dir_dialog():
    Tk().withdraw()
    pwd = askdirectory()
    return pwd




files = file_dialog()

d1 = pd.read_csv(files[0], dtype='a')
d2 = pd.read_csv(files[1], dtype='a')
d3 = pd.read_csv(files[2], dtype='a')


d1['X'] = (np.float64(d1['X']))*(np.float64(d1['Increment'][0])/1e-3)
d2['X'] = (np.float64(d2['X']))*(np.float64(d2['Increment'][0])/1e-3)
d3['X'] = (np.float64(d3['X']))*(np.float64(d3['Increment'][0])/1e-3)


d1['CH1'] = (np.float64(d1['CH1']))/1e-3
d2['CH1'] = (np.float64(d2['CH1']))/1e-3
d3['CH1'] = (np.float64(d3['CH1']))/1e-3


data = pd.DataFrame({
    'One AOM': d1['CH1'],
    'No Light': d2['CH1'],
    'Two AOMs': d3['CH1'],
    'Time': d3['X']

})

df = pd.melt(data, id_vars=['Time'], value_vars=['One AOM', 'No Light', 'Two AOMs'])


g1 = (ggplot(df, aes(x='Time', y='value', color='variable'))
    + geom_line()
    + geom_line()
    + xlab('Time (ms)')
    + ylab('Volts (mV)')
    + xlim(4, 10)
    + ggtitle('Switching of AOMs for Trap light Chopping')
)

g2 = (ggplot(df, aes(x='Time', y='value', color='variable'))
    + geom_line()
    + geom_line()
    + xlab('Time (ms)')
    + ylab('Volts (mV)')
    + xlim(4, 10)
    + ylim(5, 11)
    + ggtitle('Switching of AOMs for Trap light Chopping (zoomed)')
)


save_as_pdf_pages([g1, g2], filename='AOM_switching', path=dir_dialog())

