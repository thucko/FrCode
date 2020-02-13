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
d4 = pd.read_csv(files[3], dtype='a')
d5 = pd.read_csv(files[4], dtype='a')


d1['X'] = (np.float64(d1['X']))*(np.float64(d1['Increment'][0])/1e-3)
d2['X'] = (np.float64(d2['X']))*(np.float64(d2['Increment'][0])/1e-3)
d3['X'] = (np.float64(d3['X']))*(np.float64(d3['Increment'][0])/1e-3)
d4['X'] = (np.float64(d4['X']))*(np.float64(d4['Increment'][0])/1e-3)
d5['X'] = (np.float64(d5['X']))*(np.float64(d5['Increment'][0])/1e-3)

d1['CH1'] = (np.float64(d1['CH1']))/1e-3
d2['CH1'] = (np.float64(d2['CH1']))/1e-3
d3['CH1'] = (np.float64(d3['CH1']))/1e-3
d4['CH1'] = (np.float64(d4['CH1']))/1e-3
d5['CH1'] = (np.float64(d5['CH1']))/1e-3


data = pd.DataFrame({
    'No Light': d1['CH1'],
    'Single Pass (rf sw)': d2['CH1'],
    'Single Pass (AM sw)': d3['CH1'],
    'Two AOMs (AM sw)': d4['CH1'],
    'Double Pass (rf sw)': d5['CH1'],
    'Time': d1['X']

})

df = pd.melt(data, id_vars=['Time'], value_vars=['No Light', 'Single Pass (rf sw)', 'Single Pass (AM sw)'
    , 'Two AOMs (AM sw)', 'Double Pass (rf sw)'])


g1 = (ggplot(df, aes(x='Time', y='value', color='variable', shape='variable'))
    + geom_line()
    + geom_point()
    + xlab('Time (ms)')
    + ylab('Volts (mV)')
    #+xlim(4, 10)
    + ggtitle('Switching of AOMs for Trap light Chopping')
    + scale_color_manual(['red', 'black', 'green', 'blue', 'purple'])
)

g2= (ggplot(df, aes(x='Time', y='value', color='variable', shape='variable'))
    + geom_line()
    + geom_point()
    + xlab('Time (ms)')
    + ylab('Volts (mV)')
    + xlim(5, 10)
    + ylim(5, 12)
    + ggtitle('Switching of AOMs for Trap light Chopping (Off portion)')
    + scale_color_manual(['red', 'black', 'green', 'blue', 'purple'])
)
g3 = (ggplot(df, aes(x='Time', y='value', color='variable', shape='variable'))
    + geom_line()
    + geom_point()
    + xlab('Time (ms)')
    + ylab('Volts (mV)')
    + xlim(15, 20)
    + ylim(5, 12)
    + ggtitle('Switching of AOMs for Trap light Chopping (Off portion)')
    + scale_color_manual(['red', 'black', 'green', 'blue', 'purple'])
)
g4 = (ggplot(df, aes(x='Time', y='value', color='variable', shape='variable'))
    + geom_line()
    + geom_point()
    + xlab('Time (ms)')
    + ylab('Volts (mV)')
    + xlim(25, 30)
    + ylim(5, 12)
    + ggtitle('Switching of AOMs for Trap light Chopping (Off portion)')
    + scale_color_manual(['red', 'black', 'green', 'blue', 'purple'])
)
g5 = (ggplot(df, aes(x='Time', y='value', color='variable', shape='variable'))
    + geom_line()
    + geom_point()
    + xlab('Time (ms)')
    + ylab('Volts (mV)')
    + xlim(35, 40)
    + ylim(5, 12)
    + ggtitle('Switching of AOMs for Trap light Chopping (Off portion)')
    + scale_color_manual(['red', 'black', 'green', 'blue', 'purple'])
)

g6 = (ggplot(df, aes(x='Time', y='value', color='variable', shape='variable'))
    + geom_line()
    + geom_point()
    + xlab('Time (ms)')
    + ylab('Volts (mV)')
    + xlim(5, 5.5)
    + ylim(5, 12)
    + ggtitle('Switching of AOMs for Trap light Chopping (zoomed turn off)')
    + scale_color_manual(['red', 'black', 'green', 'blue', 'purple'])
)
g7 = (ggplot(df, aes(x='Time', y='value', color='variable', shape='variable'))
    + geom_line()
    + geom_point()
    + xlab('Time (ms)')
    + ylab('Volts (mV)')
    + xlim(15, 15.5)
    + ylim(5, 12)
    + ggtitle('Switching of AOMs for Trap light Chopping (zoomed turn off)')
    + scale_color_manual(['red', 'black', 'green', 'blue', 'purple'])
)
g8 = (ggplot(df, aes(x='Time', y='value', color='variable', shape='variable'))
    + geom_line()
    + geom_point()
    + xlab('Time (ms)')
    + ylab('Volts (mV)')
    + xlim(25, 25.5)
    + ylim(5, 12)
    + ggtitle('Switching of AOMs for Trap light Chopping (zoomed turn off)')
    + scale_color_manual(['red', 'black', 'green', 'blue', 'purple'])
)
g9 = (ggplot(df, aes(x='Time', y='value', color='variable', shape='variable'))
    + geom_line()
    + geom_point()
    + xlab('Time (ms)')
    + ylab('Volts (mV)')
    + xlim(35., 35.5)
    + ylim(5, 12)
    + ggtitle('Switching of AOMs for Trap light Chopping (zoomed turn off)')
    + scale_color_manual(['red', 'black', 'green', 'blue', 'purple'])
)
print(g1)
save_as_pdf_pages([g1, g2, g3, g4, g5, g6, g7, g8, g9], filename='AOM_switching_4.pdf', path=dir_dialog())

