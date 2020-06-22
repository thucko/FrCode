"""
To read and plot scope waveforms from CSV files

Author: Tim Hucko
Version 0

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
from tkinter.filedialog import askopenfilenames
from tkinter.filedialog import askdirectory
from tkinter import *
import os
from itertools import product



def file_dialog():
    Tk().withdraw()
    file = askopenfilenames()
    return file


def dir_dialog():
    Tk().withdraw()
    pwd = askdirectory()
    return pwd

def powermW(x):
    P = np.float64(10**(x/10))
    return P

class ScopeWaveform:
    def __init__(self):
        self.channels = {}

    def get_data(self, files):

        for x in files:
            file_name = os.path.basename(x)
            d = pd.read_csv(x, dtype='a')
            self.channels.update({file_name: d})


if __name__ == '__main__':
    files = file_dialog()
    data = ScopeWaveform()
    data.get_data(files)
    print('Done')

d1 = pd.read_csv(files[0], dtype='a')
#d2 = pd.read_csv(files[1], dtype='a')
'''d3 = pd.read_csv(files[2], dtype='a')
d4 = pd.read_csv(files[3], dtype='a')
d5 = pd.read_csv(files[4], dtype='a')
d6 = pd.read_csv(files[5], dtype='a')'''



d1['X'] = (np.float64(d1['X']))*(np.float64(d1['Increment'][0]))
#d2['X'] = (np.float64(d2['X']))*(np.float64(d2['Increment'][0]))
'''d3['X'] = (np.float64(d3['X']))*(np.float64(d3['Increment'][0]))
d4['X'] = (np.float64(d4['X']))*(np.float64(d4['Increment'][0]))
d5['X'] = (np.float64(d5['X']))*(np.float64(d5['Increment'][0]))
d6['X'] = (np.float64(d6['X']))*(np.float64(d6['Increment'][0]))'''


d1['CH1'] = (np.float64(d1['CH1']))
#d2['CH1'] = (np.float64(d2['CH1']))
'''d3['CH1'] = (np.float64(d3['CH1']))
d4['CH1'] = (np.float64(d4['CH1']))-9.8
d5['CH1'] = (np.float64(d5['CH1']))
d6['CH1'] = (np.float64(d6['CH1']))'''


'''PmW1 = powermW(d1['CH1'])/1E-6
PmW2 = powermW(d2['CH1'])/1E-6
PmW3 = powermW(d3['CH1'])/1E-6
PmW4 = powermW(d4['CH1'])/1E-6
PmW5 = powermW(d5['CH1'])/1E-6
PmW6 = powermW(d6['CH1'])/1E-6'''

'''d1 = pd.DataFrame({
    'X': d1['X'],
    'CH1': d1['CH1'],
    'PmW': PmW1
})
d2 = pd.DataFrame({
    'X': d2['X'],
    'CH1': d2['CH1'],
    'PmW': PmW2
})
d3 = pd.DataFrame({
    'X': d3['X'],
    'CH1': d3['CH1'],
    'PmW': PmW3
})
d4 = pd.DataFrame({
    'X': d4['X'],
    'CH1': d4['CH1'],
    'PmW': PmW4
})
d5 = pd.DataFrame({
    'X': d5['X'],
    'CH1': d5['CH1'],
    'PmW': PmW5
})

d6 = pd.DataFrame({
    'X': d6['X'],
    'CH1': d6['CH1'],
    'PmW': PmW6
})'''




'''diff1 = pd.DataFrame({
    'X': d1['X'],
    'diff': d1['PmW']-d5['PmW']
})

diff2 = pd.DataFrame({
    'X': d2['X'],
    'diff': d3['PmW']-d6['PmW']
})

diff3 = pd.DataFrame({
    'X': d2['X'],
    'diff': d2['PmW']-d5['PmW']
})
diff4 = pd.DataFrame({
    'X': d2['X'],
    'diff': d4['PmW']-d6['PmW']
})
diff5 = pd.DataFrame({
    'X': d2['X'],
    'diff': d2['PmW']-d1['PmW']
})
diff6 = pd.DataFrame({
    'X': d2['X'],
    'diff': d4['PmW']-d3['PmW']
})

totalp1 = sum(diff1['diff'][0:301])
totalp2 = sum(diff2['diff'][0:301])
totalp3 = sum(diff3['diff'][0:301])
totalp4 = sum(diff4['diff'][0:301])
totalp5 = sum(diff5['diff'][0:301])
totalp6 = sum(diff6['diff'][0:301])

print('Total Power (x-axis)=%.5f' % totalp1)
print('Total Power (y-axis)=%.5f' % totalp2)
'''
plt.style.use('ggplot')
gs = gridspec.GridSpec(2, 2)
fig = plt.figure(figsize=(16, 9))
ax1 = fig.add_subplot(gs[:,:])
ax1.set_title('Cavity Lock FFT')
ax1.set_ylabel(r'Amplitude (dBm)')
ax1.set_xlabel('Frequency (kHz)')
ax1.plot(d1['X'], d1['CH1'], label='Locked')
ax1.set_xlim(0, 10)
'''ax1.plot(d3['X'], d2['CH1'], label='Chamber')'''
#ax1.plot(d2['X'], d2['CH1'], label='Unlocked')
ax1.legend()

'''ax3 = fig.add_subplot(gs[0, 1])
ax3.set_title('Difference')
ax3.set_ylabel(r'Power ($\mu$W)')
ax3.set_xlabel('Frequency (kHz)')
ax3.plot(diff1['X'], diff1['diff'])
ax3.plot(diff3['X'], diff3['diff'])
ax3.plot(diff5['X'], diff5['diff'], color='green')
ax3.fill_between(diff1['X'], diff1['diff'], alpha=0.5,
                 label=r'Total Power = %.3f $\mu$W (Ref - NL, [0, 3 kHz])' % totalp1)
ax3.fill_between(diff3['X'], diff3['diff'], alpha=0.5,
                 label=r'Total Power = %.3f $\mu$W (Chm - NL, [0, 3 kHz])' % totalp3)
ax3.fill_between(diff5['X'], diff5['diff'], alpha=0.5,
                 label=r'Total Power = %.3f $\mu$W (Chm - Ref, [0, 3 kHz])' % totalp5, color='green')
ax3.set_xlim(0, 3)
ax3.legend()

ax2 = fig.add_subplot(gs[1, 0])
ax2.set_title('y-axis FFT')
ax2.set_ylabel(r'Amplitude (dBm)')
ax2.set_xlabel('Frequency (kHz)')
ax2.plot(d2['X'], d3['CH1'], label='Reflection')
ax2.plot(d4['X'], d4['CH1'], label='Chamber')
ax2.plot(d6['X'], d6['CH1'], label='No Light', color='green')
ax2.legend()

ax4 = fig.add_subplot(gs[1, 1])
ax4.set_title('Difference')
ax4.set_ylabel(r'Power ($\mu$W)')
ax4.set_xlabel('Frequency (kHz)')
ax4.plot(diff2['X'], diff2['diff'])
ax4.plot(diff4['X'], diff4['diff'])
ax4.plot(diff6['X'], diff6['diff'], color='green')
ax4.fill_between(diff2['X'], diff2['diff'], alpha=0.5,
                 label=r'Total Power = %.3f $\mu$W (Ref-NL, [0, 3 kHz])' % totalp2)
ax4.fill_between(diff4['X'], diff4['diff'], alpha=0.5,
                 label=r'Total Power = %.3f $\mu$W  (Chm-NL, [0, 3 kHz])' % totalp4)
ax4.fill_between(diff6['X'], diff6['diff'], alpha=0.5,
                 label=r'Total Power = %.3f $\mu$W  (Chm-Ref, [0, 3 kHz])' % totalp6, color='green')
ax4.set_xlim(0, 3)
ax4.legend()
fig.subplots_adjust(left=0.06, bottom=0.08, right=0.95, top=0.95, wspace=0.12, hspace=0.25)'''
plt.show()

