import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def inductor_eq (v0, t0,t, tau):
    return v0*(1-np.exp(-(t-t0)/tau))


'''Specify user current, resistance, and inductance'''
user_current = 45 # Put user current here
user_resistance = 31.6e-3 # Put user resistance here
user_inductance = 162e-6 # Put inductance here

tau = user_inductance/user_resistance # calculate time constant
print('Time constant: %.3E s' %tau)
voltage = user_current*user_resistance # DC voltage to run at user current

t_delay = 1e-3 # initial time delay
on_time = 5e-3 # On time for main pulse
off_time = 3e-3 # off time for main pulse
rise_time = 0.5e-3 # Spike pulse time

#generate time points
delay_t = np.arange(0, t_delay, 10e-6)
on_t = np.arange(0, on_time, 10e-6)
off_t = np.arange(0, off_time, 10e-6)
spike_t = np.arange(0, rise_time, 10e-6)

# full sequence of pulses points (this is a little messy!)
t = np.concatenate((delay_t, spike_t+t_delay, on_t+rise_time+t_delay, spike_t+on_time+rise_time+t_delay,
                    off_t+on_time+2*rise_time+t_delay, spike_t+off_time+on_time+2*rise_time+t_delay,
                    on_t+off_time+on_time+3*rise_time+t_delay ))


spike_voltage = voltage/(1-np.exp(-(rise_time)/tau)) # calculate voltage spike level
print('Spike Voltage: %.3f V' %spike_voltage) # print it out to see value

# get list for voltage levels
y_delay = np.zeros(len(delay_t))
y_spike_up = (spike_voltage)*np.ones(len(spike_t))
y_spike_down = -(spike_voltage-voltage)*np.ones(len(spike_t))
y_up = voltage*np.ones(len(on_t))
y_down = np.zeros(len(off_t))

#Make list of entire voltage levels
y = np.concatenate((y_delay, y_spike_up, y_up, y_spike_down, y_down,y_spike_up, y_up))


wave = pd.DataFrame({
    'x':t,
    'y':y
})

#just in case a plot is wanted to see pulse
'''plt.style.use('../matplotlib_style/stylelib/cern_root.mplstyle')
plt.plot(t, y)
plt.show()'''

#save voltage and time values
wave.to_csv('PWL%iA.txt' % user_current, header=None, index=None, sep=' ')