'''Authour: Tim Hucko
In this code I will try to calculate the atom number. If this works well I will then try to transfer the code to
LabView in an attempt to read out real time atom number instead of ROI
'''
import numpy as np
from photons import AtomicState
from photons import Transition
import matplotlib.pyplot as plt

'''Some constants'''
c = 2.99792458E10  # speed of light in cm/s
h = 6.63E-34  # Planck's constant

'''Laser parameters'''


class LaserParameters:
    def __init__(self):
        self.intensity = None

    def get_laser_intensity(self, measured_power):
        monitor_ratio = 1.6 / 100
        MOT_power = 485 * (measured_power / monitor_ratio - measured_power)  # power in the MOT beams
        beam_ar = np.pi * (2.54 / 2) ** 2  # beam area
        self.intensity = 2 * MOT_power / beam_ar  # intensity in mW/cm^2


class Saturation:
    def __int__(self):
        self.saturation_p = None
        self.scattering_rate = None

    def get_sat_params(self, transition_freq, laser_position, linewidth, lifetime, I_tot):
        detuning = c*((transition_freq*1E12)/c - laser_position)/(linewidth*1E6)
        wavelength = c*1E-2/(transition_freq*1E12)
        Is = (np.pi*h*c)/(3*lifetime*wavelength**3) # saturation intensity
        s = (1000*I_tot/100**2)/Is # saturation parameter
        gamma = linewidth*1E6
        self.saturation_p = s
        sr = (gamma/2)*(s/(1+s+(4*np.pi*detuning/gamma)**2)) # scattering rate calculationg
        self.scattering_rate = sr

class Imaging:
    def __int__(self):
        self.photon_per_count = None
        self.exposure = 0.010
        self.frac_of_photons = None

    def get_photon_count(self, wavelength, ROI):
        photon_E = h*c*1E-2/wavelength
        photon_num = (20E-9)/photon_E
        count_per_photon = ROI/(photon_num*self.exposure)
        self.photon_per_count.extend(1/count_per_photon)

    def lens_system(self, diam, focal_length):
        theta = np.arctan((diam/2)/focal_length)
        Solid_angle = 2*np.pi*(1-np.cos(theta)) # calculate the solid angle of a lens system
        self.frac_of_photons = Solid_angle/(4*np.pi)

class NumberOfAtoms:
    def __int__(self):
        self.atom_number = None

    def get_atom_number(self, ROI, scatter_rate, SAngle_frac, exposure_time, Camera_frac):
        N = 12.5*ROI/(exposure_time*SAngle_frac*scatter_rate*Camera_frac) # calculates the number of atoms
        self.atom_number = N






if __name__ == '__main__':
    Rb = LaserParameters()
    Rb.get_laser_intensity(measured_power=3.356)
    D2 = AtomicState(atom='Rb87', state='D2')
    D2.get_Qnums(I_qnum=D2.nuclear_spin, J_qnum=3 / 2)
    D2.get_gfactors(e_qnums=D2.Qnums_e, F_qnums=D2.F_i, I_qnum=D2.nuclear_spin, state='p')
    D2.get_hyperfine_splitting(F_qnums=D2.F_i, J=D2.Qnums_e[2], I=D2.nuclear_spin, A=D2.hf_A, B=D2.hf_B)
    G = AtomicState(atom='Rb87', state='G')
    G.get_Qnums(I_qnum=G.nuclear_spin, J_qnum=1 / 2)
    G.get_gfactors(e_qnums=G.Qnums_e, F_qnums=G.F_i, I_qnum=G.nuclear_spin, state='s')
    G.get_hyperfine_splitting(F_qnums=G.F_i, J=G.Qnums_e[2], I=G.nuclear_spin, A=G.hf_A, B=G.hf_B)
    TS_D2 = Transition(initial_state='G', final_state='D2')
    TS_D2.get_transition_freq(F_init=2, g_hf=G.hf_split, e_hf=D2.hf_split, freq=D2.freq)
    t_freq = TS_D2.transition_freq['2->3']
    MOT_p = Saturation()
    MOT_p.get_sat_params(transition_freq=t_freq, laser_position=12816.466120, linewidth=D2.linewidth,
                         lifetime=D2.lifetime, I_tot=Rb.intensity)
    img = Imaging()
    img.lens_system(diam=5.08, focal_length=10)
    time = np.genfromtxt('MOTchop2.txt', dtype=None,usecols=0)
    time = time - time[0]
    data = np.genfromtxt('MOTchop2.txt', dtype=None,usecols=2)
    num = NumberOfAtoms()
    num.get_atom_number(ROI=data, scatter_rate=MOT_p.scattering_rate, SAngle_frac=img.frac_of_photons,
                        exposure_time=0.01, Camera_frac=0.3)



    '''Plotting data'''
    plt.style.use('../matplotlib_style/stylelib/cern_root.mplstyle')
    plt.plot(time, num.atom_number)
    plt.ylabel('Number of Atoms (counts)', size=16)
    plt.xlabel('Time (s)', size=16)
    plt.title('MOT Size')
    plt.show()







    print('Done')
