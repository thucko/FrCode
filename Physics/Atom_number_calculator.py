'''Authour: Tim Hucko
In this code I will try to calculate the atom number. If this works well I will then try to transfer the code to
LabView in an attempt to read out real time atom number instead of ROI
'''
import numpy as np
from transition_rates import AtomicState
from transition_rates import Transition
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from iminuit import Minuit
import pandas as pd
import glob
from scipy.interpolate import interp1d

'''Some constants'''
c = 2.99792458E10  # speed of light in cm/s
h = 6.63E-34  # Planck's constant

'''Laser parameters'''


class LaserParameters:
    def __init__(self):
        self.intensity = None

    def get_laser_intensity(self,total_power, coupling):
        #beam_split = 0.33
        MOT_power = 485*total_power * coupling # power in the MOT beams
        beam_ar = np.pi * (2.54/2) ** 2  # beam area
        self.intensity = 2 * MOT_power / beam_ar  # intensity in mW/cm^2


class Saturation:
    def __int__(self):
        self.saturation_p = None
        self.scattering_rate = None

    def get_sat_params(self, transition_freq, linewidth, lifetime, I_tot):
        gamma = linewidth
        detuning = 4*linewidth
        wavelength = c*1E-2/(transition_freq*1E12)
        Is_0 = (1000/100**2)*(np.pi*h*c*1E-2)/(3*lifetime*wavelength**3) # saturation intensity
        Is = Is_0*(1+(2*detuning/gamma)**2)
        s = (I_tot)/Is # saturation parameter
        self.saturation_p = s
        sr = (gamma/2)*(s/(1+s+(2*detuning/gamma)**2)) # scattering rate calculation
        self.scattering_rate = sr

class Imaging:
    def __int__(self):
        self.photon_per_count = None
        self.frac_of_photons = None

    def get_photon_count(self, wavelength, ROI):
        photon_E = h*c*1E-2/wavelength
        photon_num = (20E-9)/photon_E
        exposure = 0.010
        count_per_photon = (ROI/.8)/(photon_num*exposure)
        print(1/count_per_photon)
        self.photon_per_count = (1/count_per_photon)

    def lens_system(self, diam, focal_length):
        theta = np.arctan((diam/2)/focal_length)
        sangle = np.pi*(diam/2)**2/(focal_length**2)
        Solid_angle = 2*np.pi*(1-np.cos(theta))# calculate the solid angle of a lens system
        self.frac_of_photons = Solid_angle/(4*np.pi)

class NumberOfAtoms:
    def __int__(self):
        self.atom_number = None

    def get_atom_number(self, ROI, scatter_rate, SAngle_frac, exposure_time, polarizer_num, ppc):
        theta = (310.0 - polarizer_num)*(np.pi/180)
        Camera_frac = np.cos(theta)**2
        N = ppc*ROI/(exposure_time*SAngle_frac*scatter_rate*Camera_frac) # calculates the number of atoms
        self.atom_number = N


def exp_decay(x, p):
    f = p[0]*np.exp(-x/p[1])+p[2]
    return f


def fun(x, a, b, c):
    f = a*np.exp(-x/b)+c
    return f

def chi2(p):
    func = exp_decay(tpt, p)
    err = pt_err  # Uncertainty in the data sqrt(N)
    # calculate the chi2
    delta = (pt - func) / err
    chi2_v = sum(pow(delta, 2))
    return chi2_v



#if __name__ == '__main__':
tds = np.array([])
N0s = np.array([])
te = np.array([])
Ne = np.array([])
for n in range(1,2):
    atom = 'Fr211'
    # for n in range(1,10):
    numa = n
    print(numa)


    Fr = LaserParameters()
    Fr.get_laser_intensity(total_power=246, coupling=0.69)
    inten = Fr.intensity
    D2 = AtomicState(atom=atom, state='D2')
    D2.get_Qnums(I_qnum=D2.nuclear_spin, J_qnum=3 / 2)
    D2.get_gfactors(e_qnums=D2.Qnums_e, F_qnums=D2.F_i, I_qnum=D2.nuclear_spin, state='p')
    D2.get_hyperfine_splitting(F_qnums=D2.F_i, J=D2.Qnums_e[2], I=D2.nuclear_spin, A=D2.hf_A, B=D2.hf_B)
    G = AtomicState(atom=atom, state='G')
    G.get_Qnums(I_qnum=G.nuclear_spin, J_qnum=1 / 2)
    G.get_gfactors(e_qnums=G.Qnums_e, F_qnums=G.F_i, I_qnum=G.nuclear_spin, state='s')
    G.get_hyperfine_splitting(F_qnums=G.F_i, J=G.Qnums_e[2], I=G.nuclear_spin, A=G.hf_A, B=G.hf_B)
    TS_D2 = Transition(initial_state='G', final_state='D2')
    TS_D2.get_transition_freq(F_init=5, g_hf=G.hf_split, e_hf=D2.hf_split, freq=D2.freq)
    t_freq = TS_D2.transition_freq['5->6.0']


    MOT_p = Saturation()
    MOT_p.get_sat_params(transition_freq=t_freq, linewidth=D2.linewidth,
                         lifetime=D2.lifetime, I_tot=inten)
    img = Imaging()
    img.lens_system(diam=5.08, focal_length=10)
    img.get_photon_count(wavelength=718E-9, ROI=4.75E8)

    path = ('~/Documents/Francium/burst_test/08252022/MOT29/')
    file_list = glob.glob(path + '/*.txt')
 
    file_list = sorted(file_list)[1::]
    
    
    
    
    for file in file_list:
        time = np.genfromtxt(file, dtype=float,usecols=0, skip_header=True)
        buffer = np.genfromtxt(file, dtype=float, usecols=1, skip_header=True)


        res = [idx for idx, item in enumerate(buffer) if item in buffer[:idx]]
        data = np.genfromtxt(file, dtype=float, usecols=2, skip_header=True)

        data = np.delete(data, res)
        time = np.delete(time, res)

        data = np.delete(data, 0)
        time = np.delete(time, 0)

        if len(time) > 100:
            time = time-time[0]
            indx = np.where(np.round(time) == 13)[0][0]
            time = time[0:indx]
          
            num = NumberOfAtoms()
            num.get_atom_number(ROI=data, scatter_rate=MOT_p.scattering_rate, SAngle_frac=img.frac_of_photons,
                                exposure_time=0.01, polarizer_num=235.0, ppc=img.photon_per_count)
            y_data = num.atom_number
            y_data = y_data[0:indx]
            err_bar = np.sqrt(y_data)*np.sqrt(800)
            ### FOR FITTING ###

            diffs = np.diff(y_data)
            mean_diffs = np.mean(diffs)
            std_diffs = np.std(diffs)

            pt = np.array([y_data[0]])
            pt_err = np.array([err_bar[0]])
            tpt = np.array([time[0]])
            for j in range(len(y_data)-1):
                if abs(y_data[j]-y_data[j+1])< mean_diffs+std_diffs: 
                    pt = np.append(pt, y_data[j+1])
                    pt_err = np.append(pt_err, err_bar[j+1])
                    tpt = np.append(tpt, time[j+1])



            pinit = [max(pt)-np.mean(pt), 0.4, np.min(pt)]
            m = Minuit(chi2, pinit)
            m.limits = ((0, None), (0, None), (0, None))
            m.errors = (0.1, 0.1, 0.1)
            m.errordef = Minuit.LEAST_SQUARES
            m.migrad()

            
            p_fit = [m.values['x0'], m.values['x1'], m.values['x2']]
            p_err = [m.errors['x0'], m.errors['x1'], m.errors['x2']]

            red = m.fval / (len(pt) - len(p_fit))
          
            ### INFLATE ERROR BARS ###
            pt_err = pt_err*np.sqrt(red)
            m.migrad()
            p_fit = [m.values['x0'], m.values['x1'], m.values['x2']]
            p_err = [m.errors['x0'], m.errors['x1'], m.errors['x2']]
            red = m.fval / (len(pt) - len(p_fit))



            # print(red)
            N0s = np.append(N0s, round(p_fit[0], 0))
            Ne = np.append(Ne, round(p_err[0], 2))
            tds = np.append(tds, round(p_fit[1], 2))
            te = np.append(te, round(p_err[1], 5))

            x_fit = np.arange(tpt[0], tpt[-1], 0.01)
            y_fit = exp_decay(x_fit, p_fit)        

            ### END FITTING ###
            
    

            '''plot data'''
            plt.style.use('../matplotlib_style/stylelib/cern_root.mplstyle')
            plt.errorbar(tpt, pt, yerr=pt_err, fmt='ro', ecolor='black',capsize=5, zorder=1)
            plt.xlabel('time (s)', fontsize=16)
            plt.ylabel('Atom Number', fontsize=16)
            plt.plot(x_fit, y_fit, 'b')
            plt.show()


    print('Done')
print(list(N0s))
print(list(tds))
print(np.average(N0s), (1/len(Ne))*np.sqrt(np.sum(Ne**2)))
print(np.average(tds), (1/len(te))*np.sqrt(np.sum(te**2)))
