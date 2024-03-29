import numpy as np
from itertools import product
from sympy.physics.quantum.cg import Wigner6j
from scipy.integrate import quad
from sympy.functions.special.tensor_functions import KroneckerDelta
from gilbertC import GilbertCoefficients
from sympy import Float
import pandas as pd
'''Some Constants'''
h = 6.62607004E-34 # Plancks Constant (J s)
h_bar = h/(2*np.pi)
eps0 = 8.8541878128e-12
c = 2.99792458e8 # speed of light m/s
J_to_ev = 6.242E18
e_charge = 1.602176634e-19
e_mass = 9.1093837015e-31
FSC = 7.29735308e-3
bohr_mag = e_charge*h_bar/(2*e_mass)
M1_AU2SI = 8.478352675883512e-30    # M1_SI = m1_AU2SI * M1_AU
# M1 ampl in SI = M1_redME2SIamp * reduced matrix element
M1_redME2SIamp = (M1_AU2SI/np.sqrt(6))*FSC/2

def lorentzian(x, p, s0):

    N = ((2/p[0])**2)*np.sqrt(1+s0)/(s0*np.pi)
    numerator = (s0/(s0+1))*(p[0] / 2)
    denominator = (1 + (2*(x - p[1])/(p[0]*np.sqrt(1+s0)))**2)
    L = N*(numerator / denominator)
    return L




class AtomicState:

    def __init__(self, atom, state):
        ''' List of atomic properties including: [Nuclear spin, hyperfine constants A & B (G, D1, D2), Linewidth
        and lifetime]'''
        list_of_atoms = {'Rb87': {'I': 3/2,
                                  'G': {'A': 3417.305, 'B': 0.0},
                                  'D1': {'A': 408.328, 'B': 0.0, 't': 27.70E-9, 'w': 377.107463},
                                  'D2': {'A': 84.7185, 'B': 12.4965, 't': 26.24E-9, 'w': 384.230484468}, # w is in THz
                                  'E': {'A': 807.66, 'B': 0.0, 't': 45.6E-9, 'w': 20132.510*c*1e2/1E12}
                                  },
                         'Fr210': {'I': 6.0,
                                   'G': {'A': 7195.1, 'B': 0.0},
                                   'D1': {'A': 408.328, 'B': 0.0, 't': 27.70E-9, 'w': 377.107463},
                                   'D2': {'A': 78.0, 'B': 51.0, 't': 26.24E-9, 'w': (c/718E-9)/1E12},  # w is in THz
                                   'E': {'A': 1577.8, 'B': 0.0, 't': 53.30E-9, 'w': (19732.523*3E10)/1E12}},
                         'Fr211': {'I': 9/2,
                                   'G': {'A': 8713.9, 'B': 0},
                                   'D1': {'A': 1142.1, 'B': 0.0, 't': 29.450E-9, 'w': (c/817E-9)/1E12},
                                   'D2': {'A': 94.9, 'B': -51.0, 't': 21.02E-9, 'w': (c/717.97E-9)/1E12},  # w is in THz
                                   'E': {'A': 1910.9, 'B': 0, 't': 53.30E-9, 'w': (19732.47999*c*1e2)/1E12}},
                         'Cs137': {'I': 7/2,
                                   'G': {'A': 2298.157, 'B': 0.0},
                                   'D1': {'A': 291.920, 'B': 0.0, 't': 34.894E-9, 'w':335.116048807},
                                   'D2': {'A': 50.275, 'B': -0.53, 't': 30.473E-9, 'w':351.72571850},  # w is in THz
                                   'E': {'A': 545.90, 'B': 0.0, 't': 48.28e-9, 'w': (18535.5286*3E10)/1E12}

                         }}
        self.nuclear_spin = list_of_atoms[atom]['I']
        self.hf_A = list_of_atoms[atom][state]['A']
        self.hf_B = list_of_atoms[atom][state]['B']
        if state != 'G':
            self.lifetime = list_of_atoms[atom][state]['t']
            self.linewidth = 1/(self.lifetime*2*np.pi)
            self.freq = list_of_atoms[atom][state]['w']
        self.F_i = []
        self.Qnums_e = []
        self.gf = []
        self.hf_split = {}


    def get_Qnums(self, I_qnum, J_qnum):
        # Let us define some quantum numbers
        S = 1/2 # spin
        J = J_qnum # total angular momentum (Electron)
        L = J-S # orbital angular momentum
        I = I_qnum # Nuclear spin
        F = J+I # total angular momentum (Atom)
        self.Qnums_e.extend([S, L, J]) # list of quantum numbers for the electrons
        self.F_i = np.arange(np.abs(J-I), J+I+1, 1)


    def get_gfactors(self, e_qnums, F_qnums, I_qnum, state):
        atom_states = {'s': 0, 'p': 1, 'd': 2, 'f': 3} # Spectroscopy notation
        L = atom_states[state]
        gj = (1 + (e_qnums[2]*(e_qnums[2]+1) + e_qnums[0]*(e_qnums[0]+1)
                 - L*(L+1))/(2*e_qnums[2]*(e_qnums[2]+1)))
        for f in F_qnums:
            if f == 0:
                self.gf.append('N/A')
            else:
                self.gf.append(gj*(f*(f+1) + e_qnums[2]*(e_qnums[2]+1)-I_qnum*(I_qnum+1))/(2*f*(f+1)))

    def get_hyperfine_splitting(self, F_qnums, J, I, A, B):
        for f in F_qnums:
            K = f*(f+1) - I*(I+1) - J*(J+1)
            if B == 0:
                self.hf_split.update({f: 0.5 * A * K})
            else:
                self.hf_split.update({f: 0.5*A*K + B*(1.5*K*(K+1)-2*I*(I+1)*J*(J+1))/(2*I*(2*I-1)*2*J*(2*J-1))})




class Transition:
    def __init__(self, initial_state, final_state):
        self.initial_state = initial_state
        self.final_state = final_state
        self.transition_str = {}
        self.transition_prob = {}
        self.transition_freq = {}
        self.oscillator_strength = []
        self.transition_rate = pd.DataFrame(columns=['F', 'F`', 'm', 'm`','Rate', 'AM1_hf(AU)','AM1_rel(AU)', 'Ratio(rel/hf)', 'Gilbert Coeff. squared'])
        self.oscillator_strength = []


    def get_transition_freq(self, F_init, g_hf, e_hf, freq):
        s_freq = freq - g_hf[F_init]/1E6
        for i in e_hf.keys():
            deltaF = i - F_init
            if deltaF >= -1 and deltaF <= 1:
                t_freq = s_freq + e_hf[i]/1E6
                self.transition_freq.update({str(F_init) + '->' + str(float(i)): t_freq})



    def get_transition_str_factor(self, Fi, Ff, Ji, Jf, I):
        Fi_to_Ff = list(product(Fi, Ff))
        i = 0
        while i in range(0, len(Fi_to_Ff)):
            delta_F = Fi_to_Ff[i][0]-Fi_to_Ff[i][1]
            if delta_F > 1 or delta_F < -1:
                Fi_to_Ff.remove(Fi_to_Ff[i])
                i = i
            else:
                i = i + 1
        s1 = self.initial_state+'->'+self.final_state
        s2 = self.final_state + '->' + self.initial_state
        SiSf ={s1: {}}
        SfSi = {s2: {}}
        for j in range(0, len(Fi_to_Ff)):
            w6j = Wigner6j(Ji, Jf, 1, Fi_to_Ff[j][1], Fi_to_Ff[j][0], I)
            w6j = w6j.doit()
            SiSf[s1].update({str(int(Fi_to_Ff[j][0]))+'->'+str(int(Fi_to_Ff[j][1])):
                                 (2*Fi_to_Ff[j][1]+1)*(2*Ji+1)*np.power(w6j, 2)})
            SfSi[s2].update({str(int(Fi_to_Ff[j][1]))+'->'+str(int(Fi_to_Ff[j][0])):
                                 (2*Fi_to_Ff[j][0]+1)*(2*Jf+1)*np.power(w6j, 2)})
        self.transition_str.update(SiSf)
        self.transition_str.update(SfSi) 

    def get_transition_prob(self, F_init, F_final, detuning, hf_vals, lifetime, linewidth, wavelength, optical_power=100E-3):
        print('The desired transition is: F=', F_init, '-> F`=', F_final)
        Is = (np.pi*h*c)/(3*lifetime*wavelength**3)  # saturation intensity
        Is = Is/1e-4
        I = optical_power/(np.pi*2.54**2)  # light intensity W/cm^2
        s1 = I / Is # Stauration parameter
        for i in hf_vals.keys():
            deltaF = i - F_init
            if deltaF >= -1 and deltaF <= 1 and i !=F_final:
                laser = np.abs(hf_vals[F_final]-hf_vals[i])-detuning
                section = np.arange(laser, laser + 0.105, 0.01)
                area = quad(lorentzian, section[0], section[len(section) - 1], args=([linewidth, 0], s1))
                prob = area[0]*100
                self.transition_prob.update({str(F_init)+'->'+ str(int(i)): prob})
                print('Transition Probability for F=', F_init, '-> F`=', int(i), ' is ', prob, '%')


    def transition_rate_forb(self, atom, hf_split, tran_freq, lifetime, k_vec, E_laser, P, beam_radius):
        # define M1_rel values from Sauvkov
        M1rel_Dict = {'Rb': 0.17, 'Cs': 10.4, 'Fr': 139.9}
        M1_rel = M1rel_Dict[atom[:2]]*M1_redME2SIamp*1e-5
        ni = AtomicState(atom=atom, state=self.initial_state)
        nf = AtomicState(atom=atom, state=self.final_state)
        ni.get_Qnums(I_qnum=ni.nuclear_spin, J_qnum=1/2)
        nf.get_Qnums(I_qnum=nf.nuclear_spin, J_qnum=1/2)
        En_hf = (abs(hf_split[0][ni.F_i[0]]) +abs(hf_split[0][ni.F_i[1]]))*1E6*h*J_to_ev
        Enprime_hf = (abs(hf_split[1][ni.F_i[0]]) +abs(hf_split[1][ni.F_i[1]]))*1E6*h*J_to_ev
        En_Enprime = tran_freq*1E12*h*J_to_ev
        M1hf =(bohr_mag/c)*pow(En_hf*Enprime_hf, 0.5)/En_Enprime
        #M1hfAU =pow(En_hf * Enprime_hf, 0.5) / En_Enprime

        I = P/(np.pi*beam_radius**2)

        R0 = 2*lifetime*I/(c*eps0*h_bar**2)

        # for hf M1 we only care about F-F' =/= 0
        m_sub = []
        Fi_to_Ff = list(product(ni.F_i, nf.F_i))
        hf_only = True #excludes deltaF = 0 when True
        i = 0
        if hf_only:
            while i in range(0, len(Fi_to_Ff)):
                delta_F = Fi_to_Ff[i][0] - Fi_to_Ff[i][1]
                if delta_F == 0:
                    Fi_to_Ff.remove(Fi_to_Ff[i])
                    i = i
                else:
                    mi = list(np.arange(-Fi_to_Ff[i][0], Fi_to_Ff[i][0]+1, 1))
                    mf = list(np.arange(-Fi_to_Ff[i][1], Fi_to_Ff[i][1] + 1, 1))
                    mlv = list([mi, mf, Fi_to_Ff[i]])
                    m_sub.append(mlv)
                    i = i + 1
        else:
            while i in range(0, len(Fi_to_Ff)):
                mi = list(np.arange(-Fi_to_Ff[i][0], Fi_to_Ff[i][0] + 1, 1))
                mf = list(np.arange(-Fi_to_Ff[i][1], Fi_to_Ff[i][1] + 1, 1))
                mlv = list([mi, mf, Fi_to_Ff[i]])
                m_sub.append(mlv)
                i = i + 1

        q = (-1, 1, 0)
        # define electric field of laser and propagation direction (x, y, z)
        k_x_e = np.cross(k_vec, E_laser)
        state_qnum = []
        for i in range(0, len(m_sub)):
            mi_to_mf = list(product(m_sub[i][0], m_sub[i][1]))
            state_qnum.append([mi_to_mf, m_sub[i][2]])

        j = 0
        amp = []
        gc = GilbertCoefficients()
        while j in range(0, len(state_qnum)):
            Fi = state_qnum[j][1][0]
            Ff = state_qnum[j][1][1]
            for k in state_qnum[j][0]:
                if KroneckerDelta(k[0], k[1]):
                    gc.coefficient(q=0, F=Fi, F_prime=Ff, m=k[0], m_prime=k[1], I=ni.nuclear_spin)
                    gcval = gc.Cval  # This is will take care of <Fm|sigma|F`m`>
                    A_xhf = 0  # no xy-comp for this condition
                    A_yhf = 0
                    A_zhf = M1hf * k_x_e[2] * gcval  # z-comp of the amplitude
                    A_xrel = 0
                    A_yrel = 0
                    A_zrel = M1_rel * k_x_e[2] * gcval  # 07/10/20: relativistic amp z comp
                elif KroneckerDelta(k[0], k[1] - 1): #This is for delta m = -1
                    gc.coefficient(q=+1, F=Fi, F_prime=Ff, m=k[0], m_prime=k[1], I=ni.nuclear_spin)
                    gcval = gc.Cval
                    A_xhf = M1hf * k_x_e[0] * gcval   # xy-comp of the amplitude for q = 1
                    A_yhf = M1hf * (1j * k_x_e[1]) * gcval
                    A_zhf = 0  # no z-comp for this condition
                    A_xrel = M1_rel * k_x_e[0] * gcval  # 07/10/20: relativistic amp xy comp
                    A_yrel = M1_rel * (1j * k_x_e[1]) * gcval
                    A_zrel = 0
                elif KroneckerDelta(k[0], k[1] + 1):
                    gc.coefficient(q=-1, F=Fi, F_prime=Ff, m=k[0], m_prime=k[1], I=ni.nuclear_spin)
                    gcval = gc.Cval
                    A_xhf = M1hf * k_x_e[0] * gcval # x-comp for q = -1
                    A_yhf = M1hf * (-1j * k_x_e[1]) * gcval
                    A_zhf = 0
                    A_xrel = M1_rel * (k_x_e[0]) * gcval
                    A_yrel = M1_rel * (-1j * k_x_e[1]) * gcval
                    A_zrel = 0
                else:
                    A_xhf = Float(0)
                    A_yhf = Float(0)
                    A_zhf = Float(0)
                    A_xrel = Float(0)
                    A_yrel = Float(0)
                    A_zrel = Float(0)

                Ahf = A_zhf + A_xhf + A_yhf  # hf-induced transition amplitude
                Arel = A_zrel + A_xrel + A_yrel  # 07/10/20: added relativistic transition amp
                if KroneckerDelta(Fi, Ff + 1):
                    Atotal = Arel + Ahf
                elif KroneckerDelta(Fi, Ff - 1):
                    Atotal = Arel - Ahf
                elif KroneckerDelta(Fi, Ff):
                    Atotal = Arel
                R = R0 * Atotal * np.conj(Atotal)  # calculate the Rate # 07/10/20: currently have it only for hf.
                R = R.evalf()
                Ahf_au = (c / bohr_mag) * Ahf.evalf()
                #Ahf_au = (Ahf_au*np.conj(Ahf_au))**0.5
                Arel_au = (c / bohr_mag) * Arel.evalf()
                #Arel_au = (Arel_au * np.conj(Arel_au)) ** 0.5
                Aratio = abs(Arel.evalf() / Ahf.evalf())
                GCsquared = float(gcval)
                if R != 0:
                    R = "{:,.3e}".format(R)
                    Ahf_au = "{:,.3e}".format(Ahf_au)
                    Arel_au = "{:,.3e}".format(Arel_au)
                    Aratio = float("{:,.3f}".format(Aratio))
                if KroneckerDelta(k[0], k[1]) or KroneckerDelta(k[0], k[1] + 1) or KroneckerDelta(k[0], k[1]-1):
                    self.transition_rate = self.transition_rate.append({
                        'F': Fi,
                        'F`': Ff,
                        'm': k[0],
                        'm`': k[1],
                        'Rate': R,
                        'AM1_hf(AU)': Ahf_au,
                        'AM1_rel(AU)': Arel_au,
                        'Ratio(rel/hf)': Aratio,
                        'Gilbert Coeff. squared': GCsquared
                    }, ignore_index=True)
            j = j + 1
    def get_oscillator_str(self, rates, wavelength):

        tr_amp = np.array(np.float64(rates['AM1_hf(AU)']) + np.float64(rates['AM1_rel(AU)']))*M1_AU2SI*1e-2
        summed_amp = sum(tr_amp)**2
        os_srt = (4*np.pi*c*e_mass)/(3*h_bar*e_charge**2*wavelength)*tr_amp**2
        self.oscillator_strength.append(os_srt)
        print('done')

    print('done') 
    

       