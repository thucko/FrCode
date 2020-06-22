import numpy as np
from itertools import product
from scipy.integrate import quad
from sympy.physics.quantum.cg import Wigner6j

'''Some Constants'''
h = 6.62607004E-34 # Plancks Constant (J s)
h_bar = h/(2*np.pi)
c = 3e8 # speed of light m/s


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
        list_of_atoms = {'Rb87': {'I':3/2,
                                  'G': {'A': 3417.305, 'B': 0},
                                  'D1': {'A': 408.328, 'B': 0, 'lw': 5.746, 't': 27.70E-9, 'w': 377.107463},
                                  'D2': {'A': 84.7185, 'B': 12.4965, 'lw': 6.065, 't': 26.24E-9, 'w': 384.230484468}
                                  },
                         'Fr210': {'I': 6.0},
                         'Fr211': {'I': 9/2}
                         }
        self.nuclear_spin = list_of_atoms[atom]['I']
        self.hf_A = list_of_atoms[atom][state]['A']
        self.hf_B = list_of_atoms[atom][state]['B']
        if state != 'G':
            self.linewidth = list_of_atoms[atom][state]['lw']
            self.lifetime = list_of_atoms[atom][state]['t']
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
                self.hf_split.update({f: 0.5*A*K})
            else:
                self.hf_split.update({f: 0.5*A*K + B*(1.5*K*(K+1)-2*I*(I+1)*J*(J+1))/(2*I*(2*I-1)*2*J*(2*J-1))})


class Transition:
    def __init__(self, initial_state, final_state):
        self.initial_state = initial_state
        self.final_state = final_state
        self.transition_str = {}
        self.transition_prob = {}
        self.transition_freq = {}

    def get_transition_freq(self, F_init, g_hf, e_hf, freq):
        t_freq = freq - g_hf[F_init]/1E12
        for i in e_hf.keys():
            deltaF = i - F_init
            if deltaF >= -1 and deltaF <= 1:
                t_freq = t_freq + e_hf[i]/1E12
                self.transition_freq.update({str(F_init) + '->' + str(int(i)): freq})



    def get_transition_str(self, Fi, Ff, Ji, Jf, I):
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

    def burst(self, F_init, F_final, hf_vals, linewidth, wavelength, lifetime, optical_power=100E-3):
        print('The desired transition is: F=', F_init, '-> F`=', F_final)
        Is = (np.pi * h * c) / (3 * lifetime * wavelength ** 3)  # saturation intensity
        Is = Is / 1e-4
        I = optical_power/ (np.pi * 2.54 ** 2)  # light intensity W/cm^2
        s1 = I / Is  # Stauration parameter
        allowed_tran =[]
        s_w = {}
        for i in hf_vals.keys():
            deltaF = i - F_init
            if deltaF >= -1 and deltaF <= 1:
                allowed_tran.append(str(F_init) + '->' + str(int(i)))
                laser = np.abs(hf_vals[F_final]-hf_vals[i])
                sw_val = (lorentzian(laser, [linewidth, 0], s1))
                s_w.update({str(F_init) + '->' + str(int(i)): sw_val})





if __name__ == '__main__':
    D2 = AtomicState(atom='Rb87', state='D2')
    D2.get_Qnums(I_qnum=D2.nuclear_spin, J_qnum=3/2)
    D2.get_gfactors(e_qnums=D2.Qnums_e, F_qnums=D2.F_i, I_qnum=D2.nuclear_spin, state='p')
    D2.get_hyperfine_splitting(F_qnums=D2.F_i, J=D2.Qnums_e[2], I=D2.nuclear_spin, A=D2.hf_A, B=D2.hf_B )

    D1 = AtomicState(atom='Rb87', state='D1')
    D1.get_Qnums(I_qnum=D1.nuclear_spin, J_qnum=1/2)
    D1.get_gfactors(e_qnums=D1.Qnums_e, F_qnums=D1.F_i, I_qnum=D1.nuclear_spin, state='p')
    D1.get_hyperfine_splitting(F_qnums=D1.F_i, J=D1.Qnums_e[2], I=D1.nuclear_spin, A=D1.hf_A, B=D1.hf_B)

    G = AtomicState(atom='Rb87', state='G')
    G.get_Qnums(I_qnum=G.nuclear_spin, J_qnum=1/2)
    G.get_gfactors(e_qnums=G.Qnums_e, F_qnums=G.F_i, I_qnum=G.nuclear_spin, state='s')
    G.get_hyperfine_splitting(F_qnums=G.F_i, J=G.Qnums_e[2], I=G.nuclear_spin, A=G.hf_A, B=G.hf_B)

    TS_D2 = Transition(initial_state='G', final_state='D2')
    TS_D2.get_transition_freq(F_init=2, g_hf=G.hf_split, e_hf=D2.hf_split, freq=D2.freq)
    TS_D2.get_transition_str(Fi=G.F_i, Ff=D2.F_i, Ji=G.Qnums_e[2], Jf=D2.Qnums_e[2], I=G.nuclear_spin)
    TS_D2.get_transition_prob(F_init=2, F_final=3, detuning=15, hf_vals=D2.hf_split, lifetime=D2.lifetime,
                              linewidth=D2.linewidth, wavelength=780E-9)
    TS_D2.burst(F_init=2, F_final=3, hf_vals=D2.hf_split, lifetime=D2.lifetime, linewidth=D2.linewidth, wavelength=780E-9)

    TS_D1 = Transition(initial_state='G', final_state='D1')
    TS_D1.get_transition_str(Fi=G.F_i, Ff=D1.F_i, Ji=G.Qnums_e[2], Jf=D1.Qnums_e[2], I=G.nuclear_spin)



    print('Done')