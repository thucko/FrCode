import numpy as np
import itertools
from atomic_info import AtomicState
from atomic_info import Transition
from tabulate import tabulate

c = 2.99792458e10 # speed of light m/s


atom = 'Fr211'
if atom == 'Rb87':
    n = 5
    n_prime = 6
    stark_shift_alpha = 0.6056
elif atom == 'Fr211':
    n = 7
    n_prime=8
    stark_shift_alpha = 0.5355


D2 = AtomicState(atom=atom, state='D2')
D2.get_Qnums(I_qnum=D2.nuclear_spin, J_qnum=3 / 2)
D2.get_gfactors(e_qnums=D2.Qnums_e, F_qnums=D2.F_i, I_qnum=D2.nuclear_spin, state='p')
D2.get_hyperfine_splitting(F_qnums=D2.F_i, J=D2.Qnums_e[2], I=D2.nuclear_spin, A=D2.hf_A, B=D2.hf_B)

G = AtomicState(atom=atom, state='G')
G.get_Qnums(I_qnum=G.nuclear_spin, J_qnum=1 / 2)
G.get_gfactors(e_qnums=G.Qnums_e, F_qnums=G.F_i, I_qnum=G.nuclear_spin, state='s')
G.get_hyperfine_splitting(F_qnums=G.F_i, J=G.Qnums_e[2], I=G.nuclear_spin, A=G.hf_A, B=G.hf_B)

E = AtomicState(atom=atom, state='E')
E.get_Qnums(I_qnum=E.nuclear_spin, J_qnum=1 / 2)
E.get_gfactors(e_qnums=E.Qnums_e, F_qnums=E.F_i, I_qnum=E.nuclear_spin, state='s')
E.get_hyperfine_splitting(F_qnums=E.F_i, J=E.Qnums_e[2], I=E.nuclear_spin, A=E.hf_A, B=E.hf_B)



TS_D2 = Transition(initial_state='G', final_state='D2')
G2E = Transition(initial_state='G', final_state='E')

Fg_L = G.F_i[0]
Fg_U = G.F_i[1]

Fe_L = E.F_i[0]
Fe_U = E.F_i[1]


E_field = 0/2.858
freq_shift = (stark_shift_alpha*E_field**2)/2/c # divide by 2 to get freq w.r.t  993 or 1012
print((freq_shift*c*2))
aom_freq = 110E6/2/c
phase_matrix_freq = 60E6
offset = 100E6/c

G2E.get_transition_freq(F_init=Fg_L, g_hf=G.hf_split, e_hf=E.hf_split, freq=E.freq)
L2U = G2E.transition_freq['%.1f->%.1f' %(Fg_L, Fe_U)]/c/2E-12-aom_freq + freq_shift + phase_matrix_freq/c + offset
L2L = G2E.transition_freq['%.1f->%.1f' %(Fg_L, Fe_L)]/c/2E-12-aom_freq + freq_shift + phase_matrix_freq/c + offset

G2E.get_transition_freq(F_init=Fg_U, g_hf=G.hf_split, e_hf=E.hf_split, freq=E.freq)


U2L = G2E.transition_freq['%.1f->%.1f' %(Fg_U, Fe_L)]/c/2E-12-aom_freq + freq_shift + phase_matrix_freq/c + offset
U2U = G2E.transition_freq['%.1f->%.1f' %(Fg_U, Fe_U)]/c/2E-12-aom_freq + freq_shift + phase_matrix_freq/c + offset

table = [['|%is, F`=%i>' %(n_prime, Fe_L), '|%is, F`=%i>' %(n_prime, Fe_U)], 
['|%is, F`=%i>' %(n, Fg_L),'%.5f cm^-1' %L2L, '%.5f cm^-1' %L2U], 
['|%is, F`=%i>' %(n, Fg_U), '%.5f cm^-1' %U2L, '%.5f cm^-1' %U2U]
]

print(tabulate(table, headers='firstrow', tablefmt='double_grid'))