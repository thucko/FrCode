import numpy as np
from itertools import product
from scipy.integrate import quad
from sympy.physics.quantum.cg import Wigner6j, CG
from sympy.functions.special.tensor_functions import KroneckerDelta
import warnings
import pandas as pd
from atomic_info import AtomicState, Transition 
warnings.simplefilter(action='ignore', category=FutureWarning)





if __name__ == '__main__':
    

    # The atom to be used:
    atom_use = 'Fr211'

    # Get ground state properties
    G = AtomicState(atom=atom_use, state='G')
    G.get_Qnums(I_qnum=G.nuclear_spin, J_qnum=1 / 2)
    G.get_hyperfine_splitting(F_qnums=G.F_i, J=G.Qnums_e[2], I=G.nuclear_spin, A=G.hf_A, B=G.hf_B)
    # Get excited state properties
    E = AtomicState(atom=atom_use, state='E')
    E.get_Qnums(I_qnum=E.nuclear_spin, J_qnum=1 / 2)
    E.get_hyperfine_splitting(F_qnums=E.F_i, J=E.Qnums_e[2], I=E.nuclear_spin, A=E.hf_A, B=E.hf_B)
    # Calculate the transition rate
    TR_ss = Transition(initial_state='G', final_state='E')
    TR_ss.get_transition_freq(F_init=5, g_hf=G.hf_split, e_hf=E.hf_split, freq=E.freq)

    # lets define the power (W) of the laser and beam radius (m)
    power = 30e-3*3980*0.74*0.5
    r = 207e-6
    # power = 30e-3
    # r=1e-3
    TR_ss.transition_rate_forb(atom=atom_use, hf_split=[G.hf_split, E.hf_split], tran_freq=E.freq,
                               lifetime=E.lifetime, P=power, beam_radius=r, k_vec=(0, 1, 0), E_laser=(0, 0, 1))
    TR_ss.get_oscillator_str(rates=TR_ss.transition_rate, wavelength=506e-9)
    TR_ss.transition_rate.to_csv('rate_table.csv')
    # pd.set_option("max_rows", None, "max_columns", None, "expand_frame_repr", False)
    print('Rate calculations for %s' % atom_use)
    print(TR_ss.transition_rate.to_string(index=False, justify='center'))
    print('Average Rate: %.3e per second per atom' %np.mean(np.float64(TR_ss.transition_rate['Rate'][0:27])))
    print('Average Rate: %.3e per second per atom' %np.mean(np.float64(TR_ss.transition_rate['Rate'][27:54])))

    


    print('Done')
