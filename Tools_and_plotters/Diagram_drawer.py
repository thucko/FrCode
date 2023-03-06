import turtle as tu
import numpy as np
# some_file.py
import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, '/home/thucko/Documents/FrCode/Physics')
import photons  as ph
import pandas as pd
from itertools import product
import tkinter
from sympy.functions.special.tensor_functions import KroneckerDelta


# The atom to be used:
atom_use = 'Fr211'

# Get ground state properties
G = ph.AtomicState(atom=atom_use, state='G')
G.get_Qnums(I_qnum=G.nuclear_spin, J_qnum=1 / 2)
G.get_hyperfine_splitting(F_qnums=G.F_i, J=G.Qnums_e[2], I=G.nuclear_spin, A=G.hf_A, B=G.hf_B)
# Get excited state properties
E = ph.AtomicState(atom=atom_use, state='E')
E.get_Qnums(I_qnum=E.nuclear_spin, J_qnum=1 / 2)
E.get_hyperfine_splitting(F_qnums=E.F_i, J=E.Qnums_e[2], I=E.nuclear_spin, A=E.hf_A, B=E.hf_B)
# Calculate the transition rate
TR_ss = ph.Transition(initial_state='G', final_state='E')

# lets define the power (W) of the laser and beam radius (m)
#power = 50e-3*3882*0.74*0.8
#r = 207e-6
power = 50e-3
r=1e-3
TR_ss.transition_rate_forb(atom=atom_use, hf_split=[G.hf_split, E.hf_split], tran_freq=E.freq,
                           lifetime=E.lifetime, P=power, beam_radius=r, k_vec=(0, 1, 0), E_laser=(-(0.5**0.5)*1j, 0, 0.5**0.5))
TR_ss.transition_rate.to_csv('rate_table.csv')
rates = TR_ss.transition_rate
pd.set_option("max_rows", None, "max_columns", None, "expand_frame_repr", False)
print('Rate calculations for %s' % atom_use)
print(TR_ss.transition_rate.to_string(index=False, justify='center'))

ground_m0 = np.arange(-G.F_i[0], G.F_i[0]+1, 1)
ground_m1 = np.arange(-G.F_i[1], G.F_i[1]+1, 1)
excited_m0 = np.arange(-E.F_i[0], E.F_i[0]+1, 1)
excited_m1 = np.arange(-E.F_i[1], E.F_i[1]+1, 1)

polarization = 'Pi'
fnt = 'Times New Roman'
fntsize = 28
tut = tu.setup(width=1920, height=1080)
tut = tu.title(polarization)
tut = tu.Pen()
tut.pen(speed=10)
tut.hideturtle()
x_start = -650
yG0 = G.hf_split[G.F_i[0]]/1000
yG1 = G.hf_split[G.F_i[1]]/1000
x_mv = 50
posiG = []
tut.pensize(10)
tut.penup()
tut.setpos(x_start-100, yG0-400-20)
tut.write('F = %i' %G.F_i[0], align="center", font=(fnt, 34, 'normal'))
for j in ground_m0:
    tut.penup()
    tut.setpos(x_start + x_mv/2, yG0 - 33-400)
    tut.write(j, align="center", font=(fnt, fntsize, 'normal'))
    posiG.append([G.F_i[0], j, tut.pos()])
    tut.setpos(x_start, yG0-400)
    tut.pendown()
    tut.fd(x_mv)
    x_start = x_start + x_mv + 15

x_start = x_start + 100

for j in ground_m1:
    tut.penup()
    tut.setpos(x_start+x_mv/2, yG1 - 33-400)
    tut.write(j, align="center", font=(fnt, fntsize, 'normal'))
    posiG.append([G.F_i[1], j, tut.pos()])
    tut.setpos(x_start, yG1-400)
    tut.pendown()
    tut.fd(x_mv)
    x_start = x_start + x_mv + 15
tut.penup()
tut.setpos(x_start+75, yG1-400-20)
tut.write('F = %i' %G.F_i[1], align="center", font=(fnt, 34, 'normal'))
x_start = -750
posiE = []
yE0 = E.freq+E.hf_split[E.F_i[0]]/1000
yE1 = E.freq+E.hf_split[E.F_i[1]]/1000

tut.penup()
tut.setpos(x_start-100, yE0-250)
tut.write('F` = %i' %E.F_i[0], align="center", font=(fnt, 34, 'normal'))
for j in excited_m0:
    tut.penup()
    tut.setpos(x_start + x_mv/2, yE0-250+5)
    tut.write(j, align="center", font=(fnt, fntsize, 'normal'))
    posiE.append([E.F_i[0], j, tut.pos()])
    tut.setpos(x_start, yE0-250)
    tut.pendown()
    tut.fd(x_mv)
    x_start = x_start + x_mv + 15

x_start = x_start + 100
for j in excited_m1:
    tut.penup()
    tut.setpos(x_start+x_mv/2, yE1-250+5)
    tut.write(j, align="center", font=(fnt, fntsize, 'normal'))
    posiE.append([E.F_i[1], j, tut.pos()])
    tut.setpos(x_start, yE1-250)
    tut.pendown()
    tut.fd(x_mv)
    x_start = x_start + x_mv + 15

tut.penup()
tut.setpos(x_start+75, yE0-250)
tut.write('F` = %i' %E.F_i[1], align="center", font=(fnt, 34, 'normal'))

Fi_to_Ff = list(product(posiG, posiE))
i = 0
while i in range(0, len(Fi_to_Ff)):
    delta_m = Fi_to_Ff[i][0][1]-Fi_to_Ff[i][1][1]
    if delta_m < -1 or delta_m > 1:
        Fi_to_Ff.remove(Fi_to_Ff[i])
        i = i
    else:
        i = i+1
j = 0
vals = []
while j in range(0, len(Fi_to_Ff)):
    check_val = rates.loc[(rates['F'] == Fi_to_Ff[j][0][0]) & (rates['F`'] == Fi_to_Ff[j][1][0])
    & (rates['m'] == Fi_to_Ff[j][0][1]) & (rates['m`'] == Fi_to_Ff[j][1][1])]

    if check_val.empty:
        Fi_to_Ff.remove(Fi_to_Ff[j])
        j = j
    else:
        delta_m = Fi_to_Ff[j][0][1] - Fi_to_Ff[j][1][1]
        if delta_m == 1:
            pol = 'sigma-' # sigma -
        elif delta_m == -1:
            pol = 'sigma+' # sigma +
        elif delta_m == 0:
            pol = 'pi' # pi

        vals.append([Fi_to_Ff[j][0][2], Fi_to_Ff[j][1][2], float(check_val['Rate']), pol])
        j = j+1
k = 0
min_rate = 1
while k in range(0, len(vals)):
    if vals[k][2] == 0:
        vals.remove(vals[k])
        k = k
    else:
        if vals[k][2] < min_rate:
            min_rate = vals[k][2]
        else:
            min_rate = min_rate
        k = k + 1
for i in vals:
    i[2] = (round(i[2]/min_rate))*0.5
    i[0] = tu.Vec2D(float(i[0][0]), float(i[0][1]+38))
    i[1] = tu.Vec2D(float(i[1][0]), float(i[1][1]-10))

tut.penup()
for j in vals:
    if polarization == 'Sigma -':
        if j[3] == 'sigma-':
            tut.pensize(j[2])
            tut.setpos(j[0])
            tut.pendown()
            tut.setpos(j[1])
            tut.penup()
    elif polarization == 'Sigma +':
        if j[3] == 'sigma+':
            tut.pensize(j[2])
            tut.setpos(j[0])
            tut.pendown()
            tut.setpos(j[1])
            tut.penup()
    elif polarization == 'Pi':
        if j[3] == 'pi':
            tut.pensize(j[2])
            tut.setpos(j[0])
            tut.pendown()
            tut.setpos(j[1])
            tut.penup()

tut.getscreen().getcanvas().postscript(file='pi.ps')
tut.screen.mainloop()

print('done')