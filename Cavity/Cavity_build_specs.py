'''Cavity_build_specs
Authour: Tim Hucko
Used to calculate optical cavity parameters.
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from sympy.solvers import solve
from sympy import Symbol


calculate_loss = True


def lorentizan (x, x0,gamma):
    L = (1/(2*np.pi))*(gamma/((x-x0)**2+(gamma/2)**2))
    return L


c = 3e10  # speed of light in cm/s
#c = 2.99702547e10
inch_to_cm = 2.54
L = 16  # length of cavity
R1 = 100.0  # Radius of curvature of M1 in cm
R2 = 100.0  # Radius of curvature of M2 in cm
wavelength = 506e-7  # Wavelength of light in cm
k = 0

Ref1 = 0.9991 # Reflection of M1
Ref2 = 0.99995 # Reflection of M2
T1 = 1-Ref1  # Transmission of M1
T2 = 1-Ref2  # Transmission of M2



vFSR = c/(2*L)
dL = (wavelength*1e-2/2)*100e3/vFSR  # Stability of the cavity
g1 = 1-(L/R1)
g2 = 1-(L/R2)
g = g1*g2
zr = np.sqrt(g*(1-g)*L**2/(g1+g2-2*g)**2)  # Rayleigh range
z1 = g2*(1-g1)*L/(g1+g2-2*g)
z2 = g1*(1-g2)*L/(g1+g2-2*g)
w0 = np.sqrt((L*wavelength/np.pi)*np.sqrt(g*(1-g)/(g1+g2-2*g)**2))  # beam waist at the center of the cavity
w1 = np.sqrt((L*wavelength/np.pi)*np.sqrt(g2/(g1*(1-g))))
w2 = np.sqrt((L*wavelength/np.pi)*np.sqrt(g1/(g2*(1-g))))

if calculate_loss is False:
    a = 0 # 4*np.pi*k/wavelength # loss coefficient
    gm = np.sqrt(Ref1 * Ref2 * np.exp(-a * 2 * L))
    F = np.pi * np.sqrt(gm) / (1 - gm)
    t = L * F / (np.pi * c)

elif calculate_loss is True:
    b = 1.11
    t = b*1E-6
    dt = 0.01e-6
    dLen = 0.1
    A = (t*c/L)
    x = Symbol('x')
    gm_cal = solve(x**0.5/(1-x)-A, x)
    gm = np.float64(gm_cal[0])
    a = -(1/(2*L))*np.log(gm**2/(Ref1*Ref2))
    F = np.pi * np.sqrt(gm) / (1 - gm)
    F2 = np.pi*t*c/L
    dF = F*((dt/t)**2+(dLen/L)**2)**0.5
    dgm = gm*(2*np.pi/(F*(4*F**2+np.pi**2)**0.5))*dF
    da = a*((dgm*(2/(gm*np.log(gm**2/(Ref1*Ref2)))))**2 + (dLen/L)**2)**0.5

vFWHM = vFSR/F  # cavity linewidth
''' General Case'''
dv = np.arange(-0.3, 0.3, 0.001)
dPhi = 2*np.pi*dv*1E6*2*L/c
g_v = gm*np.exp(-1j*dPhi)
gamma = (1-gm)**2 + 4*gm*np.sin(dPhi/2)**2  # defines |1-g(v)|^2
''' Gain equations'''
G_g = T1/((1-gm)**2*(1 + (2*F/np.pi)**2*np.sin(dPhi/2)**2))  # cavity gain
R_g = ((Ref1-(Ref1+T1)*gm)**2+4*Ref1*gm*np.sin(dPhi/2)**2)/(Ref1*gamma)  # Reflection dip/gain
T_g = T1*T2*gm/(np.sqrt(Ref1*Ref2)*gamma)  # Transmission gain
'''Phase equations'''
phi_cav = np.arctan(-gm*np.sin(dPhi)/(1-gm*np.cos(dPhi)))
phi_tran = np.arctan((-(1+gm)*np.sin(dPhi/2))/((1-gm)*np.cos(dPhi/2)))
phi_ref = np.arctan(-T1*gm*np.sin(dPhi)/(-Ref1*gamma+T1*gm*(np.cos(dPhi)-1)))


''' For on resonance '''
G = T1/(1-gm)**2
if calculate_loss is True:
    dG = G*(2/(1-gm))*dgm
RefG = (Ref1 - (Ref1+T1)*gm)**2/(Ref1*(1-np.sqrt(Ref1*Ref2))**2)
tranG = T1*T2*gm/(np.sqrt(Ref1*Ref2)*(1-gm)**2)
I = tranG+RefG

if calculate_loss is True:
    print("FSR = %.4E Hz" % vFSR)
    print("Length of Stability = %.4E m" % dL)
    print("zr = %.4f cm" % zr)
    print("z1 = %.4f cm" % z1)
    print("z2 = %.4f cm" % z2)
    print("w0 = %.4f cm" % w0)
    print("w1 = %.4f cm" % w1)
    print("w2 = %.4f cm" % w2)
    print('Cavity Loss = %.2E +/- %.2E cm^-1' % (a, da))
    print("Cavity Storage time = %.2E +/- %.2E s" % (t, dt))
    print("Linewidth = %.4f Hz" % vFWHM)
    print("g1g2 = %.4f" % g)
    print("Finesse = %.0f +/- %.0f" % (F, dF))
    print("Cavity Gain = %.0f +/- %.0f " % (G, dG))
    print("Reflection Power gain = %.4f" % RefG)
    print("Transmission Power gain = %.4f" % tranG)
    print("Reflection + Transmission gain = %.4f" % I)

else:
    print("FSR = %.4E Hz" % vFSR)
    print("Length of Stability = %.4E m" % dL)
    print("zr = %.4f cm" % zr)
    print("z1 = %.4f cm" % z1)
    print("z2 = %.4f cm" % z2)
    print("w0 = %.4f cm" % w0)
    print("w1 = %.4f cm" % w1)
    print("w2 = %.4f cm" % w2)
    print('Cavity Loss = %.2E cm ^-1' % a)
    print("Cavity Storage time = %.2E s" % t)
    print("Linewidth = %.4f Hz" % vFWHM)
    print("g1g2 = %.4f" % g)
    print("Finesse = %.0f" % F)
    print("Cavity Gain = %.0f" % G)
    print("Reflection Power gain = %.4f" % RefG)
    print("Transmission Power gain = %.4f" % tranG)
    print("Reflection + Transmission gain = %.4f" % I)
    #print("Intra-cavity Losses = %.4f" % loss)

'''Plot the general cases for the phase and gains'''
plt.style.use('../matplotlib_style/stylelib/cern_root.mplstyle')
#gs = gridspec.GridSpec(3, 2)
fig = plt.figure()
plt.title('Intra-cavity Gain Curve', size=30)
'''ax1 = fig.add_subplot(gs[0, 0])
ax1.plot(dv, T_g)
ax1.text(-0.2, 1, 'a)', transform=ax1.transAxes, size=14)
ax1.set_title("Gain")
ax1.set_ylabel(r'$G_{trans}$')'''
ax2 = fig.add_subplot()#gs[1, 0])
ax2.plot(dv, G_g)
#ax2.text(-0.2, 1, 'b)',  transform=ax2.transAxes, size=14)
ax2.set_ylabel(r'$G_{cav}$',size=20)
ax2.set_xlabel(r"Detuning $\delta$ (MHz)",size=20)
'''ax3 = fig.add_subplot(gs[2, 0])
ax3.plot(dv, R_g)
ax3.text(-0.2, 1, 'c)', transform=ax3.transAxes, size=14)
ax3.set_ylabel(r'$G_{refl}$')
ax3.set_xlabel(r"Detuning $\delta$ (MHz)")'''
'''ax4 = fig.add_subplot(gs[0, 1])
ax4.set_title('Phase')
ax4.plot(dv, phi_tran)
ax4.set_ylabel(r'$\Phi_{trans}$')
ax5 = fig.add_subplot(gs[1, 1])
ax5.plot(dv, phi_cav)
ax5.set_ylabel(r'$\Phi_{cav}$')
ax6 = fig.add_subplot(gs[2, 1])
ax6.plot(dv, phi_ref)
ax6.set_xlabel(r"Detuning $\delta$ (MHz)")
ax6.set_ylabel(r'$\Phi_{refl}$')'''
plt.xlim(-0.3,0.3)
plt.show()

print('Done')
