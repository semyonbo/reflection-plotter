# P polatisation plotter
import matplotlib
from numpy import arcsin, sin, sqrt, pi, cos, arctan, arccos
import matplotlib.pyplot as plt
import numpy as np
from numpy.lib import emath
import os

os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin'

matplotlib.rcParams.update({'font.size': 18})
plt.rcParams['text.usetex'] = True



def get_dielectic(lamb, material):
    if material == 'Au':
        mat = 'refractiveindex/esAu.txt'
    else:
        mat = 'refractiveindex/esAG1.txt'
    with open(mat, 'r+') as file:
        lambs = list(map(lambda x: float(x), file.readline().split(',')))
        epsilons = list(map(lambda x: complex(x), file.readline().replace('i', 'j').split(',')))
        epsilon = epsilons[0]
        if len(lambs) == len(epsilons):
            for i in range(0, len(lambs)):
                if lambs[i - 1] <= lamb * (10 ** (-6)) <= lambs[i]:
                    if (lambs[i - 1] + abs((lambs[i] - lambs[i - 1]) / 2)) >= lamb:
                        epsilon = epsilons[i - 1]
                    else:
                        epsilon = epsilons[i]

    return epsilon


def calc_p_polaris_r(theta, lamb, k1x, type, mat):
    k1 = 2 * pi / (lamb * (10 ** (-6)))
    if type == 'angle':
        k1x = k1 * sin(theta)

    k1z = emath.sqrt(k1 ** 2 - k1x ** 2)

    k2x = k1x

    k2 = emath.sqrt(get_dielectic(lamb, mat)) * k1

    k2z = emath.sqrt(k2 ** 2 - k2x ** 2)

    rp = (-1 * k2z + get_dielectic(lamb, mat) * k1z) / (get_dielectic(lamb, mat) * k1z + 1 * k2z)

    return rp

margins = {
    "left"   : 0.10,
    "bottom" : 0.140,
    "right"  : 0.990,
    "top"    : 1
}



# # Lambda linspace - plotting graph |R_p|^2 (lambda)
# Y = np.linspace(0.1879, 1.937, 100)
#
# # Put angle here:
# angle = 60 * 2 * pi / 360
#
# angle_deg = round(angle * 360 / (2 * pi), 0)
#
# fig2 = plt.figure(2)
#
# ax2 = plt.gca()
# fig2.subplots_adjust(**margins)
#
# R_S_Au = []
# R_S_Ag = []
#
# for i in Y:
#     R_S_Au.append(calc_p_polaris_r(angle, i, 0, 'angle', 'Au')*np.conj(calc_p_polaris_r(angle, i, 0, 'angle', 'Au')))
#     R_S_Ag.append(calc_p_polaris_r(angle, i, 0, 'angle', 'Ag')*np.conj(calc_p_polaris_r(angle, i, 0, 'angle', 'Ag')))
#
# plt.plot(Y, R_S_Au, label=r'$|R^{p}|^{2}$ for Au', color='forestgreen', lw='1.5')
# plt.plot(Y, R_S_Ag, label=r'$|R_{p}|^{2}$ for Ag', color='firebrick', lw='1.5')
# ax2.set_xlim(xmin=0)
# # ax2.set_ylim(ymin=0)
# ax2.set_xlabel(r'$\lambda$ ($\mu m$)')
# plt.legend()
# plt.grid()
# plt.savefig(f"Plot_abs_Rp2_(lamb)_{angle_deg}.pdf", format="pdf", bbox_inches="tight")
# plt.show()
#
#
# # Plotting graph |R_p|^2 (angle):
# # Angle linspace
# X = np.linspace(0, pi / 2, 100)
# # Put lambda here (in micrometre)
# lamb_Au = 0.6
# lamb_Ag = 0.35
#
# fig1 = plt.figure(1)
# fig1.subplots_adjust(**margins)
# ax1 = plt.gca()
# plt.plot(X, calc_p_polaris_r(X, lamb_Au, 0, 'angle', 'Au')*np.conj(calc_p_polaris_r(X, lamb_Au, 0, 'angle', 'Au')),
#          label=r'$|R^{p}|^{2}$ for Au:'+fr"$\lambda={lamb_Au} \;\mu m$", color='firebrick', lw='1.5')
# plt.plot(X, calc_p_polaris_r(X, lamb_Ag, 0, 'angle', 'Ag')*np.conj(calc_p_polaris_r(X, lamb_Ag, 0, 'angle', 'Ag')),
#          label=r'$|R^{p}|^{2}$ for Ag:'+fr"$\lambda={lamb_Ag} \;\mu m$", color='forestgreen', lw='1.5')
#
# plt.legend()
# plt.grid()
# ax1.set_xlabel(r'Angle (rad)')
# #ax1.set_xlim(xmin=0)
# # ax1.set_ylim(ymin=0)
# plt.savefig("Plot_abs_Rp_2_(angle).pdf", format="pdf", bbox_inches="tight")
# # plt.xticks(plt.xticks()[0],[r"$" + format(r/np.pi, ".2g")+ r"\pi$" for r in plt.xticks()[0]])
# plt.show()


# Plotting graph Im(R_p)(angle) for lamda Au and Ag:
# fig4 = plt.figure(4)
# fig4.subplots_adjust(**margins)
# ax4 = plt.gca()
# plt.plot(X, np.imag(calc_p_polaris_r(X, lamb_Au, 0, 'angle', 'Au')), label=r'$Im(R^{p})$ for Au:'+fr"$\lambda={lamb_Au} \;\mu m$", color='firebrick', lw='1.5')
# plt.plot(X, np.imag(calc_p_polaris_r(X, lamb_Ag, 0, 'angle', 'Ag')), label=r'$Im(R^{p})$ for Ag:'+fr"$\lambda={lamb_Ag} \;\mu m$", color='forestgreen', lw='1.5')
# plt.legend()
# plt.grid()
# ax4.set_xlabel(r'Angle (rad)')
# ax4.set_xlim(xmin=0)
# plt.savefig("Plot_Im(Rp)(angle).pdf", format="pdf", bbox_inches="tight")
# plt.show()



# Plot graph |R_p^2| from Array of k_x
# kx array starts from: a*k1
a = 0.5

# kx array ends with: b*k1
b = 2

# Amount of dots between a and b:
resolution = 200

# Lambda of wave:
lamb = 0.6

X1 = np.linspace(a, b, resolution)
fig3 = plt.figure(3)
fig3.subplots_adjust(**margins)
ax3 = plt.gca()
angle = 0
k1 = 2 * pi / (lamb * 10 ** (-6))
X = X1 * k1

kx_crit = (2 * pi / (lamb * 10 ** (-6))) * emath.sqrt(1 * get_dielectic(lamb,'Au') / (1 + get_dielectic(lamb,'Au')))

plt.plot(X, calc_p_polaris_r(angle, lamb, X, 'None','Au') * np.conj(calc_p_polaris_r(angle, lamb, X, 'None','Au')),
         label=r'$|R^{p}|^{2}$ for Au', color='firebrick', lw=1.5)
plt.axvline(np.real(kx_crit), color='forestgreen', linestyle='dashed', lw=1.5,
            label=r'Critical $k_{1x}$:' + fr' ${round(np.real(kx_crit / k1), 3)} \cdot k_{1}$ ')
plt.legend()
plt.grid()
ax3.set_xlabel(r'$k_{1x}$')
ax3.set_xlim(xmin=a*k1)
plt.hlines(y=0, xmin=0, xmax=k1, color='g', linestyle='-', linewidth=3)
plt.hlines(y=0, xmin=k1, xmax=5 * k1, color='r', linestyle='-', linewidth=3)
# ax2.set_ylim(ymin=0)

plt.xticks(np.arange(a*k1, b * k1 + 1, k1/2))
ax3.text((a+0.05) * k1, -25, 'propagating', style='italic')
ax3.text(1.7 * k1, -25, 'evanescent', style='italic')
plt.xticks(plt.xticks()[0], [r"$" + format(r / k1, ".2g") + r"k_1$" for r in plt.xticks()[0]])
plt.savefig("Plot_abs_Rp_2_evanicent_Au.pdf", format="pdf", bbox_inches="tight")
plt.show()



# # Lambda of wave - plotting compramation betwen two types of materials:
# lamb = 0.6
# a=0.85
# b=1.45
#
# X1 = np.linspace(a, b, resolution)
# fig5 = plt.figure(5)
# fig5.subplots_adjust(**margins)
# ax5 = plt.gca()
# angle = 0
# k1 = 2 * pi / (lamb * 10 ** (-6))
# X = X1 * k1
#
# kx_crit_Au = (2 * pi / (lamb * 10 ** (-6))) * emath.sqrt(1 * get_dielectic(lamb,'Au') / (1 + get_dielectic(lamb,'Au')))
#
# kx_crit_Ag = (2 * pi / (lamb * 10 ** (-6))) * emath.sqrt(1 * get_dielectic(lamb,'Ag') / (1 + get_dielectic(lamb,'Ag')))
#
#
# plt.plot(X, calc_p_polaris_r(angle, lamb, X, 'None','Au')*np.conj(calc_p_polaris_r(angle, lamb, X, 'None','Au')), label=r'$|R^{p}|^2$ for Au', color='firebrick', lw=1.5)
# plt.axvline(np.real(kx_crit_Au), color='royalblue', linestyle='dashed', lw=1.5,
#             label=r'Crit $k_{1x}$ for Au:' + fr' ${round(np.real(kx_crit_Au / k1), 3)} \cdot k_{1}$ ')
#
# plt.plot(X, calc_p_polaris_r(angle, lamb, X, 'None','Ag')*np.conj(calc_p_polaris_r(angle, lamb, X, 'None','Ag')), label=r'$|R^{p}|^2$ for Ag', color='forestgreen', lw=1.5)
# plt.axvline(np.real(kx_crit_Ag), color='darkorange', linestyle='dashed', lw=1.5,
#             label=r'Crit $k_{1x}$ for Ag:' + fr' ${round(np.real(kx_crit_Ag/ k1), 3)} \cdot k_{1}$ ')
#
# plt.legend()
# plt.grid()
# ax5.set_xlabel(r'$k_{1x}$')
# ax5.set_xlim(xmin=a*k1)
# plt.hlines(y=0, xmin=a*k1, xmax=k1, color='g', linestyle='-', linewidth=1.5)
# plt.hlines(y=0, xmin=k1, xmax=b*k1, color='r', linestyle='-', linewidth=1.5)
# ax5.set_ylim(ymax=41, ymin=-2)
#
# plt.xticks(np.arange(a*k1, b * k1 + 1, k1/8))
# ax5.text((a+0.01) * k1, -1.5, 'propagating', style='italic')
# ax5.text(1.26 * k1, -1.5, 'evanescent', style='italic')
# plt.xticks(plt.xticks()[0], [r"$" + format(r / k1, ".2g") + r"k_1$" for r in plt.xticks()[0]])
# plt.savefig(f"Plot_abs_Rp_2_compramation_evanicent_{lamb}.pdf", format="pdf", bbox_inches="tight")
# plt.show()
