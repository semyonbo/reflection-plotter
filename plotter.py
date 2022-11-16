#P polatisation plotter
import matplotlib
from numpy import arcsin, sin, sqrt, pi, cos, arctan, arccos
import matplotlib.pyplot as plt
import numpy as np
from numpy.lib import emath
import os

os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin'


def get_dielectic(lamb):
    with open('refractiveindex/esAu.txt','r+') as file:
        lambs = list(map(lambda x: float(x), file.readline().split(',')))
        epsilons = list(map(lambda x: complex(x), file.readline().replace('i', 'j').split(',')))
        epsilon = epsilons[0]
        if len(lambs) == len(epsilons):
            for i in range(0, len(lambs)):
                if lambs[i-1] <= lamb*(10**(-6)) <= lambs[i]:
                    if (lambs[i-1] + abs((lambs[i] - lambs[i-1]) / 2)) >= lamb:
                        epsilon = epsilons[i-1]
                    else:
                        epsilon = epsilons[i]

    return epsilon


def calc_p_polaris_r(theta, lamb, k1x, type):
    k1 = 2 * pi / lamb
    if type == 'angle':
        k1x = k1 * sin(theta)

    k1z = emath.sqrt(k1 ** 2 - k1x ** 2)

    k2x = k1x

    k2 = emath.sqrt(get_dielectic(lamb)) * k1

    k2z = emath.sqrt(k2 ** 2 - k2x ** 2)

    rp = (1 * k2z - get_dielectic(lamb) * k1z) / (get_dielectic(lamb) * k1z + 1 * k2z)

    return rp


plt.rcParams['text.usetex'] = True

# Lambda linspace - plotting graph |R_p^2| (lambda)
Y = np.linspace(0.1879, 1.937, 100)

#Put angle here:
angle = pi / 4


angle_deg = round(angle * 360 / (2 * pi), 0)

fig2 = plt.figure(2)
ax2 = plt.gca()

R_S = []

for i in Y:
    R_S.append(emath.sqrt(calc_p_polaris_r(angle, i, 0, 'angle') * np.conj(calc_p_polaris_r(angle, i, 0, 'angle'))))

plt.plot(Y, np.real(R_S), label=r'$|R_{p}^{2}|$')
ax2.set_xlim(xmin=0)
# ax2.set_ylim(ymin=0)
ax2.set_xlabel(r'$\lambda$ ($\mu m$)')
plt.title(r"Dependence $|R_{p}^{2}|$ from $\lambda$ for " + fr"$\theta_1={angle_deg}^\circ$")

plt.legend()
plt.grid()
plt.savefig("Plot2.pdf", format="pdf", bbox_inches="tight")
plt.show()


#Plotting graph |R_p^2| (angle) for lamda:
# Angle linspace
X = np.linspace(0, pi / 2, 100)

#Put lambda here (in micrometre)
lamb = 0.6

theta0 = arccos(abs(emath.sqrt(1 - (get_dielectic(lamb) / (1 + get_dielectic(lamb))))))

theta0_deg = round(theta0 * 360 / (2 * pi), 0)

fig1 = plt.figure(1)

ax1 = plt.gca()

plt.plot(X, np.real(emath.sqrt(calc_p_polaris_r(X, lamb, 0, 'angle') * np.conj(calc_p_polaris_r(X, lamb, 0, 'angle')))),
         label=r'$|R_{p}^{2}|$')
plt.axvline(theta0, color='red', lw=1, label=fr'Resonance angle: ${theta0_deg}^\circ$ ')
plt.legend()
plt.grid()
plt.title(r"Dependence $|R_{p}^{2}|$ from $\theta_{1}$ for " + fr"$\lambda={lamb} \;\mu m$")
ax1.set_xlabel(r'Angle (rad)')

plt.savefig("Plot1.pdf", format="pdf", bbox_inches="tight")

ax1.set_xlim(xmin=0)
# ax1.set_ylim(ymin=0)

# plt.xticks(plt.xticks()[0],[r"$" + format(r/np.pi, ".2g")+ r"\pi$" for r in plt.xticks()[0]])
plt.show()


# Plot graph |R_p^2| from Array of k_x
#kx array starts from: a*k1
a=0

#kx array ends with: b*k1
b=5

#Amount of dots between a and b:
resolution=200

#Lambda of wave:
lamb = 0.6

X1 = np.linspace(a, b, resolution)
fig3 = plt.figure(3)
ax3 = plt.gca()
angle = 0
k1 = 2 * pi / lamb
X = X1 * k1

kx_crit = (2 * pi / lamb) * emath.sqrt(1 * get_dielectic(lamb) / (1 + get_dielectic(lamb)))

plt.plot(X, np.real(emath.sqrt(calc_p_polaris_r(angle, lamb, X, 'None') * np.conj(calc_p_polaris_r(angle, lamb, X, 'None')))),
         label=r'$|R_{p}^{2}|$')
plt.axvline(np.real(kx_crit), color='red',linestyle='dotted', lw=1, label=r'Critical $k_{1x}$:'+fr' ${round(np.real(kx_crit/k1), 3)} \cdot k_{1}$ ')
plt.legend()
plt.grid()
plt.title(
    r"Dependence $|R_{p}^{2}|$ from $k_{1x}$ for " + fr"$\lambda={lamb} \;\mu m$" + r" where $k_1=\frac{2 \pi}{\lambda}$")
ax3.set_xlabel(r'$k_{1x}$')
ax3.set_xlim(xmin=0)
plt.hlines(y = 0, xmin = 0, xmax = k1, color='g', linestyle='-',linewidth=1.5)
plt.hlines(y = 0, xmin = k1, xmax=5*k1, color='r', linestyle='-',linewidth=1.5)
# ax2.set_ylim(ymin=0)

plt.xticks(np.arange(0, b*k1+1, k1))
ax3.text(1, -0.75, 'propagating', style='italic')
ax3.text(4.1*k1, -0.75, 'evanescent', style='italic')
plt.xticks(plt.xticks()[0], [r"$" + format(r / k1, ".2g") + r"k_1$" for r in plt.xticks()[0]])
plt.savefig("Plot3.pdf", format="pdf", bbox_inches="tight")
plt.show()
