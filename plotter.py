# S and P polatisation calculator
import matplotlib
from numpy import arcsin, sin, sqrt, pi, cos, arctan, arccos
import matplotlib.pyplot as plt
import numpy as np
from numpy.lib import emath
import os

os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin'


def get_dielectic(lamb):
    reflect_n = 0
    reflect_k = 0
    with open('refractiveindex/Au-real.csv', 'r+') as f1:
        real_part = f1.readlines()
        for i in range(2, len(real_part)):
            line_prev = list(map(float, real_part[i - 1].split(',')))
            line_cur = list(map(float, real_part[i].split(',')))
            if float(line_prev[0]) <= lamb <= float(line_cur[0]):
                if (line_prev[0] + abs((line_cur[0] - line_prev[0]) / 2)) >= lamb:
                    reflect_n = line_prev[1]
                else:
                    reflect_n = line_cur[1]

    with open('refractiveindex/Au-imag.csv', 'r+') as f2:
        imag_part = f2.readlines()
        for i in range(2, len(real_part)):
            line_prev = list(map(float, imag_part[i - 1].split(',')))
            line_cur = list(map(float, imag_part[i].split(',')))
            if float(line_prev[0]) <= lamb <= float(line_cur[0]):
                if (line_prev[0] + abs((line_cur[0] - line_prev[0]) / 2)) >= lamb:
                    reflect_k = line_prev[1]
                else:
                    reflect_k = line_cur[1]
                # reflect_k = (line_cur[1] - line_prev[1]) / (line_cur[0] - line_prev[0]) * (
                #         lamb - line_prev[0]) + line_prev[1]

    return (reflect_n + 1j * reflect_k) ** 2


def calc_p_polaris_r(theta, lamb):
    k1 = 2 * pi / lamb
    k1x = abs(k1) * sin(theta)
    k1z = abs(k1) * cos(theta)

    k2x = k1x

    abs_k2 = abs(emath.sqrt(get_dielectic(lamb))) / abs(1) * abs(k1)

    k2z = emath.sqrt(abs_k2 ** 2 - k2x ** 2)

    abs_rs_2 = (abs((get_dielectic(lamb) * k1z - 1 * k2z) / (get_dielectic(lamb) * k1z + 1 * k2z))) ** 2

    return abs_rs_2


plt.rcParams['text.usetex'] = True

# Lambda linspace
Y = np.linspace(0.1879, 1.9370, 100)

angle=pi/4
angle_deg=round(angle*360/(2*pi),0)

fig2 = plt.figure(2)
ax2 = plt.gca()



R_S = []

for i in Y:
    R_S.append(calc_p_polaris_r(angle, i))

res_lamb=max(Y[[i for i, x in enumerate(R_S) if x == max(R_S)]])


plt.plot(Y, R_S, label=r'$|R_{s}^{2}|$')
#plt.axvline(res_lamb, color='red', lw=1, label=r'$\lambda_{r}$ - resonance')
ax2.set_xlim(xmin=0)
ax2.set_ylim(ymin=0)
ax2.set_xlabel(r'$\lambda$ ($\mu m$)')
plt.title(r"Dependence $|R_{s}^{2}|$ from $\lambda$ for "+ fr"$\theta_1={angle_deg}^\circ$")

plt.legend()
plt.grid()
plt.savefig("Plot2.pdf", format="pdf", bbox_inches="tight")
plt.show()


# Angle linspace
X = np.linspace(0, pi / 2, 100)

lamb = 0.5821
theta0 = arccos(abs(emath.sqrt(1 - (get_dielectic(lamb) / (1 + get_dielectic(lamb))))))

theta0_deg=round(theta0*360/(2*pi),0)

fig1 = plt.figure(1)

ax1 = plt.gca()

plt.plot(X, calc_p_polaris_r(X, lamb), label=r'$|R_{s}^{2}|$')
plt.axvline(theta0, color='red', lw=1, label=fr'Resonance angle: ${theta0_deg}^\circ$ ')
plt.legend()
plt.grid()
plt.title(r"Dependence $|R_{s}^{2}|$ from $\theta_{1}$ for "+fr"$\lambda={lamb} \;\mu m$")
ax1.set_xlabel(r'Angle (rad)')

plt.savefig("Plot1.pdf", format="pdf", bbox_inches="tight")

ax1.set_xlim(xmin=0)
ax1.set_ylim(ymin=0)

#plt.xticks(plt.xticks()[0],[r"$" + format(r/np.pi, ".2g")+ r"\pi$" for r in plt.xticks()[0]])
plt.show()





