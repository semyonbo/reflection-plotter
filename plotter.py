# S and P polatisation calculator
import matplotlib
from numpy import arcsin, sin, sqrt, pi, cos, arctan, arccos
import matplotlib.pyplot as plt
import numpy as np
from numpy.lib import emath


def get_dielectic(lamb):
    reflect_n = 0
    reflect_k = 0
    with open('refractiveindex/Au-real.csv', 'r+') as f1:
        real_part=f1.readlines()
        for i in range(2, len(real_part)):
            line_prev = list(map(float, real_part[i - 1].split(',')))
            line_cur = list(map(float, real_part[i].split(',')))
            if float(line_prev[0]) <= lamb <= float(line_cur[0]):
                reflect_n = (line_cur[1] - line_prev[1]) / (line_cur[0] - line_prev[0]) * (
                        lamb - line_prev[0]) + line_prev[1]

    with open('refractiveindex/Au-imag.csv','r+') as f2:
        imag_part=f2.readlines()
        for i in range(2, len(real_part)):
            line_prev = list(map(float, imag_part[i - 1].split(',')))
            line_cur = list(map(float, imag_part[i].split(',')))
            if float(line_prev[0]) <= lamb <= float(line_cur[0]):
                reflect_k = (line_cur[1] - line_prev[1]) / (line_cur[0] - line_prev[0]) * (
                        lamb - line_prev[0]) + line_prev[1]

    return (reflect_n+1j*reflect_k)**2

def calc_p_polaris_r(theta,lamb):
    k1=2*pi/lamb
    k1x=abs(k1)*sin(theta)
    k1z=abs(k1)*cos(theta)

    k2x=k1x

    abs_k2=abs(emath.sqrt(get_dielectic(lamb)))/abs(1)*abs(k1)

    k2z=emath.sqrt(abs_k2**2-k2x**2)

    abs_rs_2=(abs((get_dielectic(lamb)*k1z-1*k2z)/(get_dielectic(lamb)*k1z+1*k2z)))**2

    return abs_rs_2

X=np.linspace(0,pi/2,100)

lamb=0.6595

theta0=arccos(abs(emath.sqrt(1-(get_dielectic(lamb)/(1+get_dielectic(lamb))))))


plt.plot(X, calc_p_polaris_r(X,0.6595), label='|R_p|^2')
plt.axvline(theta0)
plt.axhline(0, color='black', lw=2)
plt.axvline(0, color='black', lw=2)



Y=np.linspace(0.1879,1.9370,100)

R_S=[]

#for i in Y:
#    R_S.append(calc_p_polaris_r(pi/3,i))

#plt.plot(Y, R_S, label='|R_s|^2')
#plt.axhline(0, color='black', lw=2)
#plt.axvline(0, color='black', lw=2)

plt.legend()
plt.grid()





#fig.suptitle('Коэф-ы Френеля для e1=' + str(ep1) + ' ,e2=' + str(ep2))
plt.show()
