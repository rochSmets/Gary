
import numpy as np
from plasmapy.dispersion import plasma_dispersion_func as Z
import matplotlib.pyplot as plt



def Z_taylor(x, terms = 3):
    if (terms == 1):
        return -2*x                         +1j*(np.sqrt(np.pi)*np.exp(-x**2))
    elif (terms == 2):
        return (-2*x*(1-2*x**2/3))          +1j*(np.sqrt(np.pi)*np.exp(-x**2))
    elif (terms == 3):
        return (-2*x*(1-2*x**2/3+4*x**4/15))+1j*(np.sqrt(np.pi)*np.exp(-x**2))
    else:
        raise ValueError('terms can only be 1, 2 or 3 (default)')


def Z_asympt(x, terms = 3):
    sigma=np.empty_like(x)
    for i,val in enumerate(x):
        if (val.imag > 0):
            sigma[i] = 0
        elif (val.imag == 0):
            sigma[i] = 1
        elif (val.imag < 0):
            sigma[i] = 2
    if (terms == 1):
        return (-1/x)                          +1j*(np.sqrt(np.pi)*sigma*np.exp(-x**2))
    elif (terms == 2):
        return (-1/x*(1+1/(2*x**2)))           +1j*(np.sqrt(np.pi)*sigma*np.exp(-x**2))
    elif (terms == 3):
        return (-1/x*(1+1/(2*x**2)+3/(4*x**4)))+1j*(np.sqrt(np.pi)*sigma*np.exp(-x**2))
    else:
        raise ValueError('terms can only be 1, 2 or 3 (default)')


def main():
    x=np.linspace(-5, 5, 200)
    plt.plot(x, Z(x).real, '-k', linewidth = 4, label = 'Re(Z)')
    plt.plot(x, Z(x).imag, ':r', linewidth = 4, label = 'Im(Z)')
    plt.plot(x[80:120], Z_taylor(x[80:120], terms = 1).real, '-c', label = 'Taylor')
    plt.plot(x[80:120], Z_taylor(x[80:120], terms = 1).imag, '-c')
    plt.plot(x[120:200], Z_asympt(x[120:200], terms = 1).real, '-m', label = 'asymptotic')
    plt.plot(x[120:200], Z_asympt(x[120:200], terms = 1).imag, '-m')
    plt.xlabel('$\zeta \\ (real)$')
    plt.ylabel('$Z(\zeta$)')
    plt.axvline(x = 1, c = 'k', ls = '--')
    plt.text(0.2, -2.1, '$\zeta = 1$')
    plt.legend(loc='upper right')
    plt.savefig('Z.pdf')

if __name__=='__main__':
   main()

