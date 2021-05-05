
import numpy as np
from plasmapy.dispersion import plasma_dispersion_func as Z
import matplotlib.pyplot as plt
from Z import Z_taylor, Z_asympt


def Zprime(x):
    return -2*(1+x*Z(x))


def Zprime_taylor(x, terms = 3):
    return -2*(1+x*Z_taylor(x, terms))


def Zprime_asympt(x, terms = 3):
    return -2*(1+x*Z_asympt(x, terms))


def main():
    x=np.linspace(-5, 5, 200)
    plt.plot(x, Zprime(x).real, '-k', linewidth = 4, label = 'Re(Z)')
    plt.plot(x, Zprime(x).imag, ':r', linewidth = 4, label = 'Im(Z)')
    plt.plot(x[80:120], Zprime_taylor(x[80:120], terms = 1).real, '-c', label = 'Taylor')
    plt.plot(x[80:120], Zprime_taylor(x[80:120], terms = 1).imag, '-c')
    plt.plot(x[120:200], Zprime_asympt(x[120:200], terms = 1).real, '-m', label = 'asymptotic')
    plt.plot(x[120:200], Zprime_asympt(x[120:200], terms = 1).imag, '-m')
    plt.xlabel('$\zeta \\ (real)$')
    plt.ylabel('$Z^{\prime}(\zeta$)')
    plt.axvline(x = 1, c = 'k', ls = '--')
    plt.text(0.2, -2.1, '$\zeta = 1$')
    plt.legend(loc='upper right')
    plt.savefig('Zprime.pdf')

if __name__=='__main__':
   main()

