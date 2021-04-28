
import numpy as np
from scipy.special import iv
from plasmapy.dispersion import plasma_dispersion_func as Z
import matplotlib.pyplot as plt
from scipy.optimize import fsolve



def Z_taylor(x):
    return (-2*x*(1-2*x**2/3+4*x**4/15))+1j*(np.sqrt(np.pi)*np.exp(-x**2))


def Z_asympt(x):
    sigma=np.empty_like(x)
    for i,val in enumerate(x):
        if val.imag > 0:
            sigma[i]=0
        if val.imag == 0:
            sigma[i]=1
        if val.imag < 0:
            sigma[i]=2
    return (-1/x*(1+1/(2*x**2)+3/(4*x**4)))+1j*(np.sqrt(np.pi)*sigma*np.exp(-x**2))


def main():
    x=np.linspace(-5, 5, 200)
    z=Z(x)
    plt.plot(x, z.real, '-k', label='Re(Z)')
    plt.plot(x, z.imag, ':r', label='Im(Z)')
    plt.plot(x[80:120], Z_taylor(x[80:120]).real, '-c', label='Taylor')
    plt.plot(x[80:120], Z_taylor(x[80:120]).imag, '-c')
    plt.plot(x[120:200], Z_asympt(x[120:200]).real, '-m', label='asymptotic')
    plt.plot(x[120:200], Z_asympt(x[120:200]).imag, '-m')
    plt.xlabel('$\zeta \\ (real)$')
    plt.ylabel('$Z(\zeta$)')
    plt.axvline(x=1, c='k', ls='--')
    plt.text(0.2, -1.5, '$\zeta = 1$')
    plt.legend(loc='upper right')
    plt.savefig('FriedConte.pdf')
    #plt.show()

if __name__=='__main__':
   main()

