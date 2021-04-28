
import numpy as np
from scipy.special import iv
from plasmapy.dispersion import plasma_dispersion_func as Z
from scipy.optimize import root
import matplotlib.pyplot as plt



def susceptibility(w, k, T):
    m=1
    e=1
    B=1
    m_max=9
    w_=complex(w[0], w[1])
    vth=np.sqrt(T/m)
    rho = m*vth/(e*B)
    W=e*B/m
    l=(k*rho)**2
    sum=0
    for m in np.arange(-m_max, m_max+1, 1):
        z=(w_+m*W)/(np.sqrt(2)*np.abs(k)*vth)
        sum+= iv(m, l)*Z(z)
    return 1.0+w_*np.exp(-l)*sum/(np.sqrt(2)*np.abs(k)*vth)


def zob(w, k, T):
    #return np.abs(susceptibility(w, k, T))
    s_=susceptibility(w, k, T)
    return (s_.real, s_.imag)


def bernstein(k, T, x0):
    sol=root(zob, x0=x0, args=(k, T), method='hybr')
    print('k : {:4.2f}  -  x : {}'.format(k, sol))
    return complex(sol.x[0], sol.x[1])


def oneBranch(kk, mode, T):
    ww=np.empty_like(kk, dtype=complex)
    x0=(mode, 0)
    for i, k in enumerate(kk):
        sol=bernstein(k, T, x0=x0)
        ww[i]=sol
    return ww


def main():
    kk=np.arange(5.1, 0.1, -0.1)
    ww=oneBranch(kk=kk, mode=0.1, T=0.2)
    plt.plot(kk, ww.real, c='k')
    plt.plot(kk, ww.imag, c='r')
    plt.show()


if __name__=='__main__':
   main()

