
import numpy as np
from scipy.special import iv
from plasmapy.dispersion import plasma_dispersion_func as Z
from scipy.optimize import root
import matplotlib.pyplot as plt



def susceptibility(w, k, T, theta):
    m=1
    e=1
    B=1
    m_max=6
    k_para = k*np.cos(theta*np.pi/180)
    k_perp = k*np.sin(theta*np.pi/180)
    w_=complex(w[0], w[1])
    vth=np.sqrt(T/m)
    rho = m*vth/(e*B)
    W=e*B/m
    l=(k_perp*rho)**2
    sum=0
    for m in np.arange(-m_max, m_max+1, 1):
        z=(w_+m*W)/(np.sqrt(2)*np.abs(k_para)*vth)
        sum+= iv(m, l)*Z(z)
    s_=1.0+w_*np.exp(-l)*sum/(np.sqrt(2)*np.abs(k_para)*vth)
    return [s_.real, s_.imag]


def solver(k, T, theta, x0):
    sol=root(susceptibility, x0=x0, args=(k, T, theta), method='hybr') # hybr
    if sol.success != True :
        print('k : {:4.2f}  -  success : {} - x : {}'.format(k, sol.success, sol.x))
    return sol


def oneBranch(kk, x0, T, theta):
    ww=np.empty_like(kk, dtype=complex)
    for i, k in enumerate(kk):
        sol=solver(k, T, theta, x0=x0)
        ww[i]=complex(sol.x[0], sol.x[1])
        x0=sol.x
    return ww


def main():
    kk=np.arange(6.1, 0.1, -0.1)
    w1=oneBranch(kk=kk, x0=[0.5, -0.1], T=0.1, theta=85)
    w2=oneBranch(kk=kk, x0=[1.5, -0.1], T=0.1, theta=85)
    w3=oneBranch(kk=kk, x0=[2.5, -0.1], T=0.1, theta=85)
    w4=oneBranch(kk=kk, x0=[3.5, -0.1], T=0.1, theta=85)
    plt.plot(kk, w1.real, c='k')
    #plt.plot(kk, w1.imag, c='r')
    plt.plot(kk, w2.real, c='k')
    #plt.plot(kk, w2.imag, c='r')
    plt.plot(kk, w3.real, c='k')
    #plt.plot(kk, w3.imag, c='r')
    plt.plot(kk, w4.real, c='k')
    #plt.plot(kk, w4.imag, c='r')
    plt.xlabel('$k V_a/\Omega_p$')
    plt.ylabel('$\omega/\Omega_p$')
    plt.savefig('ionBernstein.pdf')


if __name__=='__main__':
   main()

