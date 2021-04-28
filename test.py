
import numpy as np
from scipy.special import iv
from plasmapy.dispersion import plasma_dispersion_func as Z
from scipy.optimize import root
import matplotlib.pyplot as plt



def epsilon(w, k, T, theta, cVa):
    m=1
    e=1
    B=1
    m_max=6
    k_para = k*np.cos(theta*np.pi/180)
    k_perp = k*np.sin(theta*np.pi/180)
    w_ = complex(w[0], w[1])
    vth = np.sqrt(T/m)
    rho = m*vth/(e*B)
    W = e*B/m
    k_debye = cVa*(W/vth)**2
    l = (k_perp*rho)**2
    sum = 0
    for m in np.arange(-m_max, m_max+1, 1):
        z = (w_+m*W)/(np.sqrt(2)*np.abs(k_para)*vth)
        sum += iv(m, l)*Z(z)
    s_ = 1.0+w_*np.exp(-l)*sum/(np.sqrt(2)*np.abs(k_para)*vth)
    e_ = 1.0+s_*(k_debye/k)**2
    return [e_.real, e_.imag]


def solver(k, T, theta, cVa, x0):
    sol=root(epsilon, x0=x0, args=(k, T, theta, cVa), method='hybr') # hybr
    if sol.success != True :
        print('k : {:4.2f}  -  success : {} - x : {}'.format(k, sol.success, sol.x))
    return sol


def oneBranch(kk, x0, T, theta, cVa):
    ww=np.empty_like(kk, dtype=complex)
    for i, k in enumerate(kk):
        sol=solver(k, T, theta, cVa, x0=x0)
        ww[i]=complex(sol.x[0], sol.x[1])
        x0=sol.x
    return ww


def main():
    kk=np.arange(6.1, 0.1, -0.1)
    w1=oneBranch(kk=kk, x0=[0.5, -0.1], T=0.1, theta=85, cVa=100)
    w2=oneBranch(kk=kk, x0=[1.5, -0.1], T=0.1, theta=85, cVa=100)
    w3=oneBranch(kk=kk, x0=[2.5, -0.1], T=0.1, theta=85, cVa=100)
    w4=oneBranch(kk=kk, x0=[3.5, -0.1], T=0.1, theta=85, cVa=100)
    plt.plot(kk, w1.real, c='k')
    plt.plot(kk, w2.real, c='k')
    plt.plot(kk, w3.real, c='k')
    plt.plot(kk, w4.real, c='k')
    plt.xlabel('$k V_a/\Omega_p$')
    plt.ylabel('$\omega/\Omega_p$')
    plt.savefig('ionBernstein.pdf')


if __name__=='__main__':
   main()

