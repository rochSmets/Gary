
import numpy as np
from scipy.special import iv
from scipy.optimize import fsolve
import matplotlib.pyplot as plt



def susceptibilityCold(w, k, T):
    m=1
    e=1
    B=1
    m_max=9
    rho = m*np.sqrt(T/m)/(e*B)
    W=e*B/m
    l=(k*rho)**2
    sum=0
    for m in np.arange(-m_max, m_max+1, 1):
        sum+= iv(m, l)/(w+m*W)
    s_=1.0-w*np.exp(-l)*sum
    return s_


def bernstein(k, T, x0):
    sol=fsolve(susceptibilityCold, x0=x0, args=(k, T))
    return sol


def oneBranch(kk, mode, T):
    ww=np.empty_like(kk)
    x0=mode+0.1
    for i, k in enumerate(kk):
        sol=bernstein(k, T, x0=x0)
        ww[i]=sol
        #print(k, sol)
        x0=sol
    return ww


def main():
    kk=np.arange(5.1, 0.1, -0.1)
    w1=oneBranch(kk=kk, mode=1, T=0.2)
    w2=oneBranch(kk=kk, mode=2, T=0.2)
    w3=oneBranch(kk=kk, mode=3, T=0.2)
    plt.plot(kk, w1)
    plt.plot(kk, w2)
    plt.plot(kk, w3)
    plt.savefig('ionBernsteinCold.pdf')


if __name__=='__main__':
   main()

