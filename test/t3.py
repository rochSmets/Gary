
import numpy as np
from scipy.special import iv
from plasmapy.dispersion import plasma_dispersion_func as Z
from scipy.optimize import root
import matplotlib.pyplot as plt



def susceptibility(w, k, T, theta): # w is a list or tuple  ---  return val is complex
    m =1
    e = 1
    B = 1
    m_max = 8
    k_para = k*np.cos(theta*np.pi/180)
    k_perp = k*np.sin(theta*np.pi/180)
    w_ = complex(w[0], w[1])
    vth = np.sqrt(T/m)
    rho = m*vth/(e*B)
    W = e*B/m
    l = (k_perp*rho)**2
    sum = 0
    for m in np.arange(-m_max, m_max+1, 1):
        z = (w_+m*W)/(np.sqrt(2)*np.abs(k_para)*vth)
        sum += iv(m, l)*Z(z)
    return 1.0+w_*np.exp(-l)*sum/(np.sqrt(2)*np.abs(k_para)*vth)


def oneBranch(kk, mode, T, theta):
    ww = np.empty_like(kk, dtype=complex)
    x0 = (mode, 0.1)

    for i, k in enumerate(kk):

        def objective(w):
            v_ = susceptibility(w, k, T, theta)
            return v_.real, v_.imag

        sol = root(objective, x0 = x0, method='hybr')
        ww[i] = complex(sol[0], sol[1])
        x0 = sol
        print(k, susceptibility(sol, k, T, theta))
    return ww


def main():
    n = 2 # nbr of free parameters : real frequency and growth rate
    mode = 2
    T = 0.4
    theta = 85


    kk = np.arange(5.1, 0.1, -0.1)
    ww = oneBranch(kk = kk, mode = mode, T = T, theta = theta)
    plt.plot(kk, ww.real, c='k')
    plt.plot(kk, ww.imag, c='r')
    plt.savefig('t2.pdf')

    k = 5.0
    w = 2.0
    g = 0.0

    print('this is the end...')

if __name__=='__main__':
   main()

