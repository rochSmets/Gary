
import numpy as np
from scipy.special import iv
from plasmapy.dispersion import plasma_dispersion_func as Z
import nlopt
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
        #print(m, k, rho, l, w_, w_+m*W, z)
        #print(m, iv(m, l), Z(z), sum)
        #print(m, z, Z(z), -1/z)
    return 1.0+w_*np.exp(-l)*sum/(np.sqrt(2)*np.abs(k_para)*vth)


def oneBranch(opt, kk, mode, T, theta):
    ww = np.empty_like(kk, dtype=complex)
    x0 = (mode, 0.1)

    for i, k in enumerate(kk):

        def objective(x, grad):
            v_ = np.abs(susceptibility(x, k, T, theta))
            return v_

        opt.set_min_objective(objective)
        x = opt.optimize(x0)
        ww[i] = complex(x[0], x[1])
        x0 = x
        print(k, susceptibility(x, k, T, theta))
    return ww


def main():
    n = 2 # nbr of free parameters : real frequency and growth rate
    mode = 2
    precision = 0.001
    T = 0.4
    theta = 85

    opt = nlopt.opt(nlopt.LN_NELDERMEAD, n)
    # opt = nlopt.opt(nlopt.LN_COBYLA, n)

    opt.set_lower_bounds([mode-1, -0.1])
    opt.set_upper_bounds([mode+1, +0.2])

    opt.set_xtol_rel(precision)

    kk = np.arange(5.1, 0.1, -0.1)
    ww = oneBranch(opt, kk = kk, mode = mode, T = T, theta = theta)
    plt.plot(kk, ww.real, c='k')
    plt.plot(kk, ww.imag, c='r')
    plt.savefig('t2.pdf')

    k = 5.0
    w = 2.0
    g = 0.0

    print('this is the end...')

if __name__=='__main__':
   main()

