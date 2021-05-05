
import numpy as np


def zeta(w, g, k, m, W, T):
    d_ = np.sqrt(2)*k*np.sqrt(T)
    return [(w+m*W)/d_, g/d_]


def main():
    T = 0.1
    W = 1.0
    ws = np.arange(0.8, 1.2, 0.1)
    gs = np.arange(-0.2, 0.2, 0.1)
    ks = np.arange(0.1, 4.0, 0.1)
    zeta_r, zeta_i = np.empty([0]), np.empty([0])
    wgk=[np.nan, np.nan, np.nan]
    for m in range(5):
        for w in ws:
            for g in gs:
                for k in ks:
                    zeta_r = np.append(zeta_r, zeta(w, g, k, m, W, T)[0])
                    zeta_i = np.append(zeta_i, zeta(w, g, k, m, W, T)[1])
                    wgk = np.vstack((wgk, [w, g, k]))
        print('m = {0:2d}  -  real : [{1:.2E} -- {2:.2E}] -- imag : [{3:.2E} -- {4:.2E}]'.format(m, zeta_r.min(), zeta_r.max(), zeta_i.min(), zeta_i.max()))
        mwgk_min = zeta_r.argmin()
        print('           wgk  : [{0:.2E}, {1:.2E}, {2:.2E} -- {3:.2E}, {4:.2E}, {5:.2E}]'.format(wgk[zeta_r.argmin()][0], wgk[zeta_r.argmin()][1], wgk[zeta_r.argmin()][2], wgk[zeta_r.argmax()][0], wgk[zeta_r.argmax()][1], wgk[zeta_r.argmax()][2]  ))


if __name__ == '__main__':
    main()

