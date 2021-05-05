
import numpy as np
from scipy.special import ivp
import matplotlib.pyplot as plt

#scipy.special.ivp(v, z, n=1)[source]


def main():
    l=np.arange(0.0, 4.0, 0.1)
    for m in range(5):
        plt.plot(l, ivp(m, l), label = 'm={}'.format(m))
    plt.title('Derivative of modified Bessel function of order $m$ (even w. $m$)')
    plt.xlabel('$\lambda$')
    plt.ylabel('$I_m^{\prime}(\lambda$)')
    plt.legend(loc = 'upper left')
    plt.savefig('Imprime.pdf')


if __name__ == '__main__':
    main()

