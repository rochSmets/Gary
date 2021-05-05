
import numpy as np
from scipy.special import iv
import matplotlib.pyplot as plt


multiply by exp - lambda


def main():
    l=np.arange(0.0, 4.0, 0.1)
    for m in range(5):
        plt.plot(l, iv(m, l), label = 'm={}'.format(m))
    plt.title('Modified Bessel function of order $m$ (even w. $m$)')
    plt.xlabel('$\lambda$')
    plt.ylabel('$I_m(\lambda$)')
    plt.legend(loc = 'upper left')
    plt.savefig('Im.pdf')


if __name__ == '__main__':
    main()

