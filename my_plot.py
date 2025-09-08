import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad


def term2(y, w):
    kRs = 10
    ck = 0.000791572
    pre_all = 2 * (y**((3*w-3)/2)-1)/(3*w-3)
    nuo = (np.cos(2 * np.pi / (1 + 3 * w)) - np.cos((3 * np.pi * w - np.pi) / (1 + 3 * w))) ** 2
    cmc = -2 * np.cos(
        (3 * ck * np.pi * w - 2 * kRs * y ** ((1 + 3 * w) / 2)) / (
                ck + 3 * ck * w))
    cmc2 = 2 * np.sin((ck * np.pi - 2 * kRs * y ** ((1 + 3 * w) / 2)) / (ck + 3 * ck * w))
    A = np.pi / 4 * (1 + (6 - 6 * w) / (6 * w + 2))
    Aprime = np.pi / 4 * (1 + (6 * w - 6) / (6 * w + 2))
    y2 = ((cmc ** 2 + cmc2 ** 2) - 2*cmc * cmc2* np.cos(A - Aprime))
    deo = y2
    return deo * pre_all / nuo


y_list = np.linspace(1., 100, 3000)
w = [0, 1/10, 1/4, 1/3, 2/3]
plt.plot(y_list, term2(y_list, w[0]), 'cyan', linestyle='--')
plt.plot(y_list, term2(y_list, w[1]), 'fuchsia',linestyle='--')
plt.plot(y_list, term2(y_list, w[2]), 'orange',linestyle='--')
plt.plot(y_list, term2(y_list, w[3]), 'r',linestyle='--')
plt.plot(y_list, term2(y_list, w[4]), 'b',linestyle='--')
plt.legend([r'$\theta = 0$',r'$\theta = 1/10$', r'$\theta = 1/4$', r'$\theta = 1/3$', r'$\theta = 2/3$'])
plt.xlim(1, 100)
# plt.title(r'$\Upsilon$ for arbitrary expansion rate')
plt.xlabel('y')
plt.ylabel(r'$\Upsilon(y)$')
plt.xscale('log')
plt.show()

