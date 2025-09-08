import matplotlib.pyplot as plt
import numpy as np
import scipy.special
from scipy.integrate import quad

def term1(yp, w, y, kRs, ck):
    pre_all = 2 ** (4 / (1 + 3 * w)) * (1 + 3 * w) ** (-4 / (1 + 3 * w)) * yp ** (-4 / (1 + 3 * w)) * kRs**2 / ck**2
    nuo = (np.cos(2 * np.pi / (1 + 3 * w)) - np.cos((3 * np.pi * w - np.pi) / (1 + 3 * w))) ** 2
    cmc = -2 * np.cos(
        (3 * ck * np.pi * w - 2 * kRs * y ** ((1 + 3 * w) / 2)) / (
                ck + 3 * ck * w))
    cmc2 = 2 * np.sin((ck * np.pi - 2 * kRs * y ** ((1 + 3 * w) / 2)) / (ck + 3 * ck * w))
    A = np.pi / 4 * (1 + (6 - 6 * w) / (6 * w + 2)) - kRs * yp / ck
    Aprime = np.pi / 4 * (1 + (6 * w - 6) / (6 * w + 2)) - kRs * yp / ck
    B = kRs * 0 / (2 * ck)
    y1 = 0.5 * cmc ** 2 * np.cos(2 * A) + 0.5 * cmc2 ** 2 * np.cos(2 * Aprime) - cmc * cmc2 * np.cos(A + Aprime)
    y2 = 0.5 * ((cmc ** 2 + cmc2 ** 2) - 2*cmc * cmc2* np.cos(A - Aprime)) * np.cos(2 * B)
    deo = y1
    return deo * pre_all / nuo

def term2(yp, w, y, kRs, ck):
    pre_all = 2 ** (4 / (1 + 3 * w)) * (1 + 3 * w) ** (-4 / (1 + 3 * w)) * yp ** (-4 / (1 + 3 * w)) * kRs**2 / ck**2
    nuo = (np.cos(2 * np.pi / (1 + 3 * w)) - np.cos((3 * np.pi * w - np.pi) / (1 + 3 * w))) ** 2
    cmc = -2 * np.cos(
        (3 * ck * np.pi * w - 2 * kRs * y ** ((1 + 3 * w) / 2)) / (
                ck + 3 * ck * w))
    cmc2 = 2 * np.sin((ck * np.pi - 2 * kRs * y ** ((1 + 3 * w) / 2)) / (ck + 3 * ck * w))
    A = np.pi / 4 * (1 + (6 - 6 * w) / (6 * w + 2)) - kRs * yp / ck
    Aprime = np.pi / 4 * (1 + (6 * w - 6) / (6 * w + 2)) - kRs * yp / ck
    B = kRs * 0 / (2 * ck)
    y2 = 0.5 * ((cmc ** 2 + cmc2 ** 2) - 2*cmc * cmc2* np.cos(A - Aprime)) * np.cos(2 * B)
    deo = y2
    return deo * pre_all / nuo

def intyp_term1(y, w):
    kRs = 10
    ck = 0.000791572
    intfun = quad(term1, 2 / (3 * w + 1), 2 / (3 * w + 1) * y ** ((3 * w + 1) / 2), args=(w, y, kRs, ck))
    return intfun[0] * ck ** 2 / kRs ** 2

def intyp_term2(y, w):
    kRs = 10
    ck = 0.000791572
    intfun = quad(term2, 2 / (3 * w + 1), 2 / (3 * w + 1) * y ** ((3 * w + 1) / 2), args=(w, y, kRs, ck))
    return intfun[0] * ck ** 2 / kRs ** 2


def avg_intyp_term1(w):
    intfun = quad(intyp_term1, 1, 500, args=(w))
    return intfun[0] /500



y = 50
w = 1/4
kRs = 10
ck = 0.000791572
y_list = np.linspace(1, 50, 500)
yp_list = np.linspace(2 / (3 * w + 1), 50, 2000)
term1_list = []
term2_list = []
# for i in range(len(y_list)):
#     term1_list.append(term1(yp_list[i], w, y, kRs, ck))
#     term2_list.append(term2(yp_list[i], w, y, kRs, ck))
# plt.plot(y_list, term1_list, 'b')
# plt.plot(y_list, term2_list, 'r')

plt.plot(yp_list, term1(yp_list, w, y, kRs, ck), 'b')
plt.plot(yp_list, term2(yp_list, w, y, kRs, ck), 'r')
# # plt.legend(['$\mathcal{G}_1$', '$\mathcal{G}_2$'])
plt.legend(['Integrand of $\mathcal{G}_1$', 'Integrand of $\mathcal{G}_2$'])
plt.title(r'$\theta = $' + str(w) + ', '+'y = ' + str(y))
plt.ylabel('Integrand')
plt.xlabel('$y_{+}$')
plt.xlim(1, 50)
plt.xscale('log')
plt.show()
# print(avg_intyp_term1(1/4))
# print(avg_intyp_term1(1/3))
# print(avg_intyp_term1(1/2))
# print(avg_intyp_term1(2/3))