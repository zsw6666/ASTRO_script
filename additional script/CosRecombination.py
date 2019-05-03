import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import astropy.units as u
import astropy.constants as const


def Saha(T):
    A1 = (4 * np.sqrt(2) * 1.202057) / np.sqrt(np.pi);
    A2 = (4.1e-10)
    s1 = (((const.k_B * T).to(u.J) / (const.m_e * const.c ** 2).to(u.J)).value) ** 1.5
    s2 = np.exp(((13.6 * u.eV).to(u.J) / (const.k_B * T).to(u.J)).value)

    def infunc(Xe):
        f = (1 - Xe) - A1 * A2 * s1 * s2 * (Xe ** 2)
        return f

    return infunc


def plot(x, y, xlabel, ylabel, title=None):
    plt.plot(x, y)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.show()


def Run():
    T_set = np.linspace(5000, 1000, 5000)
    Xe_set = []
    for T in T_set:
        T = T * u.K
        Xe = abs(fsolve(Saha(T), 0.5))
        Xe_set.append(Xe)
    Xe_set = np.array(Xe_set)
    plot(T_set, Xe_set, xlabel=r'$log(T) \ K$', ylabel=r'$log(X_{e})$', title=r'$X_{e}-T \ relation$')


Run()