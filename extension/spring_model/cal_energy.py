import pyMD.monte_carlo as mc
import numpy as np
import matplotlib.pyplot as plt
from pyMD.functions import fene_expand
from scipy.optimize import curve_fit


if __name__ == '__main__':

    chain_model = mc.TOChain()
    # print(chain_model.param)
    logp = list()
    rs = np.arange(20., 40., 0.5)
    r2s = np.square(rs)
    for r2 in r2s:
        logp.append(chain_model.probability(14, r2))
    logp = np.asarray(logp)
    kb = 1.9855487E-3
    T = 300
    # u = -kb * T * (logp - 2 * np.log(rs))
    u = -kb * T * logp
    u -= u.min()
    p0 = [1, 30, 30]
    bl = [0, 1, 1]
    br = [np.inf, np.inf, np.inf]
    p, _ = curve_fit(fene_expand, rs, u, p0, bounds=(bl, br))
    plt.plot(rs, u)
    plt.plot(rs, fene_expand(rs, *p))
    plt.show()
    print(p)
    # [1.20593084e-04 2.86654111e+01 1.82028604e+04]