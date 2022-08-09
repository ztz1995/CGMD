import matplotlib.pyplot as plt
import pickle as pkl
from scipy.optimize import curve_fit
from pyMD import functions as func
from sklearn.neighbors import KernelDensity
import numpy as np


if __name__ == '__main__':
    frames = 1000
    steps = 10000
    start = 0
    param = dict()
    fit_f = func.single_gaussian
    # fit_f = func.pdf_cauchy

    for n in range(1, 15):
        with open("data/EndToEnd_n=%d_%d_%d" % (n, frames, steps), "rb") as file:
            dist = pkl.load(file)
        values, bins, _ = plt.hist(dist, bins=80, range=[0, max(dist)], density=True)
        # plt.cla()
        # x = bins[:-1] + (bins[1] - bins[0]) / 2
        # popt, _ = curve_fit(fit_f, x, values, p0=[1, x[values.argmax()], 1, 0])
        # plt.plot(x, fit_f(x, *popt))
        # plt.plot(x, values)
        # print(n)
        # print(popt)
        # plt.title(str(n))
        # plt.show()
        # param[n] = list(popt)
        X = np.asarray(dist)[:, np.newaxis]
        kde = KernelDensity(kernel='gaussian', bandwidth=0.2).fit(X)
        X_plot = (bins[:-1] + (bins[1] - bins[0]) / 2)[:, np.newaxis]
        log_dens = kde.score_samples(X_plot)
        plt.plot(X_plot[:, 0], np.exp(log_dens))
        plt.show()
        param[n] = kde

    with open("data/EndToEnd_kernel.pkl" , "wb") as file:
        pkl.dump(param, file)

