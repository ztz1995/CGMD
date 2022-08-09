import pickle as pkl
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def func(x, k, t1, t2):
    return k * np.exp(- x/t1) + (1 - k) * np.exp(- x/t2)


def func2(x, k, t1):
    return k * np.exp(- x/t1) + (1 - k) * np.exp(- x/9.1391493)


def func3(x, t1, k):
    return k-x/t1


if __name__ == '__main__':

    title = "aa"
    frames = 1000
    step = 50
    nums = 100
    h_len = 3.0
    h_ang = 130

    with open("data/c_%s_%d_%d_%d_%.1f_%d.pkl" % (title, step, frames, nums, h_len, h_ang), "rb") as f:
        y = pkl.load(f)
    # y = np.log(y[500:900])

    x = np.arange(frames) * step / 1000
    # x = x[500:900]

    popt, _ = curve_fit(func, x, y, p0=[0.5,0.5,1])
    print(popt)
    plt.plot(x, y)
    plt.plot(x, func(x, *popt))
    plt.yscale("log")
    plt.show()

    # [ 9.1391493  -3.22769772]
