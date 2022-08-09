import matplotlib.pyplot as plt
import pandas as pd
import pickle as pkl
import numpy as np
from scipy.optimize import curve_fit


def cc_double_relax(x, a, b, t1, t2):
    return a * np.exp(-x / t1) + b * np.exp(-x / t2)


def linear(x, k, b):
    return k * x + b


if __name__ == '__main__':
    # plt.style.use("science")
    plt.style.use(["science", 'ieee', 'no-latex'])
    plt.rc("font",family='Times New Roman')
    # plt.rc("text", usetex=True)
    plt.figure(dpi=1000)

    y_dict = dict()
    # AA
    cut_off = 5
    title = "aa"
    label = "npt_7"
    frames = 4000
    steps = 10000
    nums = 2000
    start_delta = 40000
    with open("data/urea_%s_%d_%d_%d_%d.pkl" % (title, steps, frames, nums, start_delta), "rb") as f:
        y_dict["AA"] = pkl.load(f)

    # CGU new
    title = "cgu"
    frames = 4000
    steps = 10000
    nums = 3000
    start_delta = 40000
    with open("data/urea_%s_%d_%d_%d_%d.pkl" % (title, steps, frames, nums, start_delta), "rb") as f:
        y_dict["CGU"] = pkl.load(f)

    # # CG new
    # title = "cg"
    # frames = 4000
    # steps = 10000
    # nums = 2000
    # start_delta = 40000
    # with open("data/urea_%s_%d_%d_%d_%d.pkl" % (title, steps, frames, nums, start_delta), "rb") as f:
    #     y_dict["CG"] = pkl.load(f)

    # x = np.arange(4000) * steps / 1000000
    x = np.arange(3000) * steps / 1000000

    # plt.gca().set_yscale('log')
    ecolor = {"AA": "black", "CGU": "red", "CG":"blue"}
    marker = {"AA": "o", "CGU": "s", "CG":"v"}
    for k, v in y_dict.items():
        y = v[:len(x)]
        plt.scatter(x[::100], y[::100], label=k, color='', edgecolor=ecolor[k], s=6, marker=marker[k])
        fit_func = cc_double_relax
        p0 = [0.5, 0.5, 0.01, 100000]
        bound_l = [0, 0, 0, 0]
        bound_r = [np.inf, np.inf, np.inf, np.inf]
        popt, _ = curve_fit(fit_func, x[1:], y[1:], p0, bounds=(bound_l, bound_r))
        plt.plot(x, fit_func(x, *popt), label=k + " fit")
        print(popt)

    plt.xlabel('Time (ns)', fontsize=14)
    plt.ylabel('$C_U(t, r_C)$', fontsize=14)
    # plt.legend(loc=(0.65, 0.35))
    plt.legend(loc="upper right")
    plt.yscale('log')
    # plt.ylim(0.0015, 15)
    # plt.yticks(np.arange(2, 10)/10, np.arange(2, 10)/10)
    # from matplotlib.ticker import MultipleLocator, FuncFormatter
    # plt.gca().yaxis.set_minor_locator(MultipleLocator(0.01))
    # plt.setp(plt.gca().get_yminorticklabels(), visible=False)
    plt.show()

    # for k, v in y_dict.items():
    #     print(k)
    #     y = np.log(v[:len(x)])
    #     # y = np.log(y)
    #     fit_func = linear
    #     p0 = [-1, 1]
    #     popt, _ = curve_fit(fit_func, x[-2000:-1000], y[-2000:-1000], p0)
    #     print(popt)
    #     print(-1/popt[0])
    #     plt.plot(x, y)
    #     plt.plot(x, linear(x, *popt))
    # plt.show()