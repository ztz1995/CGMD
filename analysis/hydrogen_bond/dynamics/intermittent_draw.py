import matplotlib.pyplot as plt
import pandas as pd
import pickle as pkl
import numpy as np

def cc_double_relax(x, a, b, t1, t2):
    return a * np.exp(-x / t1) + b * np.exp(-x / t2)


if __name__ == '__main__':

    plt.style.use(["science", "ieee", "no-latex"])
    plt.rc("font", family='Times New Roman')
    inch_cm = 1 / 2.54
    fig = plt.figure(dpi=1000, figsize=(8 * inch_cm, 6 * inch_cm))

    y_dict = dict()
    # AA
    h_len = 3.
    h_ang = 120
    title = "aa"
    frames = 4000
    steps = 10000
    nums = 2000
    start_delta = 40000
    with open("data/i_%.1f_%d_%s_%d_%d_%d_%d_50chain.pkl" % (h_len, h_ang, title, steps, frames, nums, start_delta),
              "rb") as f:
        y_dict["AA"] = pkl.load(f)

    # CGU new
    h_len = 3.
    h_ang = 120
    title = "cgu"
    frames = 10000
    steps = 10000
    nums = 2000
    start_delta = 40000
    with open("data/i_%.1f_%d_%s_%d_%d_%d_%d_160chain.pkl" % (h_len, h_ang, title, steps, frames, nums, start_delta),
              "rb") as f:
        y_dict["CGU"] = pkl.load(f)

    # # CG new
    # title = "cg"
    # frames = 10000
    # steps = 10000
    # nums = 2000
    # start_delta = 40000
    # with open("data/urea_%s_%d_%d_%d_%d.pkl" % (title, steps, frames, nums, start_delta), "rb") as f:
    #     y_dict["CG"] = pkl.load(f)

    x = np.arange(3000) * steps / 1000000

    # plt.yscale("log")
    plt.gca().set_yscale('log')
    ecolor = {"AA": "black", "CGU": "red"}
    marker = {"AA": "o", "CGU": "s"}
    popt = {"AA": [2.23036644e-02, 9.01740505e-01, 2.73777525e+00, 1.85752276e+03],
            "CGU": [2.47725022e-02, 8.77766313e-01, 2.56836855e+00, 2.47526086e+02]}

    for k, v in y_dict.items():
        y = v[:len(x)]
        plt.scatter(x[::100], y[::100], label=k, color='', edgecolor=ecolor[k], s=10, marker=marker[k])
        plt.plot(x[40:], cc_double_relax(x[40:], *popt[k]), label=k+" fit")

    plt.xlabel('Time (ns)', fontsize=14)
    plt.ylabel('$C_{HB}(t)$', fontsize=14)
    plt.plot([0, 1], [100, 100])

    plt.ylim(0.7, 1.05)
    plt.xlim(0, 30)
    plt.yticks([0.7, 0.8, 0.9, 1.], [0.7, 0.8, 0.9, 1.])
    # from matplotlib.ticker import MultipleLocator, FuncFormatter
    # plt.gca().yaxis.set_minor_locator(MultipleLocator(0.01))
    # plt.setp(plt.gca().get_yminorticklabels(), visible=False)
    plt.legend(loc="upper right")
    plt.show()
