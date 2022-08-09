import pickle as pkl
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def cc_double_relax(x, a, b, t1, t2):
    return a * np.exp(-x / t1) + b * np.exp(-x / t2)


def cc_triple_relax(x, a, b, c, t1, t2, t3):
    return a * np.exp(-x / t1) + b * np.exp(-x / t2) + c * np.exp(-x / t3)


if __name__ == '__main__':
    # title = "cgu"
    # frames = 1000
    # steps = 20000
    # nums = 1000
    # with open("%s_%d_%d_%d_i.pkl" % (title, steps, frames, nums), "rb") as f:
    #     y = pkl.load(f)
    # t = np.arange(frames) * steps / 1000
    #
    # popt, _ = curve_fit(cc_triple_relax, t, y, [0.5,0.5,1.5, 500, 100000], bounds=([0.,0.,0, 100, 10000], [1,1,2.0, 1000, 1000000]))
    # # popt, _ = curve_fit(cc_double_relax, t, y, [0.5,0.01,100000], bounds=([0.1,0, 10000], [1,1.5, 1000000]))
    # print(popt)
    # plt.plot(t, y)
    # plt.plot(t, cc_triple_relax(t, *popt))
    # # plt.plot(t, cc_double_relax(t, *popt))
    # plt.ylim(0.5, 1)
    # plt.yscale("log")
    # plt.show()
    # # [2.88540699e-01 6.80250719e-02 3.39696475e-01 9.99951753e+02 1.49075543e+05]
    # # 0.288   0.349
    # # 0.0680  1000
    # # 0.644   1.49E5

    # title = "aa"
    # steps = 10000
    # frames = 1000
    # nums = 1000
    # with open("%s_%d_%d_%d_i.pkl" % (title, steps, frames, nums), "rb") as f:
    #     y = pkl.load(f)
    # t = np.arange(frames) * steps / 1000
    #
    # popt, _ = curve_fit(cc_triple_relax, t, y, [0.5,0.5,0.01, 500, 100000], bounds=([0.,0.,0, 100, 10000], [1,1,1.8, 1000, 1000000]))
    # print(popt)
    # plt.plot(t, y)
    # plt.plot(t, cc_triple_relax(t, *popt))
    # # plt.plot(t, cc_double_relax(t, *popt))
    # plt.ylim(0.5, 1)
    # plt.yscale("log")
    # plt.show()
    # # [2.34181935e-01 2.18207449e-02 1.50000000e+00 5.41253594e+02 2.72838045e+05]
    # # 0.234   1.8
    # # 0.0218  543
    # # 0.744   2.73E5
    y_dict = dict()

    # AA
    # h_len = 3.
    # h_ang = 120
    # title = "aa"
    # frames = 6000
    # steps = 10000
    # nums = 2000
    # start_delta = 40000
    # with open("data/i_%.1f_%d_%s_%d_%d_%d_%d_50chain.pkl" % (h_len, h_ang, title, steps, frames, nums, start_delta),
    #           "rb") as f:
    #     y_dict["AA"] = pkl.load(f)
    #
    # popt=[2.23036644e-02, 9.01740505e-01, 2.73777525e+00, 1.85752276e+03]


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
    # [3.32026075e-02 8.64006889e-01 6.92318022e+00 2.93338505e+02]
    # [2.07696420e-02 1.93652792e-01 6.85184881e-01 1.61498376e+00 4.11055814e+01 2.15749392e+08]
    # popt = [2.47725022e-02, 8.77766313e-01, 2.56836855e+00, 2.47526086e+02]

    t = np.arange(3000) * steps / 1000000
    y = y_dict[title.upper()][:len(t)]
    fit_func = cc_double_relax
    p0 = [0.5,0.5, 0.001, 100000]
    bound_l = [0, 0, 0, 0]
    bound_r = [np.inf, np.inf, np.inf, np.inf]
    # fit_func = cc_triple_relax
    # p0 = [0.3, 0.3, 0.3, 1, 1200, 1200]
    # bound_l = [0, 0, 0, 1, 1000, 1000]
    # bound_r = [np.inf, np.inf, np.inf, 100, 10000, np.inf]

    popt, _ = curve_fit(fit_func, t[1:], y[1:], p0, bounds=(bound_l, bound_r))
    print(popt)
    plt.plot(t, y)
    plt.plot(t, fit_func(t, *popt))
    # plt.plot(t, cc_double_relax(t, *popt))
    # plt.ylim(0.5, 1)
    plt.yscale("log")
    plt.show()
    # [2.48044032e-01 5.09156171e-02 7.58171943e+00 6.44523822e+03 7.05070822e+05]
    # 0.234   1.8
    # 0.0218  543
    # 0.744   2.73E5
