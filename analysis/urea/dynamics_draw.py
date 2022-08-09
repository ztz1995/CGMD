import matplotlib.pyplot as plt
import pandas as pd
import pickle as pkl
import numpy as np

def cc_double_relax(x, a, b, t1, t2):
    return a * np.exp(-x / t1) + b * np.exp(-x / t2)



if __name__ == '__main__':
    plt.style.use(["science", "ieee", "no-latex"])
    plt.rc("font", family='Times New Roman')
    # inch_cm = 1 / 2.54
    inch_cm = 1
    fig = plt.figure(dpi=300, figsize=(4 * inch_cm, 3 * inch_cm))

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
    nums = 2000
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

    x = np.arange(3000) * steps / 1000000

    # plt.gca().set_yscale('log')

    for k, v in y_dict.items():
        plt.plot(x, v[:len(x)], label=k)

    plt.xlabel('Time (ns)', fontsize=14)
    plt.ylabel('$C_U(t, r_C)$', fontsize=14)
    plt.legend()
    plt.yscale('log')
    # plt.ylim(0.0015, 15)
    # plt.yticks(np.arange(2, 10)/10, np.arange(2, 10)/10)
    # from matplotlib.ticker import MultipleLocator, FuncFormatter
    # plt.gca().yaxis.set_minor_locator(MultipleLocator(0.01))
    # plt.setp(plt.gca().get_yminorticklabels(), visible=False)
    plt.show()