import matplotlib.pyplot as plt
import pandas as pd
import pickle as pkl
import numpy as np


if __name__ == '__main__':
    plt.style.use(["science", "ieee", "no-latex"])
    plt.rc("font", family='Times New Roman')

    labels = ["AA", "CGU", "CG"]
    plt.figure(dpi=1000)

    for i in range(3):
        # with open("data/EndToEnd_%s" % labels[i].lower(), "rb") as file:
        with open("data/gyration_%s.pkl" % labels[i].lower(), "rb") as file:
            dists = pkl.load(file)
        # values, bins = np.histogram(dists, bins=40, range=[0, 80], density=True)
        values, bins = np.histogram(dists, bins=40, range=[5, 25], density=True)
        x = bins[:-1] + (bins[1] - bins[0])/2
        plt.plot(x, values, label=labels[i])

    plt.ylabel("Probability", fontsize=14)
    plt.xlabel("Length ($\mathrm{\AA}$)", fontsize=14)
    plt.legend(fontsize=10)
    plt.show()
