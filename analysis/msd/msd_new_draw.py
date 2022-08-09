import matplotlib.pyplot as plt
import pandas as pd
import pickle as pkl
import numpy as np


if __name__ == '__main__':

    plt.style.use(["science", "ieee"])
    # plt.rc("text", usetex=True)
    plt.figure(dpi=300)

    labels = ["AA", "CGU", "CG"]
    msd = dict()

    for i in range(3):
        with open("data/msd_cm_%s_%d.pkl" % (labels[i].lower(), 100), "rb") as f:
            msd[labels[i]] = pkl.load(f)

    # plt.plot(msd["AA"]["time"], msd["CGU"]["cm"]/msd["AA"]["cm"],label="CGU/AA")
    # plt.plot(msd["AA"]["time"], msd["CG"]["cm"]/msd["AA"]["cm"], label="CG/AA")
    _type = "cm"
    plt.plot(msd["AA"]["time"], msd["AA"][_type], label="AA", marker="o", markerfacecolor='none')
    plt.plot(msd["AA"]["time"], msd["CGU"][_type], label="CGU", marker="v", markerfacecolor='none')
    plt.plot(msd["AA"]["time"], msd["CG"][_type], label="CG", marker="^", markerfacecolor='none')
    plt.plot(msd["AA"]["time"], msd["AA"]["time"]*1E-6)

    plt.legend()
    plt.xlabel("Time (fs)")
    plt.ylabel("MSD")
    plt.xscale("log")
    plt.yscale("log")
    plt.show()