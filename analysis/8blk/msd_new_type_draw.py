import matplotlib.pyplot as plt
import pandas as pd
import pickle as pkl
import numpy as np

if __name__ == '__main__':

    plt.style.use(["science", "ieee", "no-latex"])
    plt.rc("font", family='Times New Roman')
    # inch_cm = 1 / 2.54
    inch_cm = 1
    fig = plt.figure(dpi=1000, figsize=(4 * inch_cm, 3 * inch_cm))

    labels = ["AA", "CGU", "CG"]
    aver = [500, 500, 500]

    # for i in range(3):
    #     with open("data/msd_cm_type_es_%s_%d.pkl" % (labels[i].lower(), aver[i]), "rb") as f:
    #         msd[labels[i]] = pkl.load(f)
    # _type = "hard"
    # _type2 = "soft"
    # markers = {"AA":"o", "CGU":"v", "CG":"^"}
    # for label in labels:
    #     plt.plot(msd["AA"]["time"], msd[label][_type2] / msd[label][_type], label=label + " soft/hard", markerfacecolor='none',marker=markers[label])
    #     m = (msd[label][_type2] / msd[label][_type]).max()
    #     m_id = (msd[label][_type2] / msd[label][_type]).argmax()
    #     print(label, msd["AA"]["time"][m_id], m)

    _type = "hard"
    with open("../msd/data/msd_cm_type_es_%s_%d.pkl" % ("cg", 500), "rb") as f:
        msd = pkl.load(f)
    plt.plot(msd["time"], msd[_type], label="CG-1blk " + _type, markerfacecolor='none', marker="v")

    with open("data/msd_cm_%s_%d.pkl" % ("cg", 100), "rb") as f:
        msd = pkl.load(f)
    plt.plot(msd["time"], msd[_type], label="CG-8blk " + _type, markerfacecolor='none', marker="o")

    with open("../msd/data/msd_cm_type_es_%s_%d.pkl" % ("cgu", 500), "rb") as f:
        msd = pkl.load(f)
    plt.plot(msd["time"], msd[_type], label="CGU-1blk " + _type, markerfacecolor='none', marker="^")

    with open("data/msd_cm_%s_%d.pkl" % ("cgu", 100), "rb") as f:
        msd = pkl.load(f)
    plt.plot(msd["time"], msd[_type], label="CGU-8blk " + _type, markerfacecolor='none', marker="s")




    # print(label, msd["time"][m_id], m)


    # _type = "hard"
    # plt.plot([1,2], [0.001, 0.001])
    # plt.plot(msd["AA"]["time"], msd["CGU"][_type]/msd["AA"][_type],label="CGU/AA hard", marker="v", markerfacecolor='none')
    # plt.plot(msd["AA"]["time"], msd["CG"][_type]/msd["AA"][_type], label="CG/AA hard", marker="^", markerfacecolor='none')
    # _type = "soft"
    # plt.plot(msd["AA"]["time"], msd["AA"][_type], label="AA" + " " + _type, marker="o", markerfacecolor='none')
    # plt.plot(msd["AA"]["time"], msd["CGU"][_type], label="CGU" + " " + _type, marker="v", markerfacecolor='none')
    # plt.plot(msd["AA"]["time"], msd["CG"][_type], label="CG" + " " + _type, marker="^", markerfacecolor='none')
    # plt.plot(msd["AA"]["time"], msd["AA"]["time"] * 1E-6)

    plt.legend(fontsize=12)
    plt.xlabel("Time (fs)", fontsize=14)
    plt.ylabel("MSD", fontsize=14)
    plt.xscale("log")
    plt.yscale("log")
    # plt.ylim(bottom=0.5)
    # plt.axis('equal')
    plt.show()
