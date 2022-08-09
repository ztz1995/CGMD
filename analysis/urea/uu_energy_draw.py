import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt


if __name__ == '__main__':

    plt.style.use(["science", "ieee", "no-latex"])
    plt.rc("font", family='Times New Roman')
    inch_cm = 1 / 2.54
    fig = plt.figure(dpi=1000, figsize=(8 * inch_cm, 6 * inch_cm))

    labels = ["CGU parallel", "CGU head-to-tail", "CG"]
    fps = ["data/uu_cgu_parallel.pkl", "data/uu_cgu_head-tail.pkl", "data/uu_cg.pkl"]

    for i in range(3):
        with open(fps[i], "rb") as file:
            x, y = pkl.load(file)
        plt.plot(x, y, label=labels[i])
        print(fps[i])
        print(y.min())
        print(x[y.argmin()])

    plt.ylim(-15, 25)
    plt.xlim(2, 14)
    plt.xlabel("Distance ($\mathrm{\AA}$)", fontsize=14)
    plt.ylabel("Potential (kcal/mol)", fontsize=14)
    plt.legend()
    plt.show()
