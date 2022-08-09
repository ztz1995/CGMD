import pickle as pkl
import matplotlib.pyplot as plt
import numpy as np

if __name__ == '__main__':
    plt.style.use(["science", "ieee", "no-latex"])
    plt.rc("font", family='Times New Roman')
    inch_cm = 1 / 2.54
    fig = plt.figure(dpi=1000, figsize=(8 * inch_cm, 6 * inch_cm))
    ax = fig.add_subplot(111)
    h_len = 3.
    h_ang = 120
    labels = ["AA", "CGU"]
    # with open("data/l_%s_%d_%d_%.1f_%d_500K_%d.pkl" % ("cgu", 1, 5000, h_len, h_ang, 1), "rb") as f:
    with open("data/l_%s_%d_%d_%.1f_%d_500K_%d.pkl" % ("CGU", 1, 100000, h_len, h_ang, 1), "rb") as f:
        times2 = pkl.load(f)

    with open("data/l_%s_%d_%d_%.1f_%d_500K_%d.pkl" % ("aa", 1, 100000, h_len, h_ang, 0), "rb") as f:
        times1 = pkl.load(f)

    values1, bins = np.histogram(times1, bins=1000, range=[0, 2000], density=True)
    x = bins[:-1] + (bins[1] - bins[0]) / 2
    ax.plot(x, values1, label=labels[0])

    values2, bins = np.histogram(times2, bins=1000, range=[0, 2000], density=True)
    x = bins[:-1] + (bins[1] - bins[0]) / 2
    ax.plot(x, values2, label=labels[1])
    # print(x[values.argmax()])

    plt.legend(loc="upper right")
    # plt.xlim(0, 50)
    plt.ylim(bottom=0.)
    plt.ylabel("Probability", fontsize=14)
    plt.xlabel("Lifetime (fs)", fontsize=14)

    ax.annotate("", xy=(390, 0.01), xytext=(80, 0.005),arrowprops = dict(arrowstyle="->"))
    left_inset_ax = fig.add_axes([.45, .35, .4, .4])
    left_inset_ax.plot(x, values1)
    left_inset_ax.plot(x, values2)
    left_inset_ax.set(xlim=(0, 50), xticks=[0, 25, 50], yticks=[0, 0.01, 0.02,0.03])

    plt.show()
