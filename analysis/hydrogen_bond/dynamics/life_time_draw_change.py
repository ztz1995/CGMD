import pickle as pkl
import matplotlib.pyplot as plt
import numpy as np

if __name__ == '__main__':
    plt.style.use(["science", "ieee", "no-latex"])
    plt.rc("font", family='Times New Roman')
    inch_cm = 1 / 2.54
    fig, axes = plt.subplots(3, 2, dpi=1000, figsize=(4, 4.5))

    h_len = 3.
    h_ang = 120
    count = 0
    titles = ["Bond c1-o1", "Bond h1-n1", "Angle c1-n1-h1", "Angle n1-c1-o1", "Angle Ph-n1-h1"]
    labels = ["change_co", "change_nh", "change_cnh", "change_nco", "change_Phnh"]
    for i, label in enumerate(labels):
        title = titles[i]
        # with open("data/l_%s_%d_%d_%.1f_%d_500K_%d.pkl" % ("cgu", 1, 5000, h_len, h_ang, 1), "rb") as f:
        with open("data/l_%s_%d_%d_%.1f_%d_500K_%d.pkl" % ("CGU", 1, 100000, h_len, h_ang, 1), "rb") as f:
            times1 = pkl.load(f)

        with open("data/l_%s_%d_%d_%.1f_%d_500K_%s.pkl" % ("CGU", 1, 100000, h_len, h_ang, label), "rb") as f:
            times2 = pkl.load(f)

        a, b = count // 2, count % 2
        ax = axes[a][b]
        values1, bins = np.histogram(times1, bins=1000, range=[0, 2000], density=True)
        x = bins[:-1] + (bins[1] - bins[0]) / 2
        ax.plot(x, values1)

        values2, bins = np.histogram(times2, bins=1000, range=[0, 2000], density=True)
        x = bins[:-1] + (bins[1] - bins[0]) / 2
        ax.plot(x, values2)
        ax.set_xlim(0., 50.)
        ax.set_ylim(0., 0.05)
        ax.set_title(title)

        # plt.ylim(bottom=0.)
        # plt.ylabel("Probability", fontsize=14)
        # plt.xlabel("Lifetime (fs)", fontsize=14)

        # ax.annotate("", xy=(390, 0.01), xytext=(80, 0.005),arrowprops = dict(arrowstyle="->"))
        # left_inset_ax = fig.add_axes([.45, .35, .4, .4])
        # left_inset_ax.plot(x, values1)
        # left_inset_ax.plot(x, values2)
        # left_inset_ax.set(xlim=(0, 50), xticks=[0, 25, 50], yticks=[0, 0.01, 0.02,0.03])

        count += 1
    plt.show()
