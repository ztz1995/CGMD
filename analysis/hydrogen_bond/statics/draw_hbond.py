import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt


def customize_violin(parts, color=None, label=""):
    for key in parts:
        if key == "bodies":
            for i, pc in enumerate(parts[key]):
                if color is not None:
                    pc.set_facecolor(color)  # '#D43F3A'
                pc.set_alpha(0.6)
                pc.set_edgecolor('black')
                if i == 0:
                    pc.set_label(label)
                # pc.set_edgecolor('#D43F3A')
        else:
            parts[key].set_color('black')


if __name__ == '__main__':
    plt.style.use(["science", "ieee", "no-latex"])
    plt.rc("font", family='Times New Roman')

    font_size = 12
    # plt.style.use(["science", "ieee"])
    # h1 = list()
    # hs = list()
    # h2 = list()
    # title = "aa"
    # h_cut = 2.5
    # for i in [0, 1, 2, 3, 5, 7, 8, 9]:
    #     with open("data/%s_cut_%.1f_%d.pkl" % (title, h_cut, i), "rb") as f:
    #         h_dict = pkl.load(f)
    #     h1 += list(h_dict["h1"])
    #     h2 += list(h_dict["h2"])
    #     hs += list(h_dict["hs"])
    #     # mea = np.array(hs).mean()
    #     # print
    # aa_dict = {"h1": np.asarray(h1), "h2": np.asarray(h2), "hs": np.asarray(hs)}

    title = "aa"
    with open("data/hbond_new_%s.pkl" % title, "rb") as file:
        hb1, hb2 = pkl.load(file)
    hb1 = np.asarray(hb1) / 50
    hb2 = np.asarray(hb2) / 50
    aa_dict = {"h1": hb1, "h2": hb2, "hs": hb1 + hb2 * 2}

    title = "cgu"
    with open("data/hbond_new_%s_50.pkl" % title, "rb") as file:
        hb1, hb2 = pkl.load(file)
    hb1 = np.asarray(hb1)/ 50
    hb2 = np.asarray(hb2)/ 50
    cgu_dict = {"h1": hb1, "h2": hb2, "hs": hb1 + hb2 * 2}

    # with open("data/%s_cut_%.1f_%d.pkl" % ("cgu_7blk", h_cut, 0), "rb") as f:
    #     cgu_7_dict = pkl.load(f)

    figure = plt.figure(dpi=1000, figsize=(4, 3))
    # axs = plt.subplot()
    axs = plt.gca()
    parts = axs.violinplot([aa_dict["h1"], aa_dict["h2"], aa_dict["hs"]], showmeans=True, positions=[0, 3, 6])
    # customize_violin(parts, '#ecebeb', label="AA")
    customize_violin(parts, 'dimgray', label="AA")

    parts = axs.violinplot([cgu_dict["h1"], cgu_dict["h2"], cgu_dict["hs"]], showmeans=True, positions=[1, 4, 7])
    # customize_violin(parts, '#D43F3A', label="CGU-1blk")
    # customize_violin(parts, '#7c7c7c', label="CGU")
    customize_violin(parts, 'red', label="CGU")

    # parts = axs.violinplot([cgu_7_dict["h1"], cgu_7_dict["h2"], cgu_7_dict["hs"]], showmeans=True, positions=range(0+2, 9+2, 4))
    # # customize_violin(parts, 'green', label="CGU-7blk")
    # customize_violin(parts, '#0d0d0d', label="CGU-8blk")

    # axs.yaxis.grid(True)
    axs.set_xticks([0.5, 3.5, 6.5])
    plt.setp(axs, xticks=[0.5, 3.5, 6.5], xticklabels=['disordered', 'ordered', 'total'])
    plt.xticks(fontsize=font_size)
    plt.yticks(fontsize=font_size)
    leg = plt.legend(loc="upper left", fontsize=font_size, borderpad=False)
    leg.get_frame().set_linewidth(0.0)
    plt.xlabel("Different types of HBs", fontsize=font_size+2)
    plt.ylabel("N(HBs)", fontsize=font_size+2)
    from matplotlib.pylab import minorticks_off, tick_params
    tick_params(axis="x", length=0)
    minorticks_off()
    plt.show()
