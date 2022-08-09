import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from matplotlib import cm
import pickle as pkl


# plt.style.use(["science", "ieee"])


def load_data(fp, _type=0):
    data = pd.read_csv(fp, sep="\s+", skiprows=1, header=None)
    if _type == 0:
        data.columns = ["step", "stress_x", "stress_y", "stress_z", "l_x", "l_y", "l_z", "temp"]
    elif _type == 1:
        data.columns = ["step", "stress_x", "stress_y", "stress_z", "ke_hard", "ke_soft", "virial_hard", "virial_soft",
                        "pair_hard", "pair_soft"]
    elif _type == 2:
        data.columns = ["step", "stress_x", "stress_y", "stress_z", "ke_hard", "ke_soft", "virial_hard", "virial_soft",
                        "pair_hard", "pair_soft",
                        "lx", "ly", "lz"]
    return data


def draw(x, y, window, label, normal=False, **kwargs):
    y0 = y[:1000].mean()
    y = y.rolling(window, min_periods=window, center=True).mean()
    # y /= y[window//2]
    # plt.plot((x/1E6)[::100000], y[::100000], label=label)
    # y0 = y[window // 2]
    with open("data/%s.pkl" % label, "wb") as file:
        pkl.dump(np.asarray(y), file)

    if normal:
        plt.plot((x / 1E6)[::100], y[::100] / y0, label=label, **kwargs)
    else:
        plt.plot((x / 1E6)[::100], y[::100], label=label, **kwargs)


def main():
    plt.style.use(["science", "ieee", "no-latex"])
    plt.rc("font", family='Times New Roman')
    # plt.rc("text", usetex=True)
    plt.figure(dpi=1000)

    # file_path1 = "/home/centos/work/1blk_3200/cgu_8blk_400chain/relaxation_5/8blk_400.relax.txt"
    # file_path1 = "/home/centos/work/1blk_3200/cgu_8blk_400chain/relaxation_10/8blk_400.relax.txt"
    # file_path1 = "/home/centos/work/1blk_3200/cgu_8blk_400chain/relaxation_20_%s/8blk_400.relax.txt"
    # file_path2 = "/home/centos/work/1blk_3200/cg_8blk_400chain/relaxation_20_%s/8blk_400.relax.txt"

    file_path1 = "/home/centos/work/1blk_3200/cgu_8blk_400chain_new/relaxation_20_%s/8blk_400.relax.txt"
    file_path2 = "/home/centos/work/1blk_3200/cg_8blk_400chain_new/relaxation_20_%s/8blk_400.relax.txt"
    file_path3 = "/home/centos/work/1blk_3200/cg_8blk_400chain_new/relax_npt_%s/8blk_400.relax.txt"
    file_path5 = "/home/centos/work/1blk_3200/cgu_8blk_400chain_new/relax_npt_%s/8blk_400.relax.txt"
    file_path4 = "/home/centos/work/1blk_3200/cg_8blk_400chain/relaxation_5_%s/8blk_400.relax.txt"
    file_path6 = "/home/centos/work/1blk_3200/cgu_8blk_400chain/relaxation_5_%s/8blk_400.relax.txt"
    file_path7 = "/home/centos/work/1blk_3200/cg_8blk_400chain_new/relaxation_5_%s/8blk_400.relax.txt"
    file_path8 = "/home/centos/work/1blk_3200/cgu_8blk_400chain_new/relaxation_5_%s/8blk_400.relax.txt"

    # dir_path_list = [file_path2, file_path3]
    # dir_path_list = [file_path1, file_path5]
    # dir_path_list = [file_path4]
    # dir_path_list = [file_path4, file_path6]
    dir_path_list = [file_path7, file_path8]
    # dir_path_list = [file_path6, file_path4]
    # dir_path_list = [dir_path1, dir_path3]
    # dir_path_list = [file_path1, file_path2]
    windows = [1000, 1000, 10000, 10000]
    labels = ["CG", "CGU"]
    linestyles = ["-", "-"]
    colors = ["blue", "red"]
    # colors = ["red", "blue"]
    linewidth = 0.5
    # ds = ["x"]
    # ds = ["y"]
    ds = ["x", "y", "z"]
    _types = [2, 2]
    # normal = False
    normal = True
    # plt.yscale("log")

    # dir_path_list = [dir_path29, dir_path36]
    # colors = [cm.tab10(x) for x in np.linspace(0, 1, 10)][:]
    # print(colors)
    # plt.plot([0,1],[0.01,0.01])

    base_value = {'cgu_hard': {'x': 50.765709531785014, 'y': 49.135818963697254, 'z': 47.020965637464606},
                  'cgu_soft': {'x': -52.473627764230905, 'y': -50.60273012692953, 'z': -49.99395140978971},
                  'cgu_x': -0.00017305481490258333, 'cgu_y': -0.00014863477361451398, 'cgu_z': -0.0003012377833808436,
                  'cg_hard': {'x': 51.74498946434563, 'y': 51.506702981299746, 'z': 50.08626441759919},
                  'cg_soft': {'x': -48.62960914596994, 'y': -47.78661370375951, 'z': -49.4993237109538},
                  'cg_x': 0.0003156659107594189, 'cg_y': 0.0003769380460467666, 'cg_z': 5.9471767100845146e-05}
    # base_value = {'cgu_hard': {'x': 46.71195235304606, 'y': 46.47041905655383, 'z': 45.69100164015891},
    #               'cgu_soft': {'x': -58.04333093565993, 'y': -57.27057988524728, 'z': -56.01408430923586}, 'cgu_x': -0.001148151934883353,
    #               'cgu_y': -0.0010943262959673666, 'cgu_z': -0.0010459863514442292,
    #               'cg_hard': {'x': 48.291250050677114, 'y': 49.81049457523296, 'z': 48.366125924518556},
    #               'cg_soft': {'x': -53.85276241753623, 'y': -54.353069423872206, 'z': -54.63832790574639}, 'cg_x': -0.0005635202405719991,
    #               'cg_y': -0.000460276396538369, 'cg_z': -0.0006355308657479067}

    for num, file_path in enumerate(dir_path_list):
        x = None
        y = None
        vi_hard = None
        vi_soft = None
        m = 15000000
        for d in ds:
            print(d)
            data = load_data(file_path % d, _type=_types[num])
            # data2 = load_data(file_path % d + "2")
            # data = pd.concat([data,data2], axis=0)
            base_y = base_value[labels[num].lower() + "_" + d]
            base_hard = base_value[labels[num].lower() + "_hard"][d]
            base_soft = base_value[labels[num].lower() + "_soft"][d]
            m = min(m, len(data["step"]))
            if x is None:
                x = data["step"][:m]
                y = data["stress_%s" % d][:m] - base_y
                # print(y)
                vi_hard = (data["virial_hard"] + data["ke_hard"] - base_hard)[:m]
                vi_soft = (data["virial_soft"] + data["ke_soft"] - base_soft)[:m]
            else:
                x = data["step"][:m]
                y = (y[:m] + data["stress_%s" % d][:m] - base_y)
                vi_hard = vi_hard[:m] + (data["virial_hard"] + data["ke_hard"] - base_hard)[:m]
                vi_soft = vi_soft[:m] + (data["virial_soft"] + data["ke_soft"] - base_soft)[:m]

        # draw(x, y * 1000 / len(ds), windows[num], label=labels[num], normal=normal, linestyle=linestyles[num], color=colors[num], linewidth=linewidth)
        draw(x, vi_hard / 10 * 1.01325 / len(ds), windows[num], label=labels[num] + "_hard", normal=normal,
             linestyle=linestyles[num], color=colors[num], linewidth=linewidth)
        draw(x, vi_soft / 10 * 1.01325 / len(ds), windows[num], label=labels[num] + "_soft", normal=normal,
             linestyle=linestyles[num], color=colors[num], linewidth=linewidth)

    plt.xlabel("Time (ns)", fontsize=14)
    if normal:
        plt.ylabel("G(t) normalized", fontsize=14)
    else:
        plt.ylabel("G(t) (Mpa)", fontsize=14)
    # plt.title("With dihedral")
    # plt.ylim(bottom=0.1)
    # plt.ylim(0.15, 1.0)
    # plt.xscale("log")
    # plt.yscale("log")

    plt.legend(loc="upper right")
    plt.show()

    # z = np.polyfit(data_s1.strain[190:490], data_s1_x[190:490], 1)
    # print(z)


if __name__ == '__main__':
    main()
