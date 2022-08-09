import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from matplotlib import cm


# plt.style.use(["science", "ieee"])


def load_data(fp, _type=0):
    data = pd.read_csv(fp, sep="\s+", skiprows=1, header=None)
    if _type == 0:
        data.columns = ["step", "stress_x", "stress_y", "stress_z", "l_x", "l_y", "l_z", "temp"]
    elif _type == 1:
        data.columns = ["step", "stress_x", "stress_y", "stress_z", "ke_hard", "ke_soft", "virial_hard", "virial_soft", "pair_hard", "pair_soft"]
    elif _type == 2:
        data.columns = ["step", "stress_x", "stress_y", "stress_z", "ke_hard", "ke_soft", "virial_hard", "virial_soft", "pair_hard", "pair_soft",
                        "lx", "ly", "lz"]
    return data


def draw(x, y, window, label, normal=False):
    # y_mean = y[:1000].mean()
    y = y.rolling(window, min_periods=window, center=True).mean()
    # y /= y[window//2]
    # plt.plot((x/1E6)[::100000], y[::100000], label=label)
    y0 = y[window // 2]
    if normal:
        plt.plot((x / 1E6)[::100], y[::100] / y0, label=label)
    else:
        plt.plot((x / 1E6)[::100], y[::100], label=label)


def main():
    # plt.style.use(["science", "ieee"])
    # plt.rc("text", usetex=True)
    plt.figure(dpi=300)

    # file_path1 = "/home/centos/work/1blk_3200/cgu_8blk_400chain/relaxation_5_%s/8blk_400.relax.txt"
    # file_path1 = "/home/centos/work/1blk_3200/cgu_8blk_400chain/relaxation_10/8blk_400.relax.txt"
    # file_path1 = "/home/centos/work/1blk_3200/cgu_8blk_400chain/relaxation_20_%s/8blk_400.relax.txt"
    # file_path2 = "/home/centos/work/1blk_3200/cg_8blk_400chain/relaxation_20_%s/8blk_400.relax.txt"

    # file_path1 = "/home/centos/work/1blk_3200/cgu_8blk_400chain_new/relaxation_20_%s/8blk_400.relax.txt"
    # file_path2 = "/home/centos/work/1blk_3200/cg_8blk_400chain_new/relaxation_20_%s/8blk_400.relax.txt"
    file_path3 = "/home/centos/work/1blk_3200/cg_8blk_400chain_new/relax_npt_%s/8blk_400.relax.txt"
    file_path4 = "/home/centos/work/1blk_3200/cgu_8blk_400chain/relaxation_5_%s/8blk_400.relax.txt"
    file_path6 = "/home/centos/work/1blk_3200/cg_8blk_400chain/relaxation_5_%s/8blk_400.relax.txt"
    normal = False

    # dir_path_list = [file_path1, file_path2]
    # dir_path_list = [file_path3]
    dir_path_list = [file_path4, file_path6]
    # dir_path_list = [file_path2, file_path3]
    # dir_path_list = [dir_path1, dir_path2]
    # dir_path_list = [dir_path9, dir_path6, dir_path5, dir_path10]
    windows = [1000, 1000, 10000, 10000]
    labels = ["CGU", "CG"]
    linestyles = ["--", "-"]
    colors = ["red", "blue"]
    _types = [2, 2]
    ds = ["x"]
    # plt.yscale("log")

    # dir_path_list = [dir_path29, dir_path36]
    # colors = [cm.tab10(x) for x in np.linspace(0, 1, 10)][:]
    # print(colors)
    # plt.plot([0,1],[0.01,0.01])

    base_value = {'cgu_hard': {'x': 50.765709531785014, 'y': 49.135818963697254, 'z': 47.020965637464606},
                  'cgu_soft': {'x': -52.473627764230905, 'y': -50.60273012692953, 'z': -49.99395140978971},
                  'cg_hard': {'x': 51.74498946434563, 'y': 51.506702981299746, 'z': 50.08626441759919},
                  'cg_soft': {'x': -48.62960914596994, 'y': -47.78661370375951, 'z': -49.4993237109538}}

    for num, file_path in enumerate(dir_path_list):
        x = None
        y = None
        vi_hard = None
        vi_soft = None
        # for d in ["x", "y", "z"]:
        # for d in ["x", "y"]:
        for d in ds:
            print(d)
            data = load_data(file_path % d, _type=_types[num])
            # data2 = load_data(file_path % d + "2")
            # data = pd.concat([data,data2], axis=0)
            if x is None:
                x = data["step"]
                y = data["stress_x"]
                # y = data["stress_y"]
                # vi_hard = data["virial_hard"] + data["ke_hard"]
                # vi_soft = data["virial_soft"] + data["ke_soft"]
            else:
                y += data["stress_x"]
                # vi_hard += data["virial_hard"] + data["ke_hard"]
                # vi_soft += data["virial_soft"] + data["ke_soft"]
        # print(vi_hard)
        y = y * 1000
        # y /= 2
        # draw(x, y, windows[num], label=labels[num])
        draw(x, y, windows[num], label=labels[num], normal=normal)

    plt.xlabel("Time (ns)", fontsize=14)
    # plt.ylabel("G(t) (MPa)", fontsize=14)
    if normal:
        plt.ylabel("G(t) normalized", fontsize=14)
    else:
        plt.ylabel("G(t) (Mpa)", fontsize=14)

    # plt.title("With dihedral")
    # plt.ylim(bottom=0.1)
    # plt.ylim(0.15, 1.0)
    # plt.yscale("log")

    plt.legend(loc="upper right")
    plt.show()

    # z = np.polyfit(data_s1.strain[190:490], data_s1_x[190:490], 1)
    # print(z)


if __name__ == '__main__':
    main()
