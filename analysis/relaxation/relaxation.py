import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from matplotlib import cm


def main():
    y0 = 0.
    z0 = 0.
    # plt.figure(figsize=(7, 5))
    plt.figure(figsize=(5.5, 4))
    process = "r1"

    dir_path1 = "/home/centos/work/1blk_50/cg_8blk_250chain/relaxation_5E-6_0.2/"
    dir_path2 = "/home/centos/work/1blk_50/cgu_8blk_250chain/relaxation_5E-6_0.2/"
    # fp = "/home/centos/work/1blk_3200/cgu_8blk_400chain/relaxation/8blk_400.relax.txt"

    stress_type = False
    strain_type = False

    # stress_type = True
    # strain_type = True

    dir_path_list = [dir_path1, dir_path2]
    # dir_path_list = [dir_path1, dir_path3]
    # dir_path_list = [dir_path1, dir_path2]
    # dir_path_list = [dir_path9, dir_path6, dir_path5, dir_path10]
    windows = [5000, 5000, 10000, 10000]
    labels = ["CG-8blk", "CGU-8blk"]
    linestyles = ["--", "-"]
    colors = ["blue", "red"]
    plt.yscale("log")

    # dir_path_list = [dir_path29, dir_path36]
    # colors = [cm.tab10(x) for x in np.linspace(0, 1, 10)][:]
    # print(colors)
    for num, dir_path in enumerate(dir_path_list):
        window = windows[num]
        csv_lines = list()
        file_path_1 = dir_path + "strain.%s_1.txt" % process
        if os.path.exists(file_path_1):
            with open(file_path_1, "r") as f:
                lines_1 = f.readlines()
                csv_lines += lines_1[1:]
        file_path_1 = dir_path + "strain.%s_2.txt" % process
        if os.path.exists(file_path_1):
            with open(file_path_1, "r") as f:
                lines_1 = f.readlines()
                csv_lines += lines_1[1:]
        file_path_1 = dir_path + "strain.%s.txt" % process
        if os.path.exists(file_path_1):
            with open(file_path_1, "r") as f:
                lines_1 = f.readlines()
                csv_lines += lines_1[1:]

        # print(csv_lines)
        with open("temp.csv", "w") as f:
            new_lines = ["strain stress_x stress_y stress_z l_x l_y l_z temp\n"] + csv_lines
            f.writelines(new_lines)
        data = pd.read_csv("temp.csv", delim_whitespace=True)
        y0 = data.l_y[0]
        z0 = data.l_z[0]
        y = data.l_y.rolling(window, min_periods=window, center=True).mean() / y0
        z = data.l_z.rolling(window, min_periods=window, center=True).mean() / z0

        if stress_type:
            stress = data.stress_x.rolling(window, min_periods=window, center=True).mean()
        else:
            # stress = data.stress_x.rolling(window, min_periods=window, center=True).mean() / stress_mean
            # print(stress)
            stress = data.stress_x.rolling(window, min_periods=window, center=True).mean()
            stress /= stress[window//2]

            # print(stress)
            plt.axhline(stress[len(stress)-window//2-1])

        strain = data.strain

        # label = "block-copolymer" if num ==0 else "pure soft"
        label = dir_path.split("/")[-1]
        # label = "CG_1E-6" if num==0 else "CGU_1E-7"
        if strain_type:
            # plt.plot(np.log(1 + strain), stress*1000, label=label, color=colors[num]) # , label=process
            if num == 1:
                plt.plot(np.log(1 + strain), stress * 1000, label=labels[num], color=colors[num], linewidth=2)  # , label=process
            else:
                plt.plot(np.log(1 + strain)[::500], (stress * 1000)[::500], linestyles[num], label=labels[num], linewidth=2,
                     color=colors[num])  # , label=process
        else:
            plt.plot(strain, stress, label=labels[num], color=colors[num])  # , label=process

    # plt.ylim(-0., 8)
    # plt.xlim(0.00, 2)
    # plt.ylim(bottom=0.)
    if strain_type:
        # plt.xlim(0., 1.)
        plt.ylim(bottom=0.)
        # plt.ylim(0., 200)
        strain_title = "True Strain"
        stress_title = "True Stress"
        # plt.ylim(0., 1000)
    else:
        # plt.xlim(0., 4)
        plt.ylim(bottom=0.01)
        # plt.ylim(0,150)
        strain_title = "Engineer Strain"
        stress_title = "Engineer Stress"
    # plt.vlines(0.2, 0,20)
    # plt.vlines(0.3, 0,20)
    # plt.vlines(0.4, 0,20)
    # plt.vlines(0.5, 0,20)
    # plt.vlines(0.6, 0,20)
    # plt.title(dir_path_list[-1].split("/")[-1])
    # plt.title("7blk_2000chain_1E-7")
    # plt.title("7blk_500chain_FENE_1E-6")
    # plt.title("5blk_500chain_FENE_1E-7")
    # plt.title(stress_title + " vs " + strain_title, fontsize=16)
    plt.xlabel("True Strain", fontsize=14)
    plt.ylabel("True Stress (MPa)", fontsize=14)
    # plt.title("With dihedral")
    # plt.xlim(0, 1.5)
    plt.legend()
    plt.show()

    # z = np.polyfit(data_s1.strain[190:490], data_s1_x[190:490], 1)
    # print(z)


if __name__ == '__main__':
    main()
