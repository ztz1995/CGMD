import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from matplotlib import cm

outtype = 1


def plot_data(data, y0, z0, stress_type, strain_type, label, color, window, index=0):
    y = data.l_y.rolling(window, min_periods=window, center=True).mean() / y0
    z = data.l_z.rolling(window, min_periods=window, center=True).mean() / z0

    if stress_type:
        stress = data.stress_x.rolling(window, min_periods=window, center=True).mean()
        if outtype == 2:
            virial = data.virial.rolling(window, min_periods=window, center=True).mean()
            ke = data.ke.rolling(window, min_periods=window, center=True).mean()
    else:
        stress = data.stress_x.rolling(window, min_periods=window, center=True).mean() * y * z
    strain = data.strain
    if strain_type:
        plt.plot(np.log(1 + strain), stress * 1000, label=label, color=color, zorder=20-index)  # , label=process
        # plt.plot(np.log(1 + strain), virial * 1000, label=label, color=color)  # , label=process
        # plt.plot(np.log(1 + strain), ke * 1000, label=label, color=color)  # , label=process
    else:
        plt.plot(strain, stress * 1000, label=label, color=color)  # , label=process


def plot_relax(data, label, color, window):
    stress = data.stress_x.rolling(window, min_periods=window, center=True).mean()
    # stress = (data.s_virial-data.s_pair).rolling(window, min_periods=window, center=True).mean()
    # stress = data.s_virial.rolling(window, min_periods=window, center=True).mean()
    plt.plot(range(len(stress)), stress * 1000, label=label, color=color)


def main():
    plt.figure(figsize=(7, 5))
    # dir_path = "/home/centos/work/1blk_50/cg_7blk_300chain/1E-6_cycle/"
    # dir_path = "/home/centos/work/1blk_50/cg_7blk_300chain/5E-6_cycle/"
    # dir_path = "/home/centos/work/1blk_50/cg_8blk_250chain/1E-6/"
    # dir_path = "/home/centos/work/1blk_50/cgu_7blk_300chain/1E-6_cycle/"
    # dir_path = "/home/centos/work/1blk_50/cgu_7blk_300chain/5E-6_cycle/"
    # dir_path = "/home/centos/work/1blk_50/cgu_7blk_300chain/1E-5_cycle/"
    # dir_path = "/home/centos/work/1blk_50/cgu_7blk_300chain/250K_5E-6_cycle/"
    # dir_path = "/home/centos/work/1blk_50/cgu_7blk_300chain/5E-6_cycle_changemass/"
    # dir_path = "/home/centos/work/1blk_50/cgu_7blk_300chain/5E-6_back2zero/"
    # dir_path = "/home/centos/work/1blk_50/cgu_7blk_300chain/1E-6_back2zero/"
    # dir_path = "/home/centos/work/1blk_50/cgu_8blk_250chain/5E-6_cycle/"
    # dir_path = "/home/centos/work/1blk_50/cgu_8blk_250chain/change_mass_5E-6/"
    # dir_path = "/home/centos/work/1blk_50/cgu_8blk_250chain/change_mass_1E-6/"
    # dir_path = "/home/centos/work/1blk_50/cgu_8blk_250chain/5E-6_rRESPA/"
    # dir_path = "/home/centos/work/1blk_50/cgu_8blk_250chain/5E-6/"
    dir_path = "/home/centos/work/1blk_50/cgu_8blk_250chain/1E-6/"
    # dir_path = "/home/centos/work/1blk_50/cgu_8blk_250chain/1E-6_2/"
    # dir_path = "/home/centos/work/1blk_50/cgu_8blk_250chain/cycle_loop_1E-6/"
    # dir_path = "/home/centos/work/1blk_50/cgu_8blk_250chain/cycle_loop_1E-6_change_mass/"
    # dir_path = "/home/centos/work/1blk_50/cgu_8blk_250chain/true_strain_1E-6/"
    # dir_path = "/home/centos/work/1blk_50/cgu_8blk_250chain/estrain_1E-6/"
    # dir_path = "/home/centos/work/1blk_50/cgu_8blk_250chain/fix_v_1E-6/"
    # dir_path = "/home/centos/work/1blk_50/cgu_7blk_300chain/1E-6/"
    # dir_path = "/home/centos/work/1blk_50/cgu_8blk_250chain/back2zero_1E-6/"
    # dir_path = "/home/centos/work/1blk_50/cgu_7blk_300chain/change_mass_1E-6/"
    # dir_path = "/home/centos/work/1blk_50/cgu_7blk_300chain/change_mass_1E-6_2/"
    # dir_path = "/home/centos/work/1blk_50/cgu_7blk_300chain/1E-6_no_rmap/"
    # stress_type = False
    # strain_type = Falses

    stress_type = True
    strain_type = True

    window = 100
    # dir_path_list = [dir_path29, dir_path36]
    colors = [cm.tab10(x) for x in np.linspace(0, 1, 10)][:]
    # colors = [cm.rainbow(x) for x in np.linspace(0, 1, 21)][:]
    begin = 1
    end = 1

    for i in range(begin, end + 1):
        fps = dir_path + "strain.s%d.txt" % i
        if outtype == 1:
            columns = ["strain", "stress_x", "stress_y", "stress_z", "l_x", "l_y", "l_z", "temp"]
        elif outtype == 2:
            columns = ["strain", "stress_x", "stress_y", "stress_z", "l_x", "l_y", "l_z", "temp", "virial", "ke"]
        elif outtype == 3:
            columns =  ["strain", "stress_x", "stress_y", "stress_z", "s_virial", "s_ke", "s_pair", "temp"]
            fps = dir_path + "strain.r%d.txt" % i
        else:
            raise KeyError("No such type")
        datas = pd.read_csv(fps, sep=" ", skiprows=1, names=columns, header=None)
        if outtype == 3:
            plot_relax(datas, "r%d" % i, colors[i], window)
        else:
            y0 = datas.l_y[0]
            z0 = datas.l_z[0]
            plot_data(datas, y0, z0, stress_type, strain_type, "s%d" % i, "b", window)

        k = 0
        for j in range(21):
            fpc = dir_path + "strain.b%d.txt" % j
            if os.path.exists(fpc):
                datac = pd.read_csv(fpc, sep=" ", skiprows=1, names=columns, header=None)
                plot_data(datac, y0, z0, stress_type, strain_type, "c%d" % j, colors[k], window, j)
                k+=1
                # plot_data(datac, y0, z0, stress_type, strain_type, "c%d" % i, "r", window)

    # plt.ylim(-0., 8)
    # plt.xlim(0.00, 2)
    # plt.ylim(bottom=0.)
    if strain_type:
        # plt.xscale("log")
        # plt.yscale("log")
        # plt.xlim(0., 1.0)
        # plt.ylim(bottom=0.)
        plt.ylim(0., 100)
        strain_title = "True Strain"
        stress_title = "True Stress"
        # plt.ylim(0., 1000)
    else:
        # plt.xlim(0., 4)
        # plt.ylim(bottom=0.)
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
    # plt.title(dir_path.split("/")[-3] + " " +stress_title + " vs " + strain_title, fontsize=16)
    plt.title(dir_path.split("/")[-3] + " " + dir_path.split("/")[-2], fontsize=16)
    plt.xlabel("strain", fontsize=14)
    plt.ylabel("stress", fontsize=14)
    # plt.title("With dihedral")
    # plt.xlim(0, 1.5)
    plt.legend()
    plt.show()

    # z = np.polyfit(data_s1.strain[190:490], data_s1_x[190:490], 1)
    # print(z)


if __name__ == '__main__':
    main()
