import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from matplotlib import cm


def main():

    y0 = 0.
    z0 = 0.
    plt.figure(figsize=(7,5))
    process = "s1"
    # dir_path = "/home/centos/work/200-125_1E-7_new"
    # dir_path = "/Users/ztz/model/stress/200-125_1E-7_large"
    # dir_path = "/home/centos/work/stress_9blk_2000chain/200-125_1E-7_large"
    # dir_path = "/home/centos/work/stress_7blk_2000chain/200-125_1E-7_large"
    # dir_path = "/home/centos/work/stress_7blk_2000chain/300_5E-7"
    # dir_path = "/home/centos/work/stress_7blk_2000chain/300_1E-7"     # origin
    # dir_path = "/home/centos/work/stress_7blk_2000chain/300_1E-6"
    # dir_path = "/home/centos/work/stress_7blk_2000chain/500_5E-7_change_nonbond"
    # dir_path = "/home/centos/work/stress_7blk_2000chain/500_5E-7_1"
    dir_path9 = "/home/centos/work/stress_7blk_2000chain/500_5E-7_change_bond"
    # dir_path_old1 = "/home/centos/work/7blk_500chain/500_1E-7"
    dir_path_old1 = "/home/centos/work/7blk_2000chain/200-125_1E-7"
    # dir_path = "/home/centos/work/7blk_500chain/500_1E-7_1"
    # dir_path = "/home/centos/work/stress_7blk_500chain/500_1E-6"
    dir_path_old2 = "/home/centos/work/7blk_500chain/500_5E-8"
    # dir_path = "/home/centos/work/stress_7blk_500chain/500_1E-8"
    # dir_path = "/home/centos/work/stress_7blk_500chain/500_5E-7_s1"
    dir_path4 = "/home/centos/work/7blk_500chain/1E-6_fene"  # right
    dir_path3 = "/home/centos/work/7blk_500chain/5E-8_fene"
    dir_path27 = "/home/centos/work/7blk_500chain/1E-4_fene"
    dir_path28 = "/home/centos/work/7blk_500chain/1E-5_fene"
    dir_path1 = "/home/centos/work/5blk_500chain/500_1E-7"
    dir_path2 = "/home/centos/work/pure_soft_5blk_500chain/1E-7new"   # fene, only TO
    # dir_path = "/home/centos/work/7blk_500chain/1E-7_soft"   # fene, try to reduce entanglement, but not soft enough
    dir_path5 = "/home/centos/work/7blk_500chain/5E-8_soft_1"   # fene, try to reduce entanglement, but not soft enough
    dir_path6 = "/home/centos/work/7blk_500chain/1E-6_df"   # diffusion coefficient correction
    dir_path7 = "/home/centos/work/7blk_500chain/1E-7_df"   # diffusion coefficient correction
    dir_path8 = "/home/centos/work/7blk_500chain/5E-7_df"   # diffusion coefficient correction
    dir_path10 = "/home/centos/work/7blk_500chain/2.5E-7_df"   # diffusion coefficient correction
    dir_path11 = "/home/centos/work/7blk_500chain/1E-5_df"   # diffusion coefficient correction
    dir_path12 = "/home/centos/work/7blk_500chain/1E-4_df"   # diffusion coefficient correction
    dir_path13 = "/home/centos/work/7blk_500chain/5E-5_df"   # diffusion coefficient correction
    dir_path14 = "/home/centos/work/7blk_500chain/1E-4_fene"   # fene
    dir_path15 = "/home/centos/work/pure_soft_7blk_500chain/1E-4_TO"   # only TO
    dir_path16 = "/home/centos/work/pure_soft_7blk_500chain/5E-5_TO"   # only TO
    dir_path19 = "/home/centos/work/pure_soft_7blk_500chain/1E-5_TO"   # only TO
    dir_path20 = "/home/centos/work/pure_soft_7blk_500chain/1E-6_TO"   # only TO
    dir_path17 = "/home/centos/work/7blk_500chain/1E-5_df_100"   # hard_bead, mass*100, diffusion coefficient correction
    dir_path18 = "/home/centos/work/cross_7blk_500chain/1E-6_3_df"   # cross-link, B-B, diffusion coefficient correction
    dir_path21 = "/home/centos/work/cross_7blk_500chain/1E-6_2.4_df"   # cross-link, B-B, diffusion coefficient correction
    dir_path24 = "/home/centos/work/cross_7blk_500chain/1E-5_2.4_fene"   # cross-link, B-B, diffusion coefficient correction
    dir_path22 = "/home/centos/work/cross_7blk_500chain/1E-6_2.4_fene"   # cross-link, B-B, diffusion coefficient correction
    dir_path23 = "/home/centos/work/cross_7blk_500chain/1E-7_2.4_fene"   # cross-link, B-B, diffusion coefficient correction
    dir_path25 = "/home/centos/work/cross_7blk_500chain/1E-4_2.4_fene"   # cross-link, B-B, diffusion coefficient correction
    dir_path26 = "/home/centos/work/cross_7blk_500chain/5E-8_2.4_fene"   # cross-link, B-B, diffusion coefficient correction
    dir_path29 = "/home/centos/work/fg_7blk_500chain/1E-6"   # fine grained, no cross-link
    dir_path32 = "/home/centos/work/fg_7blk_500chain/1E-6_1"   # fine grained, no cross-link
    dir_path34 = "/home/centos/work/fg_7blk_500chain/1E-6_2"   # fine grained, no cross-link
    dir_path30 = "/home/centos/work/fg_7blk_500chain/1E-7"   # fine grained, no cross-link
    dir_path31 = "/home/centos/work/fg_7blk_500chain/5E-7"   # fine grained, no cross-link
    dir_path33 = "/home/centos/work/fg_7blk_500chain/1E-5"   # fine grained, no cross-link
    dir_path35 = "/home/centos/work/fgh_7blk_200chain/init"   # fine grained, no cross-link
    dir_path36 = "/home/centos/work/fgh_7blk_200chain/2E-7"   # fine grained, no cross-link

    window = 100
    stress_type = False
    strain_type = False
    # stress_type = True
    # strain_type = True
    # dir_path_list = [dir_path1, dir_path2]
    # dir_path_list = [dir_path2]
    # dir_path_list = [dir_path3]
    # dir_path_list = [dir_path3,dir_path7, dir_path_old1, dir_path6]
    # dir_path_list = [dir_path6, dir_path3]
    # dir_path_list = [dir_path12, dir_path13,dir_path11, dir_path6, dir_path8, dir_path10, dir_path7]
    # dir_path_list = [dir_path6, dir_path8, dir_path10, dir_path7]
    # dir_path_list = [dir_path12, dir_path14]

    # dir_path_list = [dir_path17, dir_path11]    # change mass
    # dir_path_list = [dir_path27, dir_path28,dir_path4, dir_path3]    # no df, different rate
    # dir_path_list = [dir_path18, dir_path21, dir_path6]    # different cross link
    # dir_path_list = [dir_path12, dir_path15, dir_path13, dir_path16, dir_path11, dir_path19, dir_path6, dir_path20]  # df and pure soft
    # dir_path_list = [dir_path6, dir_path20, dir_path7, dir_path2, dir_path_old2]  # df and pure soft
    # dir_path_list = [dir_path12, dir_path27, dir_path11, dir_path28,dir_path6, dir_path4]  # df and no df
    # dir_path_list = [dir_path25, dir_path24, dir_path22, dir_path23]  # cross-link, different rate, no df
    # dir_path_list = [dir_path25, dir_path27, dir_path24,dir_path28, dir_path22, dir_path4, dir_path26, dir_path3]  # cross-link and no cross
    # dir_path_list = [dir_path22, dir_path4, dir_path26, dir_path3]  # cross-link and no cross
    # dir_path_list = [dir_path4, dir_path29, dir_path30]
    # dir_path_list = [dir_path33,dir_path28, dir_path29, dir_path31, dir_path30, dir_path_old1]
    # dir_path_list = [dir_path33,dir_path29, dir_path31, dir_path30]

    dir_path_list = [dir_path8, dir_path36]
    # dir_path_list = [dir_path29, dir_path36]
    colors = [cm.tab10(x) for x in np.linspace(0, 1, 10)][:]
    # print(colors)
    for num, dir_path in enumerate(dir_path_list):
        csv_lines = list()
        file_path_1 = dir_path + "/strain.%s_1.txt" % process
        if os.path.exists(file_path_1):
            with open(file_path_1, "r") as f:
                lines_1 = f.readlines()
                csv_lines += lines_1[1:]
        file_path_1 = dir_path + "/strain.%s_2.txt" % process
        if os.path.exists(file_path_1):
            with open(file_path_1, "r") as f:
                lines_1 = f.readlines()
                csv_lines += lines_1[1:]
        file_path_1 = dir_path + "/strain.%s.txt" % process
        if os.path.exists(file_path_1):
            with open(file_path_1, "r") as f:
                lines_1 = f.readlines()
                csv_lines += lines_1[1:]

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
            stress = data.stress_x.rolling(window, min_periods=window, center=True).mean() * y * z

        strain = data.strain

        # label = "block-copolymer" if num ==0 else "pure soft"
        label = dir_path.split("/")[-1]
        if strain_type:
            plt.plot(np.log(1 + strain), stress*1000, label=label, color=colors[num]) # , label=process
        else:
            plt.plot(strain, stress*1000, label=label, color=colors[num]) # , label=process

    # plt.ylim(-0., 8)
    # plt.xlim(0.00, 2)
    # plt.ylim(bottom=0.)
    if strain_type:
        plt.xlim(0., 1.6)
        plt.ylim(bottom=0.)
        # plt.ylim(0., 200)
        strain_title = "True Strain"
        stress_title = "True Stress"
        # plt.ylim(0., 1000)
    else:
        plt.xlim(0., 4)
        plt.ylim(bottom=0.)
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
    plt.title(stress_title + " vs " + strain_title)
    # plt.title("With dihedral")
    plt.legend()
    plt.show()

    # z = np.polyfit(data_s1.strain[190:490], data_s1_x[190:490], 1)
    # print(z)


if __name__ == '__main__':
    main()


