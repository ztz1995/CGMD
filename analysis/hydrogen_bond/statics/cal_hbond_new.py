from pyMD import analysis as ans
import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl
import os
# from draw_hbond import customize_violin

if __name__ == '__main__':
    # # AA
    # title = "aa"
    # hb1 = list()
    # hb2 = list()
    # for j in [0, 1, 2, 3, 5, 7, 8, 9]:
    #     print(j)
    #     fp = "/home/centos/model/aa/1blk_50/%d/" % j
    #     fp_data = fp + "1blk_50.data"
    #     fp_trj = fp + "1blk_50.lammpstrj"
    #     fp_info = fp + "cgu_struct_info.pkl"
    #     frames = 1000
    #     atom3ds = ans.Atom3dGenerator(fp_data, fp_trj, one_file=True, info_path=fp + "cgu_struct_info.pkl", pad=0)
    #
    #     start = 15000000
    #     steps = 1000
    #     auto = ans.HBAutoCorrelation("intermittent", steps, frames, atom3ds, h_len=3.0, h_ang=120)
    #     h_record = auto.create_h_record(range(start, start + frames * steps, steps))
    #
    #     for h_dict in h_record.values():
    #         h1, h2 = auto.h_cal.hb1_hb2(h_dict)
    #         hb1.append(h1)
    #         hb2.append(h2)
    # with open("data/hbond_new_%s.pkl" % title, "wb") as file:
    #     pkl.dump((hb1, hb2), file)
    #
    # axs = plt.subplot()
    # parts = axs.violinplot([hb1, hb2, np.array(hb1)+np.array(hb2)*2], showmeans=True, positions=range(0+1, 9+1, 4))
    # # customize_violin(parts, '#D43F3A', label="CGU-1blk")
    # # customize_violin(parts, '#7c7c7c', label="AA")
    # axs.yaxis.grid(True)
    # plt.show()

    # # CGU
    # title = "cgu"
    # hb1 = list()
    # hb2 = list()
    #
    # # fp = "/home/centos/work/1blk_50/cgu_1blk_50chain/long_npt/"
    # # fp_data = fp + "1blk_50.data"
    # # fp_trj = fp + "saved/300K.lammpstrj."
    # # frames = 10000
    # # start = 10000000
    # # steps = 10000
    # fp = "/home/centos/work/1blk_160/cgu_1blk_160chain/equi/"
    # fp_data = fp + "1blk_160.data"
    # fp_trj = fp + "trj/1blk_160.lammpstrj."
    # frames = 10000
    # start = 10000000
    # steps = 10000
    #
    # atom3ds = ans.Atom3dGenerator(fp_data, fp_trj, one_file=False, pad=9)
    # auto = ans.HBAutoCorrelation("intermittent", steps, frames, atom3ds, h_len=3.0, h_ang=120)
    # h_record = auto.create_h_record(range(start, start + frames * steps, steps))
    #
    # for k, h_dict in h_record.items():
    #     h1, h2, h3 = auto.h_cal.hb1_hb2(h_dict)
    #     if h3 > 0: print(k)
    #     hb1.append(h1)
    #     hb2.append(h2)
    #
    # with open("data/hbond_new_%s_160.pkl" % title, "wb") as file:
    #     pkl.dump((hb1, hb2), file)
    #
    # axs = plt.subplot()
    # parts = axs.violinplot([hb1, hb2, np.array(hb1)+np.array(hb2)*2], showmeans=True, positions=range(0+1, 9+1, 4))
    # # customize_violin(parts, '#D43F3A', label="CGU-1blk")
    # # customize_violin(parts, '#7c7c7c', label="AA")
    # axs.yaxis.grid(True)
    # plt.show()

    # CGU
    title = "cgu"
    hb1 = list()
    hb2 = list()
    frames = 1000
    start = 10000000
    steps = 1000
    for i in range(8):
        fp = "/home/centos/work/1blk_50/cgu_1blk_50chain/%d/" % i
        fp_data = fp + "1blk_50.data"
        fp_trj = fp + "trj/1blk_50.lammpstrj."
        atom3ds = ans.Atom3dGenerator(fp_data, fp_trj, one_file=False, pad=9)
        auto = ans.HBAutoCorrelation("intermittent", steps, frames, atom3ds, h_len=3.0, h_ang=120)
        h_record = auto.create_h_record(range(start, start + frames * steps, steps))

        for k, h_dict in h_record.items():
            h1, h2, h3 = auto.h_cal.hb1_hb2(h_dict)
            if h3 > 0: print(k)
            hb1.append(h1)
            hb2.append(h2)

    with open("data/hbond_new_%s_50.pkl" % title, "wb") as file:
        pkl.dump((hb1, hb2), file)

    axs = plt.subplot()
    parts = axs.violinplot([hb1, hb2, np.array(hb1)+np.array(hb2)*2], showmeans=True, positions=range(0+1, 9+1, 4))
    # customize_violin(parts, '#D43F3A', label="CGU-1blk")
    # customize_violin(parts, '#7c7c7c', label="AA")
    axs.yaxis.grid(True)
    plt.show()
    # less than 0.2%