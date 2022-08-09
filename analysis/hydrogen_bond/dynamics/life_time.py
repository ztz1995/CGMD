from pyMD import analysis as ans
import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl
import os

if __name__ == '__main__':
    # # AA
    # title = "aa"
    # step = 1
    # frames = 100000
    # h_len = 3.
    # h_ang = 120
    # i = 0
    # # i = 0
    #
    # fp = "/home/centos/fast/aa/1blk_50/short_%d/" % i
    # # fp = "/home/centos/fast/aa/1blk_50/long_8/"
    # fp_data = fp + "1blk_50.data"
    # fp_trj = fp + "trj/1blk_50.lammpstrj."
    #
    # atom3ds = ans.Atom3dGenerator(fp_data, fp_trj, one_file=False, info_path=fp + "cgu_struct_info.pkl", pad=9)
    #
    # auto = ans.HBAutoCorrelation("continuous", step, frames, atom3ds, h_len=h_len, h_ang=h_ang)
    #
    # r_path = "data/life_h_record_%.1f_%d_500K_%d.pkl" % (h_len, h_ang, i)
    # # if os.path.exists(r_path):
    # #     with open(r_path, "rb") as f:
    # #         h_record = pkl.load(f)
    # #     auto.h_record = h_record
    #
    # times = auto.life_time(0)
    #
    # h_record = auto.h_record
    # with open(r_path, "wb") as f:
    #     pkl.dump(h_record, f)
    #
    # with open("data/l_%s_%d_%d_%.1f_%d_500K_%d.pkl" % (title, step, frames, h_len, h_ang, i), "wb") as f:
    #     pkl.dump(times, f)
    #
    # print(sorted(times))
    #
    # plt.hist(times, bins=3000, range=[0, 3000])
    # # plt.xlim(0, 300)
    # plt.show()



    # # CGU
    # title = "CGU"
    # step = 1
    # # frames = 100000
    # # h_len = 3.
    # # h_ang = 120
    # frames = 5000
    # h_len = 4.
    # h_ang = 90
    # i = 1
    #
    # fp = "/home/centos/work/1blk_50/cgu_1blk_50chain/short_%d/" % i
    # # fp = "/home/centos/fast/aa/1blk_50/long_8/"
    # fp_data = fp + "1blk_50.data"
    # fp_trj = fp + "short/1blk_50.lammpstrj."
    #
    # atom3ds = ans.Atom3dGenerator(fp_data, fp_trj, one_file=False, pad=0)
    #
    # auto = ans.AutoCorrelation("continuous", step, frames, atom3ds, h_len=h_len, h_ang=h_ang)
    #
    # r_path = "data/life_%s_h_record_%.1f_%d_500K_%d.pkl" % (title, h_len, h_ang, i)
    # if os.path.exists(r_path):
    #     with open(r_path, "rb") as f:
    #         h_record = pkl.load(f)
    #     auto.h_record = h_record
    #
    # times = auto.life_time(0)
    #
    # h_record = auto.h_record
    # with open(r_path, "wb") as f:
    #     pkl.dump(h_record, f)
    #
    # with open("data/l_%s_%d_%d_%.1f_%d_500K_%d.pkl" % (title, step, frames, h_len, h_ang, i), "wb") as f:
    #     pkl.dump(times, f)
    #
    # print(sorted(times))
    #
    # plt.hist(times, bins=3000, range=[0, 3000])
    # plt.xlim(0, 300)
    # plt.show()
    #
    # # h_len = 3.
    # # h_ang = 120
    # # with open("data/l_%s_%d_%d_%.1f_%d_500K_%d.pkl" % ("cgu", 1, 5000, h_len, h_ang, 1), "rb") as f:
    # #     times1 = pkl.load(f)
    # #
    # # with open("data/l_%s_%d_%d_%.1f_%d_500K_%d.pkl" % ("aa", 1, 5000, h_len, h_ang, 8), "rb") as f:
    # #     times2 = pkl.load(f)
    # #
    # # plt.hist(times1, bins=400, range=[0, 2000], density=True)
    # # plt.hist(times2, bins=400, range=[0, 2000], density=True, color=[0.81606887,  0.31641358,  0.53254456,  0.44198844])
    # # plt.legend(["cgu", "aa"])
    # # plt.xlim(0, 300)
    # # plt.show()

    # CGU new
    title = "CGU"
    step = 1
    frames = 100000
    h_len = 3.
    h_ang = 120
    # frames = 50000
    # h_len = 4.
    # h_ang = 90
    # i = "1"
    # i = "origin"
    # i = "change_co"
    # i = "change_nh"
    # i = "change_cnh"
    # i = "change_Phnh"
    # i = "change_nco"
    # i = "change_ocnh"
    # i = "change_improper"
    #
    for i in ["change_cnh", "change_Phnh", "change_nco"]:
        print(i)
        fp = "/home/centos/work/1blk_50/cgu_1blk_50chain/short_%s/" % i
        fp_data = fp + "1blk_50.data"
        fp_trj = fp + "short/1blk_50.lammpstrj."

        atom3ds = ans.Atom3dGenerator(fp_data, fp_trj, one_file=False, pad=0)

        auto = ans.HBAutoCorrelation("continuous", step, frames, atom3ds, h_len=h_len, h_ang=h_ang)

        r_path = "data/life_%s_h_record_%.1f_%d_500K_%s.pkl" % (title, h_len, h_ang, i)
        # if os.path.exists(r_path):
        #     with open(r_path, "rb") as f:
        #         h_record = pkl.load(f)
        # else:
        h_record = auto.create_h_record(range(frames), parallel=20)
        auto.h_record = h_record

        times = auto.life_time(0)

        h_record = auto.h_record
        with open(r_path, "wb") as f:
            pkl.dump(h_record, f)

        with open("data/l_%s_%d_%d_%.1f_%d_500K_%s.pkl" % (title, step, frames, h_len, h_ang, i), "wb") as f:
            pkl.dump(times, f)

        values2, bins = np.histogram(times, bins=1000, range=[0, 2000], density=True)
        x = bins[:-1] + (bins[1] - bins[0]) / 2
        plt.plot(x, values2)
        plt.xlim(0, 50)
        plt.show()


    # # with open("data/l_%s_%d_%d_%.1f_%d_500K_%s.pkl" % (title, step, frames, h_len, h_ang, i), "rb") as f:
    # #     times = pkl.load(f)
    #
    # p = plt.hist(times, bins=3000, range=[0, 3000], density=True)[0]
    # print(p[:20])
    #
    # plt.xlim(0, 50)
    # plt.ylim(0, 0.04)
    # plt.title(i)
    # plt.show()

    # h_len = 3.
    # h_ang = 120
    # with open("data/l_%s_%d_%d_%.1f_%d_500K_%d.pkl" % ("cgu", 1, 5000, h_len, h_ang, 1), "rb") as f:
    #     times1 = pkl.load(f)
    #
    # with open("data/l_%s_%d_%d_%.1f_%d_500K_%d.pkl" % ("aa", 1, 5000, h_len, h_ang, 8), "rb") as f:
    #     times2 = pkl.load(f)
    #
    # plt.hist(times1, bins=400, range=[0, 2000], density=True)
    # plt.hist(times2, bins=400, range=[0, 2000], density=True, color=[0.81606887,  0.31641358,  0.53254456,  0.44198844])
    # plt.legend(["cgu", "aa"])
    # plt.xlim(0, 300)
    # plt.show()