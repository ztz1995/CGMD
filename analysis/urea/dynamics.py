from pyMD import analysis as ans
import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl
import os

if __name__ == '__main__':

    # # CGU
    # cut_off = 6
    # title = "cgu"
    # label = "npt"
    # frames = 10000
    # steps = 10000
    # nums = 5000
    # start_delta = 40000
    #
    # fp = "/home/centos/work/1blk_50/cgu_1blk_50chain/long_%s/" % label
    # fp_data = fp + "1blk_50.data"
    # fp_trj = fp + "saved/300K.lammpstrj."
    # starts = np.arange(nums) * start_delta
    #
    # atom3ds = ans.Atom3dGenerator(fp_data, fp_trj, one_file=False)
    # u_cal = ans.UreaCalculator(atom3ds, title, cut_off)
    # auto = ans.AutoCorrelation(u_cal, frames, steps, starts)
    # y = auto.cal_parallel(parallel=0)
    # print(y)
    #
    # with open("data/urea_record_%s_%d_%d_%d_%d_1.pkl" % (title, steps, frames, nums, start_delta), "wb") as f:
    #     pkl.dump(u_cal.urecord, f)
    # with open("data/urea_%s_%d_%d_%d_%d_1.pkl" % (title, steps, frames, nums, start_delta), "wb") as f:
    #     pkl.dump(y, f)
    #
    # plt.plot(np.arange(frames) * steps / 1000, y)
    # plt.yscale("log")
    # plt.show()

    # with open("data/urea_%s_%d_%d_%d_%d.pkl" % (title, steps, frames, nums, start_delta), "rb") as f:
    #     y = pkl.load(f)
    # t = np.arange(frames) * steps / 1000
    # plt.plot(t[:100], y[:100])
    # plt.show()


    # # CG
    # cut_off = 6
    # title = "cg"
    # label = "nvt"
    # frames = 10000
    # steps = 10000
    # nums = 4000
    # start_delta = 20000
    #
    # # fp = "/home/centos/work/1blk_50/cg_1blk_50chain/long_%s/" % label
    # # fp_data = fp + "1blk_50.data"
    # # fp_trj = fp + "trj/1blk_50.lammpstrj."
    # # starts = np.arange(nums) * start_delta
    # #
    # # atom3ds = ans.Atom3dGenerator(fp_data, fp_trj, one_file=False)
    # # u_cal = ans.UreaCalculator(atom3ds, title, cut_off)
    # # auto = ans.AutoCorrelation(u_cal, frames, steps, starts)
    # # y = auto.cal_parallel()
    # # print(y)
    # #
    # # with open("data/urea_record_%s_%d_%d_%d_%d.pkl" % (title, steps, frames, nums, start_delta), "wb") as f:
    # #     pkl.dump(u_cal.urecord, f)
    # # with open("data/urea_%s_%d_%d_%d_%d.pkl" % (title, steps, frames, nums, start_delta), "wb") as f:
    # #     pkl.dump(y, f)
    #
    # with open("data/urea_%s_%d_%d_%d_%d.pkl" % (title, steps, frames, nums, start_delta), "rb") as f:
    #     y = pkl.load(f)
    # t = np.arange(frames) * steps / 1000
    #
    # plt.style.use('science')
    # plt.plot(np.arange(frames) * steps / 1000, y)
    # plt.yscale("log")
    # plt.yticks([0.2,0.3,0.4,0.6,1.0], [0.2,0.3,0.4,0.6,1.0])
    # plt.show()

    # # AA
    # cut_off = 5
    # title = "aa"
    # label = "npt_7"
    # frames = 3000
    # steps = 10000
    # nums = 5000
    # start_delta = 20000
    #
    # fp = "/home/centos/fast/aa/1blk_50/long_%s/" % label
    # fp_data = fp + "1blk_50.data"
    # fp_trj = fp + "trj/1blk_50.lammpstrj."
    # starts = np.arange(nums) * start_delta
    #
    # atom3ds = ans.Atom3dGenerator(fp_data, fp_trj, one_file=False, info_path=fp+"cgu_struct_info.pkl", pad=9)
    # u_cal = ans.UreaCalculator(atom3ds, title, cut_off)
    # auto = ans.AutoCorrelation(u_cal, frames, steps, starts)
    # y = auto.cal_parallel(parallel=20, p2=10)
    # print(y)
    #
    # with open("data/urea_record_%s_%d_%d_%d_%d.pkl" % (title, steps, frames, nums, start_delta), "wb") as f:
    #     pkl.dump(u_cal.urecord, f)
    # with open("data/urea_%s_%d_%d_%d_%d.pkl" % (title, steps, frames, nums, start_delta), "wb") as f:
    #     pkl.dump(y, f)
    #
    # t = np.arange(frames) * steps / 1000
    # # plt.style.use('science')
    # plt.plot(np.arange(frames) * steps / 1000, y)
    # plt.yscale("log")
    # plt.show()

    # CGU new
    cut_off = 5
    title = "cgu"
    frames = 4000
    steps = 10000
    nums = 3000
    start_delta = 40000

    fp = "/home/centos/work/1blk_160/%s_1blk_160chain/equi/" % title
    fp_data = fp + "1blk_160.data"
    fp_trj = fp + "trj/1blk_160.lammpstrj."
    starts = np.arange(nums) * start_delta

    atom3ds = ans.Atom3dGenerator(fp_data, fp_trj, one_file=False, pad=9)
    u_cal = ans.UreaCalculator(atom3ds, title, cut_off)
    auto = ans.AutoCorrelation(u_cal, frames, steps, starts)
    y = auto.cal_parallel(parallel=40, p2=10)
    print(y)

    with open("data/urea_record_%s_%d_%d_%d_%d_1.pkl" % (title, steps, frames, nums, start_delta), "wb") as f:
        pkl.dump(u_cal.urecord, f)
    with open("data/urea_%s_%d_%d_%d_%d_1.pkl" % (title, steps, frames, nums, start_delta), "wb") as f:
        pkl.dump(y, f)

    plt.plot(np.arange(frames) * steps / 1000, y)
    plt.yscale("log")
    plt.show()

    # # draw
    # with open("data/urea_%s_%d_%d_%d_%d.pkl" % (title, steps, frames, nums, start_delta), "rb") as f:
    #     y = pkl.load(f)
    # t = np.arange(frames) * steps / 1000
    # plt.plot(t[:100], y[:100])
    # plt.show()

    # # CG new
    # cut_off = 5
    # title = "cg"
    # frames = 4000
    # steps = 10000
    # nums = 2000
    # start_delta = 40000
    #
    # fp = "/home/centos/work/1blk_160/%s_1blk_160chain/equi/" % title
    # fp_data = fp + "1blk_160.data"
    # fp_trj = fp + "trj/1blk_160.lammpstrj."
    # starts = np.arange(nums) * start_delta
    #
    # atom3ds = ans.Atom3dGenerator(fp_data, fp_trj, one_file=False, pad=9)
    # u_cal = ans.UreaCalculator(atom3ds, title, cut_off)
    # auto = ans.AutoCorrelation(u_cal, frames, steps, starts)
    # y = auto.cal_parallel(parallel=20, p2=8)
    # # y = auto.cal()
    #
    # print(y)
    #
    # with open("data/urea_record_%s_%d_%d_%d_%d.pkl" % (title, steps, frames, nums, start_delta), "wb") as f:
    #     pkl.dump(u_cal.urecord, f)
    # with open("data/urea_%s_%d_%d_%d_%d.pkl" % (title, steps, frames, nums, start_delta), "wb") as f:
    #     pkl.dump(y, f)
    #
    # plt.plot(np.arange(frames) * steps / 1000, y)
    # plt.yscale("log")
    # plt.show()

    # draw
    # with open("data/urea_%s_%d_%d_%d_%d.pkl" % (title, steps, frames, nums, start_delta), "rb") as f:
    #     y = pkl.load(f)
    # t = np.arange(frames) * steps / 1000
    # plt.plot(t[:100], y[:100])
    # plt.show()