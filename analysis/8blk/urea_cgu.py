from pyMD import analysis as ans
import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl
import os

if __name__ == '__main__':

    # CG
    cut_off = 6
    title = "cgu"
    label = "nvt"
    frames = 1000
    steps = 10000
    nums = 500
    start_delta = 10000

    fp = "/home/centos/ztz/1blk_50/%s_8blk_250/equi/" % title
    fp_data = fp + "strain.data"
    fp_trj = fp + "trj/strain.lammpstrj."
    starts = np.arange(nums) * start_delta

    atom3ds = ans.Atom3dGenerator(fp_data, fp_trj, one_file=False, pad=9)
    u_cal = ans.UreaCalculator(atom3ds, title, cut_off)
    auto = ans.AutoCorrelation(u_cal, frames, steps, starts)
    y = auto.cal_parallel(parallel=0)
    print(y)

    with open("data/urea_record_%s_%d_%d_%d_%d.pkl" % (title, steps, frames, nums, start_delta), "wb") as f:
        pkl.dump(u_cal.urecord, f)
    with open("data/urea_%s_%d_%d_%d_%d.pkl" % (title, steps, frames, nums, start_delta), "wb") as f:
        pkl.dump(y, f)

    # with open("data/urea_%s_%d_%d_%d_%d.pkl" % (title, steps, frames, nums, start_delta), "rb") as f:
    #     y = pkl.load(f)
    t = np.arange(frames) * steps / 1000

    plt.style.use('science')
    plt.plot(np.arange(frames) * steps / 1000, y)
    plt.yscale("log")
    plt.yticks([0.2,0.3,0.4,0.6,1.0], [0.2,0.3,0.4,0.6,1.0])
    plt.show()