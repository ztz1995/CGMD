from pyMD import analysis as ans
import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl
import os


if __name__ == '__main__':

    # AA
    fp = "/home/centos/fast/aa/1blk_50/long_8/"
    fp_data = fp + "1blk_50.data"
    fp_trj = fp + "trj/1blk_50.lammpstrj."

    title = "aa"
    frames = 1000
    step = 50
    nums = 100
    h_len = 3.0
    h_ang = 130
    atom3ds = ans.Atom3dGenerator(fp_data, fp_trj, one_file=False, info_path=fp + "cgu_struct_info.pkl", pad=9)
    auto = ans.HBAutoCorrelation("continuous", step, frames, atom3ds, h_len=h_len, h_ang=h_ang)

    record_path = "data/c_h_record_%.1f_%d_500K.pkl" % (h_len, h_ang)
    if os.path.exists(record_path):
        with open(record_path, "rb") as f:
            h_record = pkl.load(f)
        auto.h_record = h_record

    begin = 0
    y = auto.cal_multiple(np.arange(begin, begin + nums) * 5000)

    h_record = auto.h_record
    with open(record_path, "wb") as f:
        pkl.dump(h_record, f)

    with open("data/c_%s_%d_%d_%d_%.1f_%d.pkl" % (title, step, frames, nums, h_len, h_ang), "wb") as f:
        pkl.dump(y, f)
    plt.plot(np.arange(frames) * step / 1000, y)
    plt.yscale("log")
    # plt.xscale("log")
    plt.show()

    # title = "aa"
    # step = 50
    # frames = 1000
    # nums = 1000
    # with open("c_%s_%d_%d_%d.pkl" % (title, step, frames, nums), "rb") as f:
    #     y = pkl.load(f)
    # plt.plot(np.arange(frames) * step / 1000, y)
    # plt.yscale("log")
    # # plt.ylim(1E-4, 1)
    # # plt.xlim(0, 40)
    # plt.title("%s_%d_%d_%d_c" % (title, step, frames, nums))
    # plt.show()