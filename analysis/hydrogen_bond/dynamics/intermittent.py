from pyMD import analysis as ans
import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl
import os

if __name__ == '__main__':
    # # i = 369
    # # fp = "/home/centos/Projects/CGMD/data/20200722_lr=0.05_1blk_50_cgu_ex_hard/iter_%d/cal/" % i  # 369
    # # fp = "/home/centos/Projects/CGMD/data/20200826_lr=0.05_1blk_50_cgu_ex_hard_change_mass/iter_12/cal/"

    # # CGU
    # h_len = 3.
    # h_ang = 120
    # title = "cgu"
    # label = "npt"
    # frames = 4000
    # steps = 10000
    # nums = 2000
    # start_delta = 80000
    #
    # fp = "/home/centos/work/1blk_160/cgu_1blk_160chain/equi/"
    # fp_data = fp + "1blk_160.data"
    # fp_trj = fp + "trj/1blk_160.lammpstrj."
    # starts = np.arange(nums) * start_delta
    #
    # atom3ds = ans.Atom3dGenerator(fp_data, fp_trj, one_file=False, pad=9)
    # auto = ans.HBAutoCorrelation("intermittent", steps, frames, atom3ds, h_len=h_len, h_ang=h_ang)
    #
    # record_path = "data/i_h_record_%.1f_%d_%s_%s_160chain.pkl" % (h_len, h_ang, title, label)
    #
    # frame_set = set()
    # for start in starts:
    #     frame_set.update(set(np.arange(frames) * steps + start))
    #
    # if os.path.exists(record_path):
    #     with open(record_path, "rb") as f:
    #         h_record = pkl.load(f)
    #     frame_set -= set(h_record.keys())
    #     new_h_record = auto.create_h_record(frame_set)
    #     h_record.update(new_h_record)
    # else:
    #     h_record = auto.create_h_record(frame_set)
    # auto.h_record = h_record
    # with open(record_path, "wb") as f:
    #     pkl.dump(h_record, f)
    #
    # y = auto.cal_multiple_parallel(starts, parallel=10)
    # h_record = auto.h_record
    #
    # with open("data/i_%.1f_%d_%s_%d_%d_%d_%d_160chain.pkl" % (h_len, h_ang, title, steps, frames, nums, start_delta),
    #           "wb") as f:
    #     pkl.dump(y, f)
    #
    # plt.plot(np.arange(frames) * steps / 1000, y)
    # plt.yscale("log")
    # plt.show()

    # AA
    h_len = 3.
    h_ang = 120
    title = "aa"
    label = "npt"
    frames = 6000
    steps = 10000
    nums = 2000
    start_delta = 40000

    fp = "/home/centos/fast/aa/1blk_50/long_npt_7/"
    fp_data = fp + "1blk_50.data"
    fp_trj = fp + "trj/1blk_50.lammpstrj."
    starts = np.arange(nums) * start_delta

    atom3ds = ans.Atom3dGenerator(fp_data, fp_trj, one_file=False, pad=9, info_path=fp+"cgu_struct_info.pkl")
    auto = ans.HBAutoCorrelation("intermittent", steps, frames, atom3ds, h_len=h_len, h_ang=h_ang)

    record_path = "data/i_h_record_%.1f_%d_%s_%s_50chain.pkl" % (h_len, h_ang, title, label)

    frame_set = set()
    for start in starts:
        frame_set.update(set(np.arange(frames) * steps + start))

    # if os.path.exists(record_path):
    #     with open(record_path, "rb") as f:
    #         h_record = pkl.load(f)
    #     frame_set -= set(h_record.keys())
    #     new_h_record = auto.create_h_record(frame_set)
    #     h_record.update(new_h_record)
    # else:
    h_record = auto.create_h_record(frame_set)
    auto.h_record = h_record
    with open(record_path, "wb") as f:
        pkl.dump(h_record, f)

    y = auto.cal_multiple_parallel(starts, parallel=10)
    h_record = auto.h_record

    with open("data/i_%.1f_%d_%s_%d_%d_%d_%d_50chain.pkl" % (h_len, h_ang, title, steps, frames, nums, start_delta),
              "wb") as f:
        pkl.dump(y, f)

    plt.plot(np.arange(frames) * steps / 1000, y)
    plt.yscale("log")
    plt.show()
