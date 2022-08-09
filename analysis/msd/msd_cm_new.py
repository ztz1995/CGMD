import numpy as np
import pandas as pd
from pyMD import analysis as ans
from collections import OrderedDict
import os
import matplotlib.pyplot as plt
import pickle as pkl
from multiprocessing import Pool


if __name__ == '__main__':

    # # AA
    # title = "aa"
    # steps = [1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000, 2000000, 5000000, 10000000, 20000000, 50000000, 100000000]
    # # nums = 2000
    # nums = 500
    # start_delta = 40000
    # starts = np.arange(nums) * start_delta
    #
    # fp = "/home/centos/fast/aa/1blk_50/long_npt_7/"
    # fp_data = fp + "1blk_50.data"
    # fp_trj = fp + "trj/1blk_50.lammpstrj."
    # fp_info = fp + "cg_struct_info.pkl"
    # atom3d_gt = ans.Atom3dGenCM(fp_data, fp_trj, one_file=False, info_path=fp_info, pad=9)
    #
    # msd_sum = {"cm": np.zeros(len(steps)), "time": np.zeros(len(steps))}
    # ret = list()
    # pool = Pool(24)
    # for start in starts:
    #     ret.append(pool.apply_async(ans.DynamicAnalyzer.mean_squared_displacement, args=(atom3d_gt, start, steps)))
    #     # msd = ans.DynamicAnalyzer.mean_squared_displacement(atom3d_gt, start, steps)
    #     # for k, v in msd.items():
    #     #     msd_sum[k] += v
    # pool.close()
    # pool.join()
    # for r in ret:
    #     for k, v in r.get().items():
    #         msd_sum[k] += v
    # for k in msd_sum:
    #     msd_sum[k] /= nums
    # print(msd_sum)
    # with open("data/msd_cm_%s_%d.pkl" % (title, nums), "wb") as f:
    #     pkl.dump(msd_sum, f)
    # plt.plot(msd_sum["time"], msd_sum["cm"])
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.show()

    # # CG
    # title = "cg"
    # steps = [1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000, 300000, 1000000, 3000000, 10000000, 30000000, 100000000]
    # # nums = 2000
    # nums = 100
    # start_delta = 40000
    # starts = np.arange(nums) * start_delta
    #
    # fp = "/home/centos/work/1blk_160/%s_1blk_160chain/equi/" % title
    # fp_data = fp + "1blk_160.data"
    # fp_trj = fp + "trj/1blk_160.lammpstrj."
    # atom3d_gt = ans.Atom3dGenCM(fp_data, fp_trj, one_file=False, pad=9)
    #
    # msd_sum = {"cm": np.zeros(len(steps)), "time": np.zeros(len(steps))}
    # for start in starts:
    #     print(start)
    #     msd = ans.DynamicAnalyzer.mean_squared_displacement(atom3d_gt, start, steps)
    #     for k, v in msd.items():
    #         msd_sum[k] += v
    # for k in msd_sum:
    #     msd_sum[k] /= nums
    # print(msd_sum)
    # with open("data/msd_cm_%s_%d.pkl" % (title, nums), "wb") as f:
    #     pkl.dump(msd_sum, f)
    # plt.plot(msd_sum["time"], msd_sum["cm"])
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.show()
    #
    # CGU
    title = "cgu"
    # steps = [1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000, 300000, 1000000, 3000000, 10000000, 30000000, 100000000]
    # steps = [100000, 200000, 600000, 1000000, 2000000, 6000000, 10000000, 20000000,60000000, 100000000]
    steps = [1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000, 2000000, 5000000, 10000000, 20000000, 50000000, 100000000]
    # nums = 2000
    nums = 500
    start_delta = 40000
    starts = np.arange(nums) * start_delta

    fp = "/home/centos/work/1blk_160/%s_1blk_160chain/equi/" % title
    fp_data = fp + "1blk_160.data"
    fp_trj = fp + "trj/1blk_160.lammpstrj."
    atom3d_gt = ans.Atom3dGenCM(fp_data, fp_trj, one_file=False, pad=9)

    msd_sum = {"cm": np.zeros(len(steps)), "time": np.zeros(len(steps))}
    ret = list()
    pool = Pool(24)
    for start in starts:
        ret.append(pool.apply_async(ans.DynamicAnalyzer.mean_squared_displacement, args=(atom3d_gt, start, steps)))
        # msd = ans.DynamicAnalyzer.mean_squared_displacement(atom3d_gt, start, steps)
        # for k, v in msd.items():
        #     msd_sum[k] += v
    pool.close()
    pool.join()
    for r in ret:
        for k, v in r.get().items():
            msd_sum[k] += v
    for k in msd_sum:
        msd_sum[k] /= nums
    print(msd_sum)
    with open("data/msd_cm_%s_%d.pkl" % (title, nums), "wb") as f:
        pkl.dump(msd_sum, f)
    plt.plot(msd_sum["time"], msd_sum["cm"])
    plt.xscale('log')
    plt.yscale('log')
    plt.show()