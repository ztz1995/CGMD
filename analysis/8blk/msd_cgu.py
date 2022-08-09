
import numpy as np
import pandas as pd
from pyMD import analysis as ans
from collections import OrderedDict
import os
import matplotlib.pyplot as plt
import pickle as pkl
from multiprocessing import Pool


if __name__ == '__main__':

    # CGU
    title = "cgu"
    steps = [1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000, 300000, 1000000, 3000000, 10000000]
    # nums = 2000
    nums = 100
    start_delta = 40000
    starts = np.arange(nums) * start_delta

    fp = "/home/centos/ztz/1blk_50/%s_8blk_250/equi/" % title
    fp_data = fp + "strain.data"
    fp_trj = fp + "trj/strain.lammpstrj."
    atom3d_gt = ans.Atom3dGenSegCMV2(fp_data, fp_trj, one_file=False, pad=9)

    msd_sum = {"soft": np.zeros(len(steps)), "hard": np.zeros(len(steps)), "time": np.zeros(len(steps))}

    # for start in starts:
    #     print(start)
    #     msd = ans.DynamicAnalyzer.mean_squared_displacement(atom3d_gt, start, steps)
    #     print(msd)
    #     for k, v in msd.items():
    #         msd_sum[k] += v

    ret = list()
    pool = Pool(24)
    for start in starts:
        ret.append(pool.apply_async(ans.DynamicAnalyzer.mean_squared_displacement, args=(atom3d_gt, start, steps)))
    pool.close()
    pool.join()
    for r in ret:
        for k, v in r.get().items():
            msd_sum[k] += v

    for k in msd_sum:
        msd_sum[k] /= nums
    # print(msd_sum)
    with open("data/msd_cm_%s_%d.pkl" % (title, nums), "wb") as f:
        pkl.dump(msd_sum, f)
    plt.plot(msd_sum["time"], msd_sum["soft"])
    plt.plot(msd_sum["time"], msd_sum["hard"])
    plt.xscale('log')
    plt.yscale('log')
    plt.show()