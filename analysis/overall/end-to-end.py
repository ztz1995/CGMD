from pyMD import analysis as ans
import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl
import os
from multiprocessing import Pool


def cal_end_end(starts, _i):
    fp = "/home/centos/model/aa/1blk_50/%d/" % _i
    fp_data = fp + "1blk_50.data"
    fp_trj = fp + "1blk_50.lammpstrj"
    atom3ds = ans.Atom3dGenerator(fp_data, fp_trj, one_file=True, info_path=fp + "cg_struct_info.pkl")

    eed = ans.EndToEndDistanceType(atom3ds)
    dists = eed.cal_dists(starts)
    return dists


if __name__ == '__main__':

    # frames = 1000
    frames = 500
    steps = 10000
    # start = 0
    start = 5000000

    fp = "/home/centos/work/1blk_50/EsTO13Es_1000chain/init/"
    fp_data = fp + "1blk_1000.data"
    fp_trj = fp + "trj/1blk_1000.lammpstrj."
    atom3ds = ans.Atom3dGenerator(fp_data, fp_trj, one_file=False, pad=9)

    ns = range(14, 15)
    eed = ans.EndToEndDistance(atom3ds, chain_len=15, chain_num=1000, ns=ns)
    dists = eed.cal_dists(range(start, start + frames * steps, steps))
    for n in ns:
        # with open("data/EndToEnd_n=%d_%d_%d" % (n, frames, steps), "wb") as file:
        #     pkl.dump(dists[n], file)

        plt.hist(dists[n], bins=50, range=[0, max(dists[n])], density=True)
        plt.title("EndToEnd_n=%d" % n)
        plt.show()

    # # AA
    # title = "aa"
    # frames = 10
    # steps = 1000
    # start = 15000000
    # dist_all = list()
    # pool = Pool(8)
    # ret = list()
    # for i in [0, 1, 2, 3, 5, 7, 8, 9]:
    #     ret.append(pool.apply_async(cal_end_end, args=(i,)))
    #
    # pool.close()
    # pool.join()
    # for r in ret:
    #     dists = r.get()
    #     dist_all += dists
    #
    # with open("data/EndToEnd_%s" % title, "wb") as file:
    #     pkl.dump(dist_all, file)
    #
    # plt.hist(dist_all, bins=50, range=[0, max(dist_all)], density=True)
    # plt.title("EndToEnd_%s" % title)
    # plt.show()

    # # CGU
    # title = "cgu"
    # frames = 10000
    # steps = 10000
    # start = 0
    # starts = range(start, start + frames * steps, steps)
    #
    # fp = "/home/centos/work/1blk_160/cgu_1blk_160chain/equi/"
    # fp_data = fp + "1blk_160.data"
    # fp_trj = fp + "trj/1blk_160.lammpstrj."
    # atom3ds = ans.Atom3dGenerator(fp_data, fp_trj, one_file=False, pad=9)
    #
    # eed = ans.EndToEndDistanceType(atom3ds)
    # dist_all = eed.cal_dists_parallel(starts)
    #
    # with open("data/EndToEnd_%s" % title, "wb") as file:
    #     pkl.dump(dist_all, file)
    #
    # plt.hist(dist_all, bins=50, range=[0, max(dist_all)], density=True)
    # plt.title("EndToEnd_%s" % title)
    # plt.show()

    # # CG
    # title = "cg"
    # frames = 10000
    # steps = 10000
    # start = 0
    # starts = range(start, start + frames * steps, steps)
    #
    # fp = "/home/centos/work/1blk_160/cg_1blk_160chain/equi/"
    # fp_data = fp + "1blk_160.data"
    # fp_trj = fp + "trj/1blk_160.lammpstrj."
    # atom3ds = ans.Atom3dGenerator(fp_data, fp_trj, one_file=False, pad=9)
    #
    # eed = ans.EndToEndDistanceType(atom3ds)
    # dist_all = eed.cal_dists_parallel(starts)
    #
    # with open("data/EndToEnd_%s" % title, "wb") as file:
    #     pkl.dump(dist_all, file)
    #
    # plt.hist(dist_all, bins=50, range=[0, max(dist_all)], density=True)
    # plt.title("EndToEnd_%s" % title)
    # plt.show()
