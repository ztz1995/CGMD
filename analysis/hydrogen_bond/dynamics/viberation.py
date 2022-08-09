from pyMD import analysis as ans
from pyMD import functions as f
import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl
import os


def cal_period(co_lens, start, steps):
    x = np.array(range(steps))
    y = np.array(co_lens[start: start + steps])

    y -= y.mean()
    ft = np.fft.rfft(y)
    freq = np.fft.rfftfreq(len(y))

    f = freq[abs(ft).argmax()]
    if (1/f) > 700:
        print(freq)
        print(ft)
        plt.plot(x, y)
        plt.show()
        plt.plot(freq, abs(ft))
        plt.show()
    # bs = 0.01 / (3.0E8 / (f * 1E15))
    return 1 / f


def cal_all_period(co_dict, starts, steps):
    h0 = []
    for start in starts:
        for k, v in co_dict.items():
            period = cal_period(v, start, steps)
            h0.append(period)
    return h0


if __name__ == '__main__':

    # CGU
    title = "CGU"
    step = 1
    frames = 2000
    h_len = 3.
    h_ang = 120
    i = "2_change_co"

    fp = "/home/centos/work/1blk_50/cgu_1blk_50chain/short_%s/" % i
    # fp = "/home/centos/fast/aa/1blk_50/long_8/"
    fp_data = fp + "1blk_50.data"
    fp_trj = fp + "short/1blk_50.lammpstrj."

    atom3ds = ans.Atom3dGenerator(fp_data, fp_trj, one_file=False, pad=0)

    uans = ans.UreaAnalyzer(atom3ds)

    co_dict, nh_dict = uans.cal_all_bonds(range(frames), atom3ds)

    r_path = "data/vi_co_%d_%s.pkl" % (frames, i)
    with open(r_path, "wb") as f:
        pkl.dump(co_dict, f)
    r_path = "data/vi_nh_%d_%s.pkl" % (frames, i)
    with open(r_path, "wb") as f:
        pkl.dump(nh_dict, f)

    # r_path = "data/vi_co_%d.pkl" % frames
    # with open(r_path, "rb") as f:
    #     co_dict = pkl.load(f)
    # r_path = "data/vi_nh_%d.pkl" % frames
    # with open(r_path, "rb") as f:
    #     nh_dict = pkl.load(f)

    steps = 800
    starts = range(0, frames-steps, 600)
    cop = cal_all_period(co_dict, starts, steps)
    nhp = cal_all_period(nh_dict, starts, steps)
    # print(cop)
    # print(nhp)
    # plt.hist(cop, range=[15, 21], bins=40)
    # plt.show()
    # plt.hist(nhp)
    # plt.show()
    nhp_c = list()
    for p in nhp:
        if p < 10:
            nhp_c.append(p)
    cop_c = list()
    for p in cop:
        if p < 20:
            cop_c.append(p)
    print(np.array(nhp_c).mean())
    print(np.array(cop_c).mean())

