from pyMD import analysis as ans
import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl
import os


def get_h_h_dict(atom3d):
    h_h = dict()
    for atom_id, atom in atom3d.Atoms.items():
        if atom.type == "c1":
            h_list = list()
            connects = atom3d.get_connected_atom_id(atom_id)
            for a_id in connects:
                a = atom3d.Atoms[a_id]
                if a.type == "n1":
                    con = atom3d.get_connected_atom_id(a_id)
                    for n_id in con:
                        nc = atom3d.Atoms[n_id]
                        if nc.type == "h1":
                            h_list.append(n_id)
            h_h[h_list[0]] = h_list[1]
            h_h[h_list[1]] = h_list[0]
    return h_h


if __name__ == '__main__':
    # # AA
    # title = "aa"
    # step = 1
    # frames = 5000
    # h_len = 3.
    # h_ang = 120
    # i = 8
    # # i = 0
    #
    # fp = "/home/centos/fast/aa/1blk_50/short_%d/" % i
    # # fp = "/home/centos/fast/aa/1blk_50/long_8/"
    # fp_data = fp + "1blk_50.data"
    # fp_trj = fp + "trj/1blk_50.lammpstrj."
    #
    # atom3ds = ans.Atom3dGenerator(fp_data, fp_trj, one_file=False, info_path=fp + "cgu_struct_info.pkl", pad=9)
    #
    # auto = ans.AutoCorrelation("continuous", step, frames, atom3ds, h_len=h_len, h_ang=h_ang)
    #
    # r_path = "data/life_h_record_%.1f_%d_500K_%d.pkl" % (h_len, h_ang, i)
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
    # plt.hist(times, bins=100, range=[0, 300])
    # # plt.xlim(0, 300)
    # plt.show()

    # with open("l_%s_%d_%d.pkl" % (title, step, frames), "rb") as f:
    #     times = pkl.load(f)
    #
    # print(sorted(times))
    # plt.hist(times, bins=150)
    # plt.show()

    # CGU
    title = "CGU"
    step = 1

    # frames = 20000
    # h_len = 2.5
    # h_ang = 130
    frames = 100000
    h_len = 3.
    h_ang = 120
    # frames = 20000
    # h_len = 4.
    # h_ang = 90
    i = 1

    fp = "/home/centos/work/1blk_50/cgu_1blk_50chain/short_%d/" % i
    # fp = "/home/centos/fast/aa/1blk_50/long_8/"
    fp_data = fp + "1blk_50.data"
    fp_trj = fp + "short/1blk_50.lammpstrj."

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
    # times = auto.life_time_da(0)
    #
    # h_record = auto.h_record
    # with open(r_path, "wb") as f:
    #     pkl.dump(h_record, f)
    #
    # with open("data/l_da_%s_%d_%d_%.1f_%d_500K_%d.pkl" % (title, step, frames, h_len, h_ang, i), "wb") as f:
    #     pkl.dump(times, f)

    with open("data/l_da_%s_%d_%d_%.1f_%d_500K_%d.pkl" % (title, step, frames, h_len, h_ang, i), "rb") as f:
        times = pkl.load(f)
    all_h = list()
    for k, v in times.items():
        all_h += v
    #     plt.hist(v, bins=2000, range=[0, 2000], density=False)
    #     plt.xlim(0, 300)
    #     plt.title(k.__str__())
    #     plt.show()
    a = plt.hist(all_h, bins=10000, range=[0, 10000], density=True)[0]
    print(a[:20])
    plt.xlim(0, 50)
    plt.ylim(0, 0.04)
    plt.show()


    # from collections import OrderedDict
    # ave_time = dict()
    # for k, v in times.items():
    #     ave_time[k] = np.average(v)
    # ave_time = OrderedDict(sorted(ave_time.items(), key=lambda d: d[1], reverse=True))
    # print(ave_time)


    # h_h_dict = get_h_h_dict(atom3ds())
    # hb1 = list()
    # hb2 = list()
    # covered = set()
    # pairs = list(times.keys())
    # for da in pairs:
    #     if da in covered:
    #         continue
    #     h2 = h_h_dict[da[0]]
    #     if (h2, da[1]) in pairs:
    #         hb2 += times[da]
    #         hb2 += times[(h2, da[1])]
    #         covered.add(da)
    #         covered.add((h2, da[1]))
    #     else:
    #         hb1 += times[da]
    #         covered.add(da)
    # # print(hb1)
    # # print(hb2)
    #
    # hb1_h = plt.hist(hb1, bins=2000, range=[0, 2000], density=True)[0]
    # hb2_h = plt.hist(hb2, bins=2000, range=[0, 2000], density=True, color=[0.81606887,  0.31641358,  0.53254456,  0.44198844])[0]
    # plt.legend(["hb1", "hb2"])
    # plt.xlim(0, 50)
    # plt.show()
    #
    # print(hb1_h)
    # hb1_h = hb1_h * np.arange(0, len(hb1_h))
    # plt.bar(np.arange(0, len(hb1_h))[:300], hb1_h[:300])
    # plt.show()
