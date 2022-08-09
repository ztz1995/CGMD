from pyMD import analysis as ans
import matplotlib.pyplot as plt
import numpy as np
from pyMD import cartesian_operation as cso
from scipy.spatial import cKDTree
import pickle as pkl
import os
import random


class U:
    def __init__(self, c3_id):
        self.o1 = 0
        self.c3 = c3_id
        self.h1 = list()
        self.id_co = dict()

    @property
    def c3_co(self):
        return self.id_co[self.c3]

    def get_copy(self):
        nu = U(self.c3)
        nu.o1 = self.o1
        nu.h1 = self.h1[:]
        for k, v in self.id_co.items():
            nu.id_co[k] = v.copy()
        return nu


def get_u_list(a3d):
    _u_list = list()
    for atom_id, atom in a3d.Atoms.items():
        if atom.type == "c1":
            u = U(atom_id)
            u.id_co[atom_id] = atom.coordinate.copy() 
            connects = a3d.get_connected_atom_id(atom_id)
            for a_id in connects:
                a = a3d.Atoms[a_id]
                if a.type == "o1":
                    u.o1 = a_id
                    u.id_co[a_id] = a.coordinate.copy()
                elif a.type == "n1":
                    con = a3d.get_connected_atom_id(a_id)
                    for n_id in con:
                        nc = a3d.Atoms[n_id]
                        if nc.type == "h1":
                            u.h1.append(n_id)
                            u.id_co[n_id] = nc.coordinate.copy()
            _u_list.append(u)
    return _u_list


def get_extended_u_in_range(a3d, c_u_list, extend_range):
    def get_moved_u(_u, x, y, z):
        new_u = _u.get_copy()
        d = np.asarray([x, y, z], dtype=int)
        for k in new_u.id_co.keys():
            new_u.id_co[k] -= d * a3d.lattice_parameter
        return new_u

    extended_list = list()
    upper = a3d.lattice_parameter + a3d.box_h - extend_range
    lower = a3d.box_l + extend_range
    for u in c_u_list:
        record = (u.c3_co > upper).astype(int) + (u.c3_co <= lower).astype(int) * -1
        if record[0] != 0:
            extended_list.append(get_moved_u(u, record[0], 0, 0))
            if record[1] != 0:
                extended_list.append(get_moved_u(u, record[0], record[1], 0))
                if record[2] != 0:
                    extended_list.append(get_moved_u(u, record[0], record[1], record[2]))
            if record[2] != 0:
                extended_list.append(get_moved_u(u, record[0], 0, record[2]))
        if record[1] != 0:
            extended_list.append(get_moved_u(u, 0, record[1], 0))
            if record[2] != 0:
                extended_list.append(get_moved_u(u, 0, record[1], record[2]))
        if record[2] != 0:
            extended_list.append(get_moved_u(u, 0, 0, record[2]))
    return extended_list


def cal_h_bond(u1, u2, h_cut=2.5):
    h = 0
    for h_id in u1.h1:
        d = cso.calculate_distance(u1.id_co[h_id], u2.id_co[u2.o1])
        if d <= h_cut:
            h += 1
    for h_id in u2.h1:
        d = cso.calculate_distance(u2.id_co[h_id], u1.id_co[u1.o1])
        if d <= h_cut:
            h += 1
    return h


def cal_h_bond_frames(_atom3d_gen, h_cut=3., max_u=None):
    h1_list = list()
    h2_list = list()
    h_sum_list = list()
    # strain = list()
    for i, atom3d in enumerate(_atom3d_gen):
        atom3d.convert_to_in_cell()
        cen_u_list = get_u_list(atom3d)
        # print(cen_u_list)
        ext_u_list = get_extended_u_in_range(atom3d, cen_u_list, 10.)
        # for u in ext_u_list:
        #     print(u.id_co)
        random.shuffle(cen_u_list)
        random.shuffle(ext_u_list)
        cu_co_list = list(map(lambda u: u.c3_co, cen_u_list))
        eu_co_list = list(map(lambda u: u.c3_co, ext_u_list))
        all_u_list = cen_u_list + ext_u_list
        au_co_list = cu_co_list + eu_co_list
        tree = cKDTree(au_co_list)
        cut_off = 8.
        h1 = 0
        h2 = 0
        for cu_id, cu in enumerate(cen_u_list):
            if max_u is not None:
                if cu_id >= max_u:
                    break
            idx_list = tree.query_ball_point(cu.c3_co, cut_off)
            for idx in idx_list:
                if idx != cu_id:
                    h_num = cal_h_bond(cu, all_u_list[idx], h_cut=h_cut)
                    if h_num == 1:
                        h1 += 1
                    # elif h_num == 2:
                    #     h2 += 1
                    # elif h_num > 2:
                    #     print(i)
                    #     print(cu.c3)
                    #     print(all_u_list[idx].c3)
                    #     raise ValueError("More than 2 h-bond!")
                    elif h_num > 1:
                        h2 += 1
            # print(h1, h2)
        # strain.append(i * 2E-7)
        h1_list.append(h1 / 2)
        h2_list.append(h2 / 2)
        h_sum_list.append(h1 / 2 + h2)

    # plt.plot(strain, h1_list, label="h1")
    # plt.plot(strain, h2_list, label="h2")
    # plt.plot(strain, h_sum_list, label="h_sum")
    # plt.legend()
    # plt.show()
    return h1_list, h2_list, h_sum_list


if __name__ == '__main__':
    # # i = 0
    # for i in [1, 2, 3, 5, 7, 8, 9]:
    #     fp = "/home/centos/model/aa/1blk_50/%d/" % i
    #     title = "aa"
    #     fp_data = fp + "1blk_50.data"
    #     fp_trj = fp + "1blk_50.lammpstrj"
    #     fp_info = fp + "cgu_struct_info.pkl"
    #     frames = 1000
    #     atom3ds = ans.atom3d_generator_from_aa(fp_info, fp_data, fp_trj, start=15000000, step=1000, frames=frames,
    #                                            one_file=True, fake=False)

    # i = 369
    # fp = "/home/centos/Projects/CGMD/data/20200722_lr=0.05_1blk_50_cgu_ex_hard/iter_%d/cal/" % i  # 353
    # title = "cgu"
    # fp_data = fp + "1blk_40.data"
    # fp_trj = fp + "300K.lammpstrj."
    # frames = 1000
    # atom3ds = ans.atom3d_generator(fp_data, fp_trj, 0, 500, frames, one_file=False, _type="fg", fake=False)

    i = 0
    fp = "/home/centos/work/1blk_50/cgu_7blk_300chain/init/"
    title = "cgu_7blk"
    fp_data = fp + "strain.data"
    fp_trj = fp + "record.lammpstrj."
    frames = 15
    chains = 300
    blocks = 7
    max_u_num = 100
    atom3ds = ans.atom3d_generator(fp_data, fp_trj, 5000000, 100000, frames, one_file=False, _type="fg", fake=False)

    h_cut = 2.5
    h1, h2, hs = cal_h_bond_frames(atom3ds, h_cut=h_cut, max_u=max_u_num)
    # h1 = np.asarray(h1) / chains / blocks / 2
    # h2 = np.asarray(h2) / chains / blocks / 2
    # hs = np.asarray(hs) / chains / blocks / 2
    h1 = np.asarray(h1) / max_u_num
    h2 = np.asarray(h2) / max_u_num
    hs = np.asarray(hs) / max_u_num
    h_dict = {"h1": h1, "h2": h2, "hs": hs, "h_cut": h_cut, "frames": frames, "i": i}
    os.makedirs("data", exist_ok=True)
    with open("data/%s_cut_%.1f_%d.pkl" % (title, h_cut, i), "wb") as f:
        pkl.dump(h_dict, f)

    axs = plt.subplot()
    axs.violinplot([h1, h2, hs], showmeans=True)
    axs.yaxis.grid(True)
    axs.set_xticks([y + 1 for y in range(3)])

    plt.setp(axs, xticks=[y + 1 for y in range(3)], xticklabels=['h1', 'h2', 'hs'])
    plt.title("%s cut_off=%.1f frames=%d" % (title, h_cut, frames))
    plt.show()
