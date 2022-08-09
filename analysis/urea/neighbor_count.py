from pyMD import analysis as ans
import matplotlib.pyplot as plt
import numpy as np
from pyMD import cartesian_operation as cso
from scipy.spatial import cKDTree
import pickle as pkl
import os


# def cal_neighbor(_atom3d_gen, atom_type="c1", cut_off=5.):
#     u_num_list = np.zeros(7)
#     u_sum = 0
#     for i, atom3d in enumerate(_atom3d_gen):
#         atom3d.convert_to_in_cell()
#         cen_list, extended_list = atom3d.get_central_and_extended_in_range(atom_type, cut_off + 0.5)
#         c_co_list = list(map(lambda atom: atom.coordinate, cen_list))
#         e_co_list = list(map(lambda atom: atom.coordinate, extended_list))
#         # all_list = cen_list + extended_list
#         a_co_list = c_co_list + e_co_list
#         tree = cKDTree(a_co_list)
#         for ca_id, ca in enumerate(cen_list):
#             idx_list = tree.query_ball_point(ca.coordinate, cut_off)
#             num = len(idx_list) - 1
#             if num < 0:
#                 raise IndexError("num must > 1")
#             if num > 6:
#                 continue
#             u_num_list[num] += 1
#             u_sum += 1
#     return u_num_list / u_sum


def cal_neighbor_1(_atom3d_gen, _frames, atom_type="c1", cut_off=5.):
    u_num_list = np.zeros(7)
    u_sum = 0
    for _frames in _frames:
        atom3d = _atom3d_gen(_frames, only_type=[atom_type])
        atom3d.convert_to_in_cell()
        cen_list, extended_list = atom3d.get_central_and_extended_in_range(atom_type, cut_off + 0.5)
        c_co_list = list(map(lambda atom: atom.coordinate, cen_list))
        e_co_list = list(map(lambda atom: atom.coordinate, extended_list))
        # all_list = cen_list + extended_list
        a_co_list = c_co_list + e_co_list
        tree = cKDTree(a_co_list)
        for ca_id, ca in enumerate(cen_list):
            idx_list = tree.query_ball_point(ca.coordinate, cut_off)
            num = len(idx_list) - 1
            if num < 0:
                raise IndexError("num must > 1")
            if num > 6:
                continue
            u_num_list[num] += 1
            u_sum += 1
    return u_num_list / u_sum


if __name__ == '__main__':

    # AA
    # i = 0
    # for i in [1, 2, 3, 5, 7, 8, 9]:
    # for i in [0, 8, 9]:
    #     fp = "/home/centos/model/aa/1blk_50/%d/" % i
    #     # title = "aa"
    #     fp_data = fp + "1blk_50.data"
    #     fp_trj = fp + "1blk_50.lammpstrj"
    #     fp_info = fp + "cgu_struct_info.pkl"
    #     frames = 1000
    #     atom3ds = ans.atom3d_generator_from_aa(fp_info, fp_data, fp_trj, start=15000000, step=1000, frames=frames,
    #                                            one_file=True, fake=True)
    #     nums = cal_neighbor(atom3ds, atom_type="c1", cut_off=5.)
    #     with open("data/aa_u_5.0_%d.pkl" % i, "wb") as f:
    #         pkl.dump(nums, f)

    # plt.bar(range(7), nums)
    # # plt.show()
    #
    # # CGU
    # i = 369
    # fp = "/home/centos/Projects/CGMD/data/20200722_lr=0.05_1blk_50_cgu_ex_hard/iter_%d/cal/" % i  # 353
    # # title = "cgu"
    # fp_data = fp + "1blk_40.data"
    # fp_trj = fp + "300K.lammpstrj."
    # frames = 1000
    # atom3ds = ans.atom3d_generator(fp_data, fp_trj, 0, 500, frames, one_file=False, _type="fg", fake=True)
    # nums = cal_neighbor(atom3ds, atom_type="c1", cut_off=5.)
    # with open("data/cgu_u_5.0_%d.pkl" % i, "wb") as f:
    #     pkl.dump(nums, f)
    # plt.bar(range(8, 15), nums)
    #
    # # CG
    # i = 52
    # fp = "/home/centos/Projects/CGMD/data/20200720_lr=0.05_1blk_50_cg_ex_hard/iter_%d/cal/" % i
    # title = "cgu"
    # fp_data = fp + "1blk_40.data"
    # fp_trj = fp + "300K.lammpstrj."
    # frames = 100
    # atom3ds = ans.atom3d_generator(fp_data, fp_trj, 0, 500, frames, one_file=False, _type="fg", fake=True)
    # nums = cal_neighbor(atom3ds, atom_type="U", cut_off=6.)
    # with open("data/cg_u_5.0_%d.pkl" % i, "wb") as f:
    #     pkl.dump(nums, f)
    # plt.bar(range(16, 23), nums)
    #
    #
    # plt.show()

    # CGU-7blk
    # i = 0
    # fp = "/home/centos/work/1blk_50/cgu_7blk_300chain/init/"
    # i = 1
    # fp = "/home/centos/work/1blk_50/cgu_7blk_300chain/SH6S_anneal_400/"
    # # title = "cgu"
    # fp_data = fp + "strain.data"
    # fp_trj = fp + "record.lammpstrj."
    # frames = 2
    # atom3ds = ans.atom3d_generator(fp_data, fp_trj, 900000, 100000, frames, one_file=False, _type="fg", fake=True)
    # nums = cal_neighbor(atom3ds, atom_type="c1", cut_off=5.)
    # with open("data/cgu_7blk_u_5.0_%d.pkl" % i, "wb") as f:
    #     pkl.dump(nums, f)

    # # CG-7blk
    # i = 0
    # fp = "/home/centos/work/1blk_50/cg_7blk_300chain/init/"
    # # title = "cgu"
    # fp_data = fp + "strain.data"
    # fp_trj = fp + "record.lammpstrj."
    # frames = 20
    # atom3ds = ans.atom3d_generator(fp_data, fp_trj, 8000000, 100000, frames, one_file=False, _type="fg", fake=True)
    # nums = cal_neighbor(atom3ds, atom_type="U", cut_off=5.)
    # with open("data/cg_7blk_u_5.0_%d.pkl" % i, "wb") as f:
    #     pkl.dump(nums, f)

    # CGU new

    # CGU new
    sum_res = np.zeros(7)
    for i in range(8):
        print(i)
        fp = "/home/centos/work/1blk_50/cgu_1blk_50chain/%d/" % i
        # title = "cgu"
        fp_data = fp + "1blk_50.data"
        fp_trj = fp + "trj/1blk_50.lammpstrj."
        frames = 1000
        start = 10000000
        steps = 1000
        atom3ds = ans.Atom3dGenerator(fp_data, fp_trj, one_file=False, pad=9)
        nums = cal_neighbor_1(atom3ds, range(start, start+frames*steps, steps), atom_type="c1", cut_off=5.)
        sum_res += nums
    sum_res /= 8
    with open("data/cgu_u_5.0_new_10ns.pkl", "wb") as f:
        pkl.dump(sum_res, f)
    print(sum_res)

    # # CG new
    # title = "cg"
    # sum_res = np.zeros(7)
    # for i in range(8):
    #     print(i)
    #     fp = "/home/centos/work/1blk_50/cg_1blk_50chain/%d/" % i
    #     fp_data = fp + "1blk_50.data"
    #     fp_trj = fp + "trj/1blk_50.lammpstrj."
    #     frames = 1000
    #     start = 10000000
    #     steps = 1000
    #     atom3ds = ans.Atom3dGenerator(fp_data, fp_trj, one_file=False, pad=9)
    #     nums = cal_neighbor_1(atom3ds, range(start, start + frames * steps, steps), atom_type="U", cut_off=5.)
    #     sum_res += nums
    # sum_res /= 8
    # with open("data/%s_u_5.0_new_10ns.pkl" % title, "wb") as f:
    #     pkl.dump(sum_res, f)
    # print(sum_res)
