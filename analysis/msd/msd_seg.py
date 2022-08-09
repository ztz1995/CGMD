import numpy as np
import pandas as pd
from pyMD import analysis as ans
from collections import OrderedDict
import os

if __name__ == '__main__':

    # fp = "/home/centos/model/aa/1blk_40/long/1/"
    # fp_data = fp + "1blk_40.data"
    # fp_trj = fp + "trj/1blk_40.lammpstrj."
    # fp_info = fp + "cg_struct_info.pkl"
    # msd_name = "aa_seg_msd.csv"
    #
    # atom3d_gt = ans.Atom3dGenSegCM(fp_data, fp_trj, one_file=False, info_path=fp_info)
    # a3d = atom3d_gt()
    #
    # # print(atom3d_gt.seg_dict)
    # # a3d = atom3d_gt(0)
    # # for k, v in atom3d_gt.seg_dict.items():
    # #     print(k)
    # #     for a_id in v:
    # #         print(atom3d_gt.atom3d.Atoms[a_id], end=" ")
    # #         print()
    #
    # if os.path.exists(msd_name):
    #     df = pd.read_csv(msd_name)
    # else:
    #     data_dict = OrderedDict()
    #     data_dict["time"] = [0]
    #     for at in a3d.Atom_type_dict:
    #         data_dict[at] = [0.]
    #     df = pd.DataFrame(data_dict)
    #
    # # for t_p in [1, 10, 100, 1000, 10000, 100000]:
    # # for t_p in [1, 10, 100, 1000, 10000, 100000, 1000000, 10000000]:
    # # for t_p in [1, 10, 100, 1000, 10000, 100000, 1000000]:
    # # for t_p in [1, 10, 100, 1000, 1000 0, 100000]:
    # for t_p in [1000000, 10000000]:
    #     if t_p == 1:
    #         times = np.arange(1, 11) * t_p
    #     elif t_p == 10000000:
    #         times = np.arange(2, 10) * t_p
    #     else:
    #         times = np.arange(2, 11) * t_p
    #     starts = np.arange(50) * t_p
    #     max_step = 100 * t_p
    #     if t_p >= 1000000:
    #         starts = np.arange(100) * t_p / 10
    #         max_step = 10 * t_p
    #     for t in times:
    #         msd = ans.DynamicAnalyzer.mean_squared_displacement(atom3d_gt, starts, t, max_step=max_step)
    #         print(msd)
    #         df.loc[df.shape[0]] = msd
    # # df.sort_values("time", inplace=True)
    # print(df)
    # df.to_csv(msd_name, index=False)

    # # aa_cgu
    # fp = "/home/centos/model/aa/1blk_40/long/1/"
    # fp_data = fp + "1blk_40.data"
    # fp_trj = fp + "trj/1blk_40.lammpstrj."
    # fp_info = fp + "cgu_struct_info.pkl"
    # msd_name = "aa_cgu_msd.csv"
    #
    # atom3d_gt = ans.Atom3dGenerator(fp_data, fp_trj, one_file=False, info_path=fp_info)
    # a3d = atom3d_gt()
    # if os.path.exists(msd_name):
    #     df = pd.read_csv(msd_name)
    # else:
    #     data_dict = OrderedDict()
    #     data_dict["time"] = [0]
    #     for at in a3d.Atom_type_dict:
    #         data_dict[at] = [0.]
    #     df = pd.DataFrame(data_dict)
    #
    # for t_p in [1, 10, 100, 1000, 10000, 100000]:
    # # for t_p in [1000000]:
    #     if t_p == 1:
    #         times = np.arange(1, 11) * t_p
    #     else:
    #         times = np.arange(2, 11) * t_p
    #     starts = np.arange(100) * t_p
    #     for t in times:
    #         msd = ans.DynamicAnalyzer.mean_squared_displacement(atom3d_gt, starts, t, max_step=100 * t_p)
    #         print(msd)
    #         df.loc[df.shape[0]] = msd
    # # df.sort_values("time", inplace=True)
    # print(df)
    # df.to_csv(msd_name, index=False)

    # cgu
    fp = "/home/centos/Projects/CGMD/data/20200506_1blk_40_cgu/iter_02/long/"
    # fp = "/home/centos/Projects/CGMD/data/20200506_1blk_40_cgu/iter_02/long2/"
    fp_data = fp + "1blk_40.data"
    fp_trj = fp + "trj/cgu_1blk_40.lammpstrj."
    # msd_name = "cgu_msd.csv"
    msd_name = "cgu_seg_msd.csv"

    # atom3d_gt = ans.Atom3dGenerator(fp_data, fp_trj, one_file=False)
    atom3d_gt = ans.Atom3dGenSegCM(fp_data, fp_trj, one_file=False)
    a3d = atom3d_gt()
    if os.path.exists(msd_name):
        df = pd.read_csv(msd_name)
    else:
        data_dict = OrderedDict()
        data_dict["time"] = [0]
        for at in a3d.Atom_type_dict:
            data_dict[at] = [0.]
        df = pd.DataFrame(data_dict)

    # for t_p in [1, 10, 100, 1000, 10000, 100000]:
    for t_p in [1, 10, 100, 1000, 10000, 100000, 1000000, 10000000]:
    # for t_p in [1, 10, 100, 1000, 10000, 100000, 1000000]:
        if t_p == 1:
            times = np.arange(1, 11) * t_p
        elif t_p == 10000000:
            times = np.arange(2, 10) * t_p
        else:
            times = np.arange(2, 11) * t_p
        starts = np.arange(50) * t_p
        max_step = 100 * t_p
        if t_p == 10000000:
            starts = np.arange(100) * t_p / 10
            max_step = 10 * t_p
        for t in times:
            msd = ans.DynamicAnalyzer.mean_squared_displacement(atom3d_gt, starts, t, max_step=max_step)
            print(msd)
            df.loc[df.shape[0]] = msd
    # df.sort_values("time", inplace=True)
    print(df)
    df.to_csv(msd_name, index=False)

    # cg
    fp = "/home/centos/Projects/CGMD/data/20200506_1blk_40_cg/iter_02/long/"
    fp_data = fp + "1blk_40.data"
    fp_trj = fp + "trj/cgu_1blk_40.lammpstrj."
    # msd_name = "cg_msd.csv"
    msd_name = "cg_seg_msd.csv"

    # atom3d_gt = ans.Atom3dGenerator(fp_data, fp_trj, one_file=False)
    atom3d_gt = ans.Atom3dGenSegCM(fp_data, fp_trj, one_file=False)
    a3d = atom3d_gt()
    if os.path.exists(msd_name):
        df = pd.read_csv(msd_name)
    else:
        data_dict = OrderedDict()
        data_dict["time"] = [0]
        for at in a3d.Atom_type_dict:
            data_dict[at] = [0.]
        df = pd.DataFrame(data_dict)

    # for t_p in [1, 10, 100, 1000, 10000, 100000]:
    for t_p in [1, 10, 100, 1000, 10000, 100000, 1000000, 10000000]:
        if t_p == 1:
            times = np.arange(1, 11) * t_p
        elif t_p == 10000000:
            times = np.arange(2, 10) * t_p
        else:
            times = np.arange(2, 11) * t_p
        starts = np.arange(50) * t_p
        max_step = 100 * t_p
        if t_p == 10000000:
            starts = np.arange(100) * t_p / 10
            max_step = 10 * t_p
        for t in times:
            msd = ans.DynamicAnalyzer.mean_squared_displacement(atom3d_gt, starts, t, max_step=max_step)
            print(msd)
            df.loc[df.shape[0]] = msd
    # df.sort_values("time", inplace=True)
    print(df)
    df.to_csv(msd_name, index=False)
