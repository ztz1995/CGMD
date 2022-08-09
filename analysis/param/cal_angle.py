from pyMD.file_parser import LmpParser
from pyMD.ibi import IBI
import os
import copy as cp
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pyMD import analysis as ans
from multiprocessing import Pool


def cal_all_param_single(x, y, cg):
    for joint_type in ["Bond", "Angle", "Dihedral"]:
        for tt in cg.__dict__[joint_type + "_type_dict"].keys():
            y[joint_type][tt] += cg.cal_intra_para(joint_type, tt, x[joint_type])
    for i in range(len(ats)):
        for j in range(i, len(ats)):
            att = tuple(sorted([ats[i], ats[j]]))
            y["non_bond"][att] += cg.cal_inter_para(att, x["non_bond"])[:-200]
    print("finish")
    return y


if __name__ == '__main__':
    frames = 80
    steps = 10000

    fp = "/home/centos/work/1blk_160/cgu_1blk_160chain/equi/"
    fp_data = fp + "1blk_160.data"
    fp_trj = fp + "trj/1blk_160.lammpstrj."
    atom3ds = ans.Atom3dGenerator(fp_data, fp_trj, one_file=False, pad=9)

    cg_atom3d = cp.deepcopy(atom3ds())
    cg_atom3d.cgu_to_cg()
    # cg_atom3d.sort_atom_bond_id()

    x = IBI.cal_points
    y = dict()
    for joint_type in ["Bond", "Angle", "Dihedral"]:
        y[joint_type] = dict()
        for tt in cg_atom3d.__dict__[joint_type + "_type_dict"].keys():
            y[joint_type][tt] = np.zeros(len(x[joint_type]))

    ats = list(cg_atom3d.Atom_type_dict.keys())
    y["non_bond"] = dict()
    for i in range(len(ats)):
        for j in range(i, len(ats)):
            y["non_bond"][tuple(sorted([ats[i], ats[j]]))] = np.zeros(len(x["non_bond"]))

    # x = IBI.cal_points["Angle"]
    # y = np.zeros(len(x))
    pool = Pool(10)
    ret = list()
    volume = 0.
    for i in range(frames):
        cgu_atom3d = atom3ds(i * steps)
        cg_atom3d.renew_from_cgu(cgu_atom3d)
        volume += np.power(cgu_atom3d.lattice_parameter[0], 3)
        ret.append(pool.apply_async(cal_all_param_single, args=(x, y, cg_atom3d)))
    pool.close()
    pool.join()

    for r in ret:
        new_y = r.get()
        for joint_type in y.keys():
            for tt in y[joint_type].keys():
                y[joint_type][tt] += new_y[joint_type][tt]

    for joint_type in y.keys():
        for tt in y[joint_type].keys():
            z = y[joint_type][tt] / frames
            if joint_type == "non_bond":
                # print(len(cg_atom3d.Atom_type_dict[tt[1]]))
                # print(volume/frames)
                z /= (len(cg_atom3d.Atom_type_dict[tt[1]]) / volume * frames)
                pass
            else:
                z /= z.sum() * (x[joint_type][1] - x[joint_type][0])
            plt.plot(x[joint_type], z)
            plt.title(joint_type + " " + "-".join(tt))
            plt.show()
            data = pd.DataFrame({"x": x[joint_type], "y": z})
            data.to_csv("data/%s_%s_cgu2cg.csv" % (joint_type, "_".join(tt)))

    # x = IBI.cal_points["Angle"]
    # y = np.zeros(len(x))
    #
    # for i in range(frames):
    #     print(i * steps)
    #     cgu_atom3d = atom3ds(i * steps, only_type=["Ph", "c1", "h1", "o1", "n1"])
    #     cg_atom3d.renew_from_cgu(cgu_atom3d)
    #     y += cg_atom3d.cal_intra_para("Angle", ("Ph", "U", "Ph"), x)
    # y /= y.sum() * (x[1] - x[0])
    # data = pd.DataFrame({"x": x, "y":y})
    # data.to_csv("data/Angle_Ph_U_Ph_cgu.csv")
    #
    # plt.plot(x, y)
    # plt.show()

    # fp = "/home/centos/work/1blk_50/cg_1blk_160chain/equi/"
    # fp_data = fp + "1blk_160.data"
    # fp_trj = fp + "trj/1blk_160.lammpstrj."
    # atom3ds = ans.Atom3dGenerator(fp_data, fp_trj, one_file=False, pad=9)
    # cg_atom3d = atom3ds()
    # # print(cg_atom3d)
    # angles = cp.deepcopy(cg_atom3d.Angles)
    # angle_type_dict = cp.deepcopy(cg_atom3d.Angle_type_dict)
    # # print(angles)
    #
    # x = IBI.cal_points["Angle"]
    # y = np.zeros(len(x))
    #
    # for i in range(frames):
    #     print(i * steps)
    #     cg_atom3d = atom3ds(i * steps, only_type=["Ph", "U"])
    #     cg_atom3d.Angles = angles
    #     cg_atom3d.Angle_type_dict = angle_type_dict
    #     y += cg_atom3d.cal_intra_para("Angle", ("Ph", "U", "Ph"), x)
    # y /= y.sum() * (x[1] - x[0])
    # data = pd.DataFrame({"x": x, "y":y})
    # data.to_csv("data/Angle_Ph_U_Ph_cg.csv")
    #
    # plt.plot(x, y)
    # plt.show()
