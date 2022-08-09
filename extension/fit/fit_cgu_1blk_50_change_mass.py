# -*- coding: UTF-8 -*-
import sys

project_path = "/home/centos/Projects/CGMD/"
sys.path.append(project_path)
from pyMD.ibi import IBI
from pyMD import collective_structure_class as csc


def atom3d_gen():
    from pyMD.ibi import Record
    fp1 = "/home/centos/Projects/CGMD/data/20200719_lr=0.05_1blk_40_cg_ex_hard/record/"
    r1 = Record(fp1, restart=True)
    param1 = r1.get_param(62, 1)
    pattern = ["TO(2)"] * 6 + ["Es", "Ph", "U", "Ph", "Me", "Ph", "U", "Ph", "Es"] + ["TO(2)"] * 7
    pattern[0] = "TO(1)"
    pattern[-1] = "TO(1)"
    atom3d = csc.Atom3D()
    # atom3d.create_from_pattern(pattern, 40, param=param1, density=1.)
    atom3d.create_from_pattern(pattern, 50, param=param1, density=1.)
    atom3d.cg_to_cgu()
    return atom3d


if __name__ == '__main__':

    # ibi = IBI("20200722_lr=0.05_1blk_50_cgu_ex_hard", atom3d_gen=atom3d_gen)
    # ibi.iterative_fit_params(500, project_path + "data/target/1blk_50_cgu_ex_hard/aver/",
    #                          restart=True, learning_rate=0.025,
    #                          lmp_parallel=24, num_parallel=48, improper=True)
    #   # iter 353 is good

    # # get param from 20200506_1blk_40_cgu/iter_02
    # ibi = IBI("20200728_lr=0.05_1blk_50_cgu_ex_hard", atom3d_gen=atom3d_gen)
    # ibi.iterative_fit_params(500, project_path + "data/target/1blk_50_cgu_ex_hard/aver/",
    #                          restart=True, learning_rate=0.05,
    #                          lmp_parallel=24, num_parallel=48, improper=True)

    # ibi = IBI("20200730_lr=0.05_1blk_50_cgu_ex_hard", atom3d_gen=atom3d_gen)
    # ibi.iterative_fit_params(500, project_path + "data/target/1blk_50_cgu_ex_hard/aver/",
    #                          restart=False, learning_rate=0.1,
    #                          lmp_parallel=24, num_parallel=48, improper=False)

    # get param from 20200722_lr=0.05_1blk_50_cgu_ex_hard/iter_353
    csc.Collection.type_mass_dict = {"TO(1)": 998.0477637143156, "TO(2)": 915.5593296405326, "TO(3)": 73.1142,
                                     "Es": 628.6932966999788, "Me": 14.0267,
                                     "Ph": 754.0837864388974, "U": 340.94080802767706,
                                     "c1": 12.0107, "o1": 15.9994, "n1": 14.0067, "h1": 1.008}
    ibi = IBI("20200826_lr=0.05_1blk_50_cgu_ex_hard_change_mass", atom3d_gen=atom3d_gen)
    ibi.iterative_fit_params(500, project_path + "data/target/1blk_50_cgu_ex_hard/aver/",
                             restart=True, learning_rate=0.025,
                             lmp_parallel=24, num_parallel=48, improper=True, with_pre=False)
