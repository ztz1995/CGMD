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

    # # get param from 20200506_1blk_40_cgu/iter_02
    # ibi = IBI("20200728_lr=0.05_1blk_50_cgu_ex_hard", atom3d_gen=atom3d_gen)
    # ibi.iterative_fit_params(500, project_path + "data/target/1blk_50_cgu_ex_hard/aver/",
    #                          restart=True, learning_rate=0.05,
    #                          lmp_parallel=24, num_parallel=48, improper=True)

    # ibi = IBI("20200730_lr=0.05_1blk_50_cgu_ex_hard", atom3d_gen=atom3d_gen)
    # ibi.iterative_fit_params(500, project_path + "data/target/1blk_50_cgu_ex_hard/aver/",
    #                          restart=False, learning_rate=0.1,
    #                          lmp_parallel=24, num_parallel=48, improper=False)

    # ibi = IBI("20200915_lr=0.05_1blk_50_600K_cgu_ex_hard", atom3d_gen=atom3d_gen)
    # ibi.iterative_fit_params(500, project_path + "data/target/1blk_50_600K_cgu_ex_hard/aver/",
    #                          restart=True, learning_rate=0.5,
    #                          lmp_parallel=24, num_parallel=48, improper=True)

    # 20200914_lr=0.05_1blk_50_600K_cg_ex_hard/iter 39
    ibi = IBI("20200917_lr=0.05_1blk_50_600K_cgu_ex_hard", atom3d_gen=atom3d_gen)
    ibi.iterative_fit_params(500, project_path + "data/target/1blk_50_600K_cgu_ex_hard/aver/",
                             restart=True, learning_rate=0.05,
                             lmp_parallel=24, num_parallel=48, improper=True)