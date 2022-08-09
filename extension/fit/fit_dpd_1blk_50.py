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
    return atom3d


if __name__ == '__main__':


    # get param from data/target/1blk_50_cgu_ex_hard/aver/20200722_lr=0.05_1blk_50_cg_ex_hard/iter_53
    # refine bonded potetial for dpd simulation
    ibi = IBI("20211120_lr=0.05_1blk_50_600K_dpd", atom3d_gen=atom3d_gen)
    ibi.iterative_fit_params(500, project_path + "data/target/1blk_50_600K_cg_ex_hard/aver/",
                             restart=True, learning_rate=0.05,
                             non_bond=True, ensemble="dpd", 
                             lmp_parallel=24, num_parallel=48, 
                             improper=False, with_pre=False)
 
