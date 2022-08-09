# -*- coding: UTF-8 -*-
import sys

project_path = "/home/centos/Projects/CGMD/"
sys.path.append(project_path)
from pyMD.ibi import IBI
from pyMD import collective_structure_class as csc
from pyMD.file_parser import LmpParser
import numpy as np


def atom3d_gen():
    folder_list = [0,1,2,3,5,7,8,9]
    rand_num = np.random.randint(0, 1)
    fp = "/home/centos/ztz/model/aa/1blk_50/%d/" % (folder_list[rand_num])
    atom3d = csc.Atom3D()
    atom3d.create_from_info(fp + "struct_info.pkl", improper="ignore")
    aa = csc.AA3D()
    aa.get_mass(fp + "1blk_50.data")
    aa.renew_coordinate(fp + "1blk_50.lammpstrj", step=int(15000000))
    atom3d.renew_coordinate(aa)
    return atom3d


if __name__ == '__main__':

    ibi = IBI("20210401_lr=0.05_1blk_50_fg_ex_hard", atom3d_gen=atom3d_gen)
    ibi.iterative_fit_params(500, project_path + "data/target/1blk_50_fg_ex_hard_new/aver/",
                             restart=True, learning_rate=0.02,
                             lmp_parallel=24, num_parallel=40, improper=True)

    # iter 73 is good
    # fix error and rerun (recal target, and correct neighbor)
    # iter 169 is good