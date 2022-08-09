# -*- coding: UTF-8 -*-
import sys

project_path = "/home/centos/Projects/CGMD/"
sys.path.append(project_path)
from pyMD.ibi import IBI
from pyMD import collective_structure_class as csc

if __name__ == '__main__':
    # fp = "/home/centos/model/aa/SH6S_gaff/0/"
    # atom3d = csc.Atom3D()
    # atom3d.create_from_info(fp + "cgu_struct_info.pkl")
    # aa = csc.AA3D()
    # aa.get_mass(fp + "SH6S_gaff.data")
    # aa.renew_coordinate(fp + "SH6S_gaff.lammpstrj", step=int(16000000))
    # atom3d.renew_coordinate(aa)
    # ibi = IBI("20200419_lr=0.05_SH6S_cgu_ex_hard", atom3d)
    # ibi.iterative_fit_params(500, project_path + "data/target/SH6S_gaff_cgu_ex_hard/aver/",
    #                          restart=True, learning_rate=0.025,
    #                          lmp_parallel=24, num_parallel=48, improper=True)
    #
    # # get param from cg_ex_hard, iter_184
    # # Ph_n1_c1_o1 dihedral from fgh
    # # iter_247 is good

    # # fp = "/home/centos/model/aa/SH6S_gaff/0/"
    # fp = "/home/centos/model/aa/1blk_40/0/"
    # atom3d = csc.Atom3D()
    # atom3d.create_from_info(fp + "cgu_struct_info.pkl")
    # aa = csc.AA3D()
    # aa.get_mass(fp + "1blk_40.data")
    # aa.renew_coordinate(fp + "1blk_40.lammpstrj", step=int(15000000))
    # atom3d.renew_coordinate(aa)
    # ibi = IBI("20200422_lr=0.05_1blk_40_cgu_ex_hard", atom3d)
    # ibi.iterative_fit_params(500, project_path + "data/target/1blk_40_cgu_ex_hard/aver/",
    #                          restart=True, learning_rate=0.02,
    #                          lmp_parallel=24, num_parallel=48, improper=True)
    # # get param from SH6S cgu iter_247
    # # change target and charge to origin COMPASS
    # # take iter 183

    # fp = "/home/centos/model/aa/SH6S/0_1/"
    # atom3d = csc.Atom3D()
    # atom3d.create_from_info(fp + "cgu_struct_info.pkl")
    # aa = csc.AA3D()
    # aa.get_mass(fp + "SH6S.data")
    # aa.renew_coordinate(fp + "SH6S.lammpstrj", step=int(16000000))
    # atom3d.renew_coordinate(aa)
    # ibi = IBI("20200427_lr=0.05_SH6S_cgu_ex_hard", atom3d)
    # ibi.iterative_fit_params(500, project_path + "data/target/SH6S_1_cgu_ex_hard/0/",
    #                          restart=True, learning_rate=0.025,
    #                          lmp_parallel=24, num_parallel=48, improper=True)
    #
    # # get param from 20200420_lr=0.05_SH6S_cgu_ex_hard, iter_274
    # # iter_32 acceptable

    # # fp = "/home/centos/model/aa/SH6S_gaff/0/"
    # fp = "/home/centos/model/aa/1blk_40/0/"
    # atom3d = csc.Atom3D()
    # atom3d.create_from_info(fp + "cgu_struct_info.pkl")
    # aa = csc.AA3D()
    # aa.get_mass(fp + "1blk_40.data")
    # aa.renew_coordinate(fp + "1blk_40.lammpstrj", step=int(15000000))
    # atom3d.renew_coordinate(aa)
    # ibi = IBI("20200427_lr=0.05_1blk_40_cgu_ex_hard", atom3d)
    # ibi.iterative_fit_params(500, project_path + "data/target/1blk_40_cgu_ex_hard/aver/",
    #                          restart=True, learning_rate=0.02, with_pre=False,
    #                          lmp_parallel=24, num_parallel=48, improper=True, non_bond=False)
    # # get param from 20200427_lr=0.05_SH6S_cgu_ex_hard iter_32
    # #

    # fp = "/home/centos/model/aa/SH6S/0/"
    # atom3d = csc.Atom3D()
    # atom3d.create_from_info(fp + "cgu_struct_info.pkl")
    # aa = csc.AA3D()
    # aa.get_mass(fp + "SH6S.data")
    # aa.renew_coordinate(fp + "SH6S.lammpstrj", step=int(24000000))
    # atom3d.renew_coordinate(aa)
    # ibi = IBI("20200502_SH6S_cgu", atom3d)
    # ibi.iterative_fit_params(500, project_path + "data/target/SH6S_cgu_ex_hard/aver/",
    #                          restart=True, learning_rate=0.025,
    #                          lmp_parallel=24, num_parallel=48, improper=True)
    #
    # # get param from 20200427_lr=0.05_SH6S_cgu_ex_hard/iter_32
    # # iter 234 is good

    fp = "/home/centos/model/aa/1blk_40/0/"
    atom3d = csc.Atom3D()
    atom3d.create_from_info(fp + "cgu_struct_info.pkl")
    aa = csc.AA3D()
    aa.get_mass(fp + "1blk_40.data")
    aa.renew_coordinate(fp + "1blk_40.lammpstrj", step=int(15000000))
    atom3d.renew_coordinate(aa)
    ibi = IBI("20200506_1blk_40_cgu", atom3d)
    ibi.iterative_fit_params(500, project_path + "data/target/1blk_40_cgu_ex_hard/aver/",
                             restart=True, learning_rate=0.02, with_pre=False,
                             lmp_parallel=24, num_parallel=48, improper=True, non_bond=False)
    # get param from 20200502_SH6S_cgu/iter_234, 20200427_lr=0.05_1blk_40_cgu_ex_hard/iter_02
    # iter 02 is good
