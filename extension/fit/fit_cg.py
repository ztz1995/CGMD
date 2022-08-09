# -*- coding: UTF-8 -*-
import sys

project_path = "/home/centos/Projects/CGMD/"
sys.path.append(project_path)
from pyMD.ibi import IBI
from pyMD import collective_structure_class as csc

if __name__ == '__main__':
    # fp = "/home/centos/model/aa/SH6S_gaff/0/"
    # atom3d = csc.Atom3D()
    # atom3d.create_from_info(fp + "cg_struct_info.pkl")
    # aa = csc.AA3D()
    # aa.get_mass(fp + "SH6S_gaff.data")
    # aa.renew_coordinate(fp + "SH6S_gaff.lammpstrj", step=int(16000000))
    # atom3d.renew_coordinate(aa)
    # # ibi = IBI("20200418_lr=0.05_SH6S_cg", atom3d)
    # ibi = IBI("20200418_lr=0.05_SH6S_cg_ex_hard", atom3d)
    # ibi.iterative_fit_params(500, project_path + "data/target/SH6S_gaff_cg_ex_hard/aver/",
    #                          restart=True, learning_rate=0.02,
    #                          lmp_parallel=24, num_parallel=48, improper=False)
    # # iter29 change soft exclude
    # # iter56 change lr to 0.025
    # # iter158 a good param, iter_184
    # # iter_289 even better

    # fp = "/home/centos/model/aa/1blk_40/0/"
    # atom3d = csc.Atom3D()
    # atom3d.create_from_info(fp + "cg_struct_info.pkl")
    # aa = csc.AA3D()
    # aa.get_mass(fp + "1blk_40.data")
    # aa.renew_coordinate(fp + "1blk_40.lammpstrj", step=int(15000000))
    # atom3d.renew_coordinate(aa)
    # # ibi = IBI("20200418_lr=0.05_SH6S_cg", atom3d)
    # ibi = IBI("20200420_lr=0.05_1blk_40_cg_ex_hard", atom3d)
    # ibi.iterative_fit_params(500, project_path + "data/target/1blk_40_cg_ex_hard/aver/",
    #                          restart=True, learning_rate=0.05,
    #                          lmp_parallel=24, num_parallel=48, improper=False, with_pre=False)
    # # get param from 20200418_lr=0.05_SH6S_cg_ex_hard/iter_289
    # # do not change non_bond param
    # # iter_22 is good
    # # take iter_32

    # fp = "/home/centos/model/aa/1blk_40/0/"
    # atom3d = csc.Atom3D()
    # atom3d.create_from_info(fp + "cg_struct_info.pkl")
    # aa = csc.AA3D()
    # aa.get_mass(fp + "1blk_40.data")
    # aa.renew_coordinate(fp + "1blk_40.lammpstrj", step=int(15000000))
    # atom3d.renew_coordinate(aa)
    # # ibi = IBI("20200418_lr=0.05_SH6S_cg", atom3d)
    # ibi = IBI("20200420_lr=0.05_1blk_40_cg_ex_hard", atom3d)
    # ibi.iterative_fit_params(500, project_path + "data/target/1blk_40_cg_ex_hard/aver/",
    #                          restart=True, learning_rate=0.05,
    #                          lmp_parallel=24, num_parallel=48, improper=False, with_pre=False)
    # # get param from 20200418_lr=0.05_SH6S_cg_ex_hard/iter_289
    # # do not change non_bond param
    # # iter_22 is good
    # # take iter_32

    # fp = "/home/centos/model/aa/SH6S/0_1/"
    # atom3d = csc.Atom3D()
    # atom3d.create_from_info(fp + "cg_struct_info.pkl")
    # aa = csc.AA3D()
    # aa.get_mass(fp + "SH6S.data")
    # aa.renew_coordinate(fp + "SH6S.lammpstrj", step=int(16000000))
    # atom3d.renew_coordinate(aa)
    # ibi = IBI("20200429_lr=0.05_SH6S_cg_ex_hard", atom3d)
    # # # ibi = IBI("20200420_lr=0.05_SH6S_cgu_ex_hard_1", atom3d)
    # ibi.iterative_fit_params(500, project_path + "data/target/SH6S_1_cg_ex_hard/0/",
    #                          restart=True, learning_rate=0.025,
    #                          lmp_parallel=24, num_parallel=48, improper=False)
    # # get param from 20200418_lr=0.05_SH6S_cg_ex_hard/iter_289
    # # iter 108 is good

    # fp = "/home/centos/model/aa/1blk_40/0/"
    # atom3d = csc.Atom3D()
    # atom3d.create_from_info(fp + "cg_struct_info.pkl")
    # aa = csc.AA3D()
    # aa.get_mass(fp + "1blk_40.data")
    # aa.renew_coordinate(fp + "1blk_40.lammpstrj", step=int(15000000))
    # atom3d.renew_coordinate(aa)
    # # ibi = IBI("20200418_lr=0.05_SH6S_cg", atom3d)
    # ibi = IBI("20200429_lr=0.05_1blk_40_cg_ex_hard", atom3d)
    # ibi.iterative_fit_params(500, project_path + "data/target/1blk_40_cg_ex_hard/aver/",
    #                          restart=True, learning_rate=0.05, non_bond=False,
    #                          lmp_parallel=24, num_parallel=48, improper=False, with_pre=False)
    # # get param from 20200429_lr=0.05_SH6S_cg_ex_hard/iter_108
    # # iter29 is good

    # fp = "/home/centos/model/aa/SH6S/0/"
    # atom3d = csc.Atom3D()
    # atom3d.create_from_info(fp + "cg_struct_info.pkl")
    # aa = csc.AA3D()
    # aa.get_mass(fp + "SH6S.data")
    # aa.renew_coordinate(fp + "SH6S.lammpstrj", step=int(24000000))
    # atom3d.renew_coordinate(aa)
    # ibi = IBI("20200502_SH6S_cg", atom3d)
    # ibi.iterative_fit_params(500, project_path + "data/target/SH6S_cg_ex_hard/aver/",
    #                          restart=True, learning_rate=0.02,
    #                          lmp_parallel=24, num_parallel=48, improper=False)
    # # get param from 20200429_lr=0.05_SH6S_cg_ex_hard/iter_108
    # # iter 284 is good

    fp = "/home/centos/model/aa/1blk_40/0/"
    atom3d = csc.Atom3D()
    atom3d.create_from_info(fp + "cg_struct_info.pkl")
    aa = csc.AA3D()
    aa.get_mass(fp + "1blk_40.data")
    aa.renew_coordinate(fp + "1blk_40.lammpstrj", step=int(15000000))
    atom3d.renew_coordinate(aa)
    ibi = IBI("20200506_1blk_40_cg", atom3d)
    ibi.iterative_fit_params(500, project_path + "data/target/1blk_40_cg_ex_hard/aver/",
                             restart=True, learning_rate=0.05, non_bond=False,
                             lmp_parallel=24, num_parallel=48, improper=False, with_pre=False)
    # get param from 20200502_SH6S_cg/iter_284, 20200429_lr=0.05_1blk_40_cg_ex_hard/iter_29
    # iter 2 is good
