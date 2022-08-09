# -*- coding: UTF-8 -*-
import sys

sys.path.append("/home/centos/Projects/CGMD/")
from pyMD.ibi import IBI
from pyMD import collective_structure_class as csc


if __name__ == '__main__':
    # fp = "/home/centos/model/aa/1blk_40/0/"
    # atom3d = csc.create_atom3d(fp, 40)
    #
    # ibi = IBI("20200329_lr=0.1", atom3d)   # iter 36 change to 0.05 lr, iter_119 change to 0.03 lr, iter_139 change dihedral
    # ibi.iterative_fit_params(500, "../1blk_40_target_exclude_hard_convolve/aver/", restart=True, learning_rate=0.01,
    #                          lmp_parallel=24, num_parallel=48)

    # fp = "/home/centos/model/aa/SH6S/0/"
    # atom3d = csc.Atom3D()
    # atom3d.create_from_info(fp + "fg_struct_info.pkl")
    # aa = csc.AA3D()
    # aa.get_mass(fp + "SH6S.data")
    # # aa.renew_coordinate(fp + "SH6S.lammpstrj", step=int(24005000))
    # aa.renew_coordinate(fp + "SH6S.lammpstrj", step=int(16000000))
    # atom3d.renew_coordinate(aa)
    # # atom3d.convert_to_in_cell()
    # ibi = IBI("20200329_lr=0.05_SH6S", atom3d)  # take iter_284
    # ibi.iterative_fit_params(500, "../data/target/SH6S_target/0/", restart=True, learning_rate=0.03,
    #                          lmp_parallel=24, num_parallel=48)
    # # print(atom3d.Bond_type_dict)

    # fp = "/home/centos/model/aa/1blk_40/0/"
    # atom3d = csc.Atom3D()
    # atom3d.create_from_info(fp + "fgh_struct_info.pkl")
    # aa = csc.AA3D()
    # aa.get_mass(fp + "1blk_40.data")
    # # aa.renew_coordinate(fp + "SH6S.lammpstrj", step=int(24005000))
    # aa.renew_coordinate(fp + "1blk_40.lammpstrj", step=int(16000000))
    # atom3d.renew_coordinate(aa)
    # atom3d.cal_charge()
    #
    # # atom3d.convert_to_in_cell()
    # ibi = IBI("20200408_lr=0.05_hard_aa", atom3d)
    # # ibi = IBI("20200408_lr=0.05_hard_aa_pppm", atom3d)
    #
    # ibi.iterative_fit_params(500, "../data/target/1blk_40_hard_aa/aver/", restart=True, learning_rate=0.04,
    #                          lmp_parallel=24, num_parallel=48, q=False)
    # # print(atom3d.Bond_type_dict)

    fp = "/home/centos/model/aa/SH6S_gaff/0/"
    atom3d = csc.Atom3D()
    atom3d.create_from_info(fp + "fgh_struct_info.pkl")
    aa = csc.AA3D()
    aa.get_mass(fp + "SH6S_gaff.data")
    # aa.renew_coordinate(fp + "SH6S.lammpstrj", step=int(24005000))
    aa.renew_coordinate(fp + "SH6S_gaff.lammpstrj", step=int(16000000))
    atom3d.renew_coordinate(aa)
    atom3d.cal_charge()
    # atom3d.convert_to_in_cell()
    ibi = IBI("20200410_lr=0.05_SH6S_fgh", atom3d)
    # ibi.iterative_fit_params(500, "../data/target/SH6S_fgh/0/", restart=True, learning_rate=0.02,
    #                          lmp_parallel=24, num_parallel=48)
    # iter_43 change target
    ibi.iterative_fit_params(500, "../data/target/SH6S_fgh_gaff/0/", restart=True, learning_rate=0.025,
                             lmp_parallel=24, num_parallel=48, q=False)

    # fp = "/home/centos/model/aa/1blk_40/0/"
    # atom3d = csc.Atom3D()
    # atom3d.create_from_info(fp + "fgh_struct_info.pkl")
    # aa = csc.AA3D()
    # aa.get_mass(fp + "1blk_40.data")
    # # aa.renew_coordinate(fp + "SH6S.lammpstrj", step=int(24005000))
    # aa.renew_coordinate(fp + "1blk_40.lammpstrj", step=int(16000000))
    # atom3d.renew_coordinate(aa)
    # # atom3d.cal_charge()
    #
    # # atom3d.convert_to_in_cell()
    # ibi = IBI("test", atom3d)
    #
    # ibi.iterative_fit_params(500, "../data/target/1blk_40_hard_aa/aver/", restart=True, learning_rate=0.04,
    #                          lmp_parallel=24, num_parallel=48, q=False)
    # # print(atom3d.Bond_type_dict)