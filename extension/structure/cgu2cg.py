from pyMD.file_parser import LmpParser
import os
import copy


if __name__ == '__main__':
    parser = LmpParser()

    # fp1 = "/home/centos/Projects/CGMD/data/20200722_lr=0.05_1blk_50_cgu_ex_hard/iter_353/cal/"
    # cgu_atom3d = parser.load_data_file(fp1 + "1blk_40.data")
    # cgu_atom3d.renew_coordinate_file(fp1 + "300K.lammpstrj.0")
    # cg_atom3d = copy.deepcopy(cgu_atom3d)
    # cg_atom3d.cgu_to_cg()
    # print(cg_atom3d)
    # cg_atom3d.renew_from_cgu(cgu_atom3d)
    # print(cg_atom3d)
    # cg_atom3d.convert_to_in_cell()
    # print(cg_atom3d)

    fp1 = "/home/centos/ztz/true_strain/x/1E-6/"
    cgu_atom3d = parser.load_data_file(fp1 + "data.8blk_2000")
    cgu_atom3d.renew_coordinate_file(fp1 + "trj/8blk_2000.lammpstrj.s1.0")
    cg_atom3d = cgu_atom3d
    cg_atom3d.cgu_to_cg()
    # print(cg_atom3d)
    # cg_atom3d.renew_from_cgu(cgu_atom3d)
    # print(cg_atom3d)
    cg_atom3d.convert_to_default()
    fp2 = "/home/centos/ztz/true_strain_cg/x/1E-6/"
    LmpParser.create_data_file(cg_atom3d, fp2 + "data.8blk_2000", q=False, improper=False)
