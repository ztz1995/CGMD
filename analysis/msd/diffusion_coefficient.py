# -*- coding: UTF-8 -*-
import sys

sys.path.append("/home/centos/Projects/CGMD/")
from pyMD import analysis as ans


if __name__ == '__main__':

    # fp_data = "/home/centos/Projects/CGMD/data/20200329_lr=0.1/iter_59/cal/1blk_40.data"
    # fp_trj = "/home/centos/Projects/CGMD/data/20200329_lr=0.1/iter_59/cal/300K.lammpstrj."

    # fp = "/home/centos/Projects/CGMD/data/20200506_1blk_40_cg/iter_02/cal/"
    # fp = "/home/centos/Projects/CGMD/data/20200506_1blk_40_cgu/iter_02/cal/"
    # fp = "/home/centos/Projects/CGMD/data/20200502_SH6S_cgu/iter_234/cal/"
    # fp_data = fp + "1blk_40.data"
    # fp_trj = fp + "300K.lammpstrj."
    #
    # atom3ds = ans.atom3d_generator(fp_data, fp_trj, start=0, step=500, frames=1000)
    # dc = ans.DynamicAnalyzer.diffusion_coefficient(atom3ds, 500, 1)
    # print(dc)

    # fp = "/home/centos/work/SH6S/cg_7blk_300chain/anneal_400/"
    # start = 16000000
    #
    # fp_data = fp + "strain.data"
    # fp_trj = fp + "record.lammpstrj."
    #
    # atom3ds = ans.atom3d_generator(fp_data, fp_trj, start=start, step=1000, frames=1000)
    # dc = ans.DynamicAnalyzer.diffusion_coefficient(atom3ds, 1000, 1)

    fp = "/home/centos/model/aa/SH6S/1/"
    fp_data = fp + "SH6S.data"
    fp_trj = fp + "SH6S.lammpstrj"
    fp_info = fp + "cgu_struct_info.pkl"

    atom3ds = ans.atom3d_generator_from_aa(fp_info, fp_data, fp_trj, start=0, step=1000, frames=1000, one_file=False)
    dc = ans.DynamicAnalyzer.diffusion_coefficient(atom3ds, 1000, 1)
    print(dc)