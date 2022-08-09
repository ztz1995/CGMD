from pyMD import ibi
from pyMD import collective_structure_class as csc
from pyMD.file_parser import LmpParser
import os

if __name__ == '__main__':
    # fp1 = "../../data/20200420_lr=0.05_1blk_40_cg_ex_hard/record/"
    # fp1 = "/home/centos/Projects/CGMD/data/20200502_SH6S_cg/record/"
    fp1 = "/home/centos/Projects/CGMD/data/20200720_lr=0.05_1blk_50_cg_ex_hard/record/"
    r1 = ibi.Record(fp1, restart=True)
    # param1 = r1.get_param(22, 1)
    param1 = r1.get_param(52, 1)

    blocks = 1
    # pattern = (["TO(2)"] * 6 + ["Es", "Ph", "U", "Ph", "Me", "Ph", "U", "Ph", "Es"] + ["TO(2)"] * 7) * blocks
    # pattern[0] = "TO(1)"
    # pattern[-1] = "TO(1)"
    H = ["Ph", "U", "Ph"]
    pattern =  ["TO(2)"] * 6 + ["Es"] +  H + ["Me"] + H + ["Me"]+ H + ["Me"]+ H + ["Es"] + ["TO(2)"] * 7
    # pattern = pattern + ["Es"] + pattern
    pattern[-1] = "TO(1)"
    pattern[0] = "TO(1)"


    # pattern = ["TO(2)"] * 13 + ["Es"] + H + ["Me"] + H + ["Es"] + ["TO(2)"] * 3 + ["Es"] + H + ["Me"] + H + ["Es"]+ ["TO(2)"] * 13
    # pattern[0] = "TO(1)"
    # pattern[-1] = "TO(1)"

    # pattern = ["TO(2)"] * 6 + ["Es"] + H + ["Me"] + H + ["Es"] + ["TO(2)"] * 3 + ["Es"] + H + ["Me"] + H + ["Es"]+ ["TO(2)"] * 7
    # pattern[0] = "TO(1)"
    # pattern[-1] = "TO(1)"

    chains = 800
    atom3d = csc.Atom3D()
    atom3d.create_from_pattern(pattern, chains, param=param1, density=1.)
    fp = "/home/centos/ztz/dpd/AB/S7_H15_S8_800/"
    os.makedirs(fp, exist_ok=True)
    LmpParser.create_data_file(atom3d, fp + "data.%dblk_%d" % (blocks, chains), q=False, improper=False, dihedral=False)

    # LmpParser.create_data_file(atom3d, "/home/centos/work/SH6S/cg_7blk_300chain/init/strain.data", q=False, improper=False)
    # LmpParser.create_data_file(atom3d, "/home/centos/work/SH6S/PDI/cg_7blk_200chain/init/strain.data", q=False, improper=False)
    # LmpParser.create_data_file(atom3d, "/home/centos/work/SH6S/cg_1blk_2000chain/init/strain.data", q=False, improper=False)
    # LmpParser.create_data_file(atom3d, "../../pure.data", q=False, improper=False)

    # atom3d.create_from_pattern(pattern, 2, param=param1, density=1., seed=1)
    # atom3d.create_from_patterns([pattern, pattern2], [2, 2], param=param1, density=1., seed=1)
    # LmpParser.create_data_file(atom3d, "./test2.data", q=False, improper=False)

    # atom3d.cg_to_cgu()
    # fp = "/home/centos/ztz/dpd/cgu/2blk_100/dpd_hybrid/"
    # os.makedirs(fp, exist_ok=True)
    # LmpParser.create_data_file(atom3d, fp + "data.%dblk_%d" % (blocks, chains), q=False, improper=True)
