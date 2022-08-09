from pyMD import ibi
from pyMD import collective_structure_class as csc
from pyMD.file_parser import LmpParser
import os
import numpy as np

if __name__ == '__main__':
    # fp1 = "../../data/20200420_lr=0.05_1blk_40_cg_ex_hard/record/"
    # fp1 = "/home/centos/Projects/CGMD/data/20200502_SH6S_cg/record/"
    fp1 = "../../data/20200720_lr=0.05_1blk_50_cg_ex_hard/record/"
    r1 = ibi.Record(fp1, restart=True)
    # param1 = r1.get_param(22, 1)
    param1 = r1.get_param(52, 1)
    # param1["Bond"][("sv1", "sv2")] = np.array([2.68282353, 2.0, 8.20071071])

    blocks = 8
    chains = 2000
    # sn = 8000
    # pattern = (["TO(2)"] * 6 + ["Es", "Ph", "U", "Ph", "U", "Ph", "Es"] + ["TO(2)"] * 7) * 7
    pattern = (["TO(2)"] * 6 + ["Es", "Ph", "U", "Ph", "Me", "Ph", "U", "Ph", "Es"] + ["TO(2)"] * 7) * blocks
    # pattern = ["TO(2)"] * 1 + ["Es", "Ph", "U", "Ph", "Me", "Ph", "U", "Ph", "Es"] + ["TO(2)"] * 1
    # pattern = ["TO(2)"] * 6 + ["Es", "Ph", "U", "Ph", "Me", "Ph", "U", "Ph", "Es"] + ["TO(2)"] * 7
    # for i in range(5, 16):
    # pattern = ["TO(2)"] * 15
    # pattern[0] = "Es"
    # pattern[-1] = "Es"
    # pattern2 = ["TO(2)"] * 11
    pattern[0] = "TO(1)"
    pattern[-1] = "TO(1)"
    atom3d = csc.Atom3D()

    # atom3d.create_from_patterns([pattern, ["sv1", "sv2"]], [chains, sn], param=param1, density=0.5)
    atom3d.create_from_patterns([pattern], [chains], param=param1, density=0.5)

    # fp = "/home/centos/ztz/solvent/test2/"
    # os.makedirs(fp, exist_ok=True)

    atom3d.cg_to_cgu()
    # fp = "/home/centos/ztz/stress/cgu/1blk_%dchain/init/" % chains
    # fp = "/home/centos/ztz/change_h_bond/%dblk_%d/" % (blocks, chains)
    fp = "/home/centos/ztz/stress/phase_structure/%dblk_%d/init/" % (blocks, chains)

    os.makedirs(fp, exist_ok=True)
    # LmpParser.create_data_file(atom3d, fp + "%dblk_%d_Sv_%d.data" % (blocks, chains, sn), q=False, improper=True)
    LmpParser.create_data_file(atom3d, fp + "data.%dblk_%d" % (blocks, chains), q=False, improper=True)

