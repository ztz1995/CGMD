from pyMD import ibi
from pyMD import collective_structure_class as csc
from pyMD.file_parser import LmpParser
import os

if __name__ == '__main__':
    fp1 = "../../data/20200420_lr=0.05_1blk_40_cg_ex_hard/record/"
    # fp1 = "/home/centos/Projects/CGMD/data/20200502_SH6S_cg/record/"
    r1 = ibi.Record(fp1, restart=True)
    param1 = r1.get_param(22, 1)
    # param1 = r1.get_param(284, 1)

    pattern = (["TO(2)"] * 6 + ["Es", "Ph", "U", "Ph", "Me", "Ph", "U", "Ph", "Es"] + ["TO(2)"] * 7) * 7
    # pattern = ["TO(2)"] * 1 + ["Es", "Ph", "U", "Ph", "Me", "Ph", "U", "Ph", "Es"] + ["TO(2)"] * 1
    # # pattern = ["TO(2)"] * 6 + ["Es", "Ph", "U", "Ph", "Me", "Ph", "U", "Ph", "Es"] + ["TO(2)"] * 7
    pattern[0] = "TO(1)"
    pattern[-1] = "TO(1)"

    atom3d = csc.Atom3D()
    atom3d.create_from_pattern_wall(pattern, 300, param=param1, density=1., wall=(1, 0, 0))
    out_path = "/home/centos/work/1blk_50/cg_7blk_300chain_wallx/init/"
    os.makedirs(out_path, exist_ok=True)
    LmpParser.create_data_file(atom3d, out_path + "strain.data", q=False, improper=False)


    # atom3d.cg_to_cgu()
    # LmpParser.create_data_file(atom3d, "/home/centos/work/SH6S/cgu_7blk_300chain/init/strain.data", q=False, improper=True)
    # LmpParser.create_data_file(atom3d, "/home/centos/work/SH6S/PDI/cgu_7blk_200chain/init/strain.data", q=False, improper=True)
    # LmpParser.create_data_file(atom3d, "/home/centos/work/SH6S/cgu_1blk_2000chain/init/strain.data", q=False, improper=True)

    # LmpParser.create_data_file(atom3d, "cgu_7blk_200chain.data", q=False, improper=True)
    # LmpParser.create_data_file(atom3d, "/home/centos/work/cgu_2blk_800chain/init/strain.data", q=False, improper=True)
    # LmpParser.create_data_file(atom3d, "/home/centos/work/cgu_2blk_800chain_2/init/strain.data", q=False, improper=True)
    # LmpParser.create_data_file(atom3d, "/home/centos/work/cgu_1blk_1000chain/init/strain.data", q=False, improper=True)
    # LmpParser.create_data_file(atom3d, "/home/centos/work/cg_1blk_1000chain/init/strain.data", q=False, improper=False)
