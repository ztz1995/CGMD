from pyMD import ibi
from pyMD import collective_structure_class as csc
from pyMD.file_parser import LmpParser

if __name__ == '__main__':
    # fp1 = "../../data/20200420_lr=0.05_1blk_40_cg_ex_hard/record/"
    fp1 = "/home/centos/Projects/CGMD/data/20200502_SH6S_cg/record/"
    r1 = ibi.Record(fp1, restart=True)
    # param1 = r1.get_param(22, 1)
    param1 = r1.get_param(284, 1)

    # pattern = (["TO(2)"] * 6 + ["Es", "Ph", "U", "Ph", "U", "Ph", "Es"] + ["TO(2)"] * 7) * 7
    # pattern = (["TO(2)"] * 6 + ["Es", "Ph", "U", "Ph", "Me", "Ph", "U", "Ph", "Es"] + ["TO(2)"] * 7) * 7
    pattern = ["TO(2)"] * 1 + ["Es", "Ph", "U", "Ph", "Me", "Ph", "U", "Ph", "Es"] + ["TO(2)"] * 1
    # # pattern = ["TO(2)"] * 6 + ["Es", "Ph", "U", "Ph", "Me", "Ph", "U", "Ph", "Es"] + ["TO(2)"] * 7
    pattern[0] = "TO(1)"
    pattern[-1] = "TO(1)"
    pattern2 = ["TO(2)"] * 11
    pattern2[0] = "TO(1)"
    pattern2[-1] = "TO(1)"
    atom3d = csc.Atom3D()
    # atom3d.create_from_pattern(pattern, 2000, param=param1, density=1.)
    # LmpParser.create_data_file(atom3d, "/home/centos/work/new/cg_7blk_500chain/test/strain.data", q=False, improper=False)
    # LmpParser.create_data_file(atom3d, "/home/centos/work/SH6S/cg_7blk_300chain/init/strain.data", q=False, improper=False)
    # LmpParser.create_data_file(atom3d, "/home/centos/work/SH6S/PDI/cg_7blk_200chain/init/strain.data", q=False, improper=False)
    # LmpParser.create_data_file(atom3d, "/home/centos/work/SH6S/cg_1blk_2000chain/init/strain.data", q=False, improper=False)
    # LmpParser.create_data_file(atom3d, "../../pure.data", q=False, improper=False)

    # atom3d.create_from_pattern(pattern, 2, param=param1, density=1., seed=1)
    atom3d.create_from_patterns([pattern, pattern2], [2, 2], param=param1, density=1., seed=1)
    LmpParser.create_data_file(atom3d, "./test2.data", q=False, improper=False)


    # atom3d.cg_to_cgu()
    # LmpParser.create_data_file(atom3d, "/home/centos/work/SH6S/cgu_7blk_300chain/init/strain.data", q=False, improper=True)
    # LmpParser.create_data_file(atom3d, "/home/centos/work/SH6S/PDI/cgu_7blk_200chain/init/strain.data", q=False, improper=True)
    # LmpParser.create_data_file(atom3d, "/home/centos/work/SH6S/cgu_1blk_2000chain/init/strain.data", q=False, improper=True)

    # LmpParser.create_data_file(atom3d, "cgu_7blk_200chain.data", q=False, improper=True)
    # LmpParser.create_data_file(atom3d, "/home/centos/work/cgu_2blk_800chain/init/strain.data", q=False, improper=True)
    # LmpParser.create_data_file(atom3d, "/home/centos/work/cgu_2blk_800chain_2/init/strain.data", q=False, improper=True)
    # LmpParser.create_data_file(atom3d, "/home/centos/work/cgu_1blk_1000chain/init/strain.data", q=False, improper=True)
    # LmpParser.create_data_file(atom3d, "/home/centos/work/cg_1blk_1000chain/init/strain.data", q=False, improper=False)
