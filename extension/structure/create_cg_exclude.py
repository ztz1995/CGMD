from pyMD import ibi
from pyMD import collective_structure_class as csc
from pyMD.file_parser import LmpParser

if __name__ == '__main__':
    fp1 = "../../data/20200420_lr=0.05_1blk_40_cg_ex_hard/record/"
    # fp1 = "/home/centos/Projects/CGMD/data/20200502_SH6S_cg/record/"
    r1 = ibi.Record(fp1, restart=True)
    param1 = r1.get_param(22, 1)
    # param1 = r1.get_param(284, 1)

    # pattern = (["TO(2)"] * 6 + ["Es", "Ph", "U", "Ph", "U", "Ph", "Es"] + ["TO(2)"] * 7) * 7
    pattern = (["TO(2)"] * 6 + ["Es", "Ph", "U", "Ph", "Me", "Ph", "U", "Ph", "Es"] + ["TO(2)"] * 7) * 7
    # pattern = ["TO(2)"] * 1 + ["Es", "Ph", "U", "Ph", "Me", "Ph", "U", "Ph", "Es"] + ["TO(2)"] * 1
    # pattern = ["TO(2)"] * 6 + ["Es", "Ph", "U", "Ph", "Me", "Ph", "U", "Ph", "Es"] + ["TO(2)"] * 7
    pattern[0] = "TO(1)"
    pattern[-1] = "TO(1)"
    pattern2 = ["TO(2)"] * 11
    pattern2[0] = "TO(1)"
    pattern2[-1] = "TO(1)"
    atom3d = csc.Atom3D()

    # atom3d.create_from_patterns_exclude([pattern, pattern2], [20, 20], param=param1, density=1., seed=1)
    atom3d.create_from_patterns_exclude([pattern], [200], param=param1, density=1.)
    # LmpParser.create_data_file(atom3d, "./test2.data", q=False, improper=False)
    LmpParser.create_data_file(atom3d, "/home/centos/work/electric/cg_7blk_200chain/init/strain.data", q=False, improper=False)


    atom3d.cg_to_cgu()
    LmpParser.create_data_file(atom3d, "/home/centos/work/electric/cgu_7blk_200chain/init/strain.data", q=False, improper=True)
