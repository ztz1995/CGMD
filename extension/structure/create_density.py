from pyMD import ibi
from pyMD import collective_structure_class as csc
from pyMD.file_parser import LmpParser
from pyMD import functions as f
import os


def gaussian_closet(mu, s2):
    def g(_x):
        return f.pdf_gaussian(_x, mu, s2)

    return g


def gaussian_2(_x):
    return (f.pdf_gaussian(_x, 0.25, 0.01) + f.pdf_gaussian(_x, 0.75, 0.01)) / 2


if __name__ == '__main__':
    # fp1 = "/home/centos/Projects/CGMD/data/20200502_SH6S_cg/record/"
    # r1 = ibi.Record(fp1, restart=True)
    # param1 = r1.get_param(284, 1)

    fp1 = "../../data/20200420_lr=0.05_1blk_40_cg_ex_hard/record/"
    r1 = ibi.Record(fp1, restart=True)
    param1 = r1.get_param(22, 1)

    pattern = ["TO(2)"] * 6 + ["Es", "Ph", "U", "Ph", "Me", "Ph", "U", "Ph", "Es"] + ["TO(2)"] * 7
    pattern[0] = "TO(1)"
    pattern[-1] = "TO(1)"

    atom3d = csc.Atom3D()
    atom3d.create_from_pattern_center_p(pattern, 200, param=param1, density=1., center=10,
                                        # p=[f.pdf_uniform] + [gaussian_closet(0.5, 0.02)] * 2)
                                        p=[f.pdf_uniform] + [gaussian_2] * 2)
    cg_fp = "/home/centos/work/change_density/cg_1blk_200chain/init_0_2_2/"
    if not os.path.exists(cg_fp):
        os.makedirs(cg_fp)
    LmpParser.create_data_file(atom3d, cg_fp + "strain.data", q=False, improper=False)

    atom3d.cg_to_cgu()
    cgu_fp = "/home/centos/work/change_density/cgu_1blk_200chain/init_0_2_2/"
    if not os.path.exists(cgu_fp):
        os.makedirs(cgu_fp)
    LmpParser.create_data_file(atom3d, cgu_fp + "strain.data", q=False, improper=True)
