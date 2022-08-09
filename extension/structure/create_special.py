from pyMD import ibi
from pyMD import collective_structure_class as csc
from pyMD.file_parser import LmpParser
import os
import numpy as np

if __name__ == '__main__':
    fp1 = "../../data/20200720_lr=0.05_1blk_50_cg_ex_hard/record/"
    r1 = ibi.Record(fp1, restart=True)
    # param1 = r1.get_param(22, 1)
    param1 = r1.get_param(52, 1)
    pattern = (["TO(2)"] * 6 + ["Es", "Ph", "U", "Ph", "Me", "Ph", "U", "Ph", "Es"] + ["TO(2)"] * 7) * 1
    pattern[0] = "TO(1)"
    pattern[-1] = "TO(1)"

    atom3d = csc.Atom3D()
    # atom3d.create_line(pattern, 40, param1, density=1.1, lines=2)
    atom3d.create_line_tilt2(pattern, 14, param1, density=1.1)
    LmpParser.create_data_file(atom3d, "test_cg.data", q=False, improper=False)

    atom3d.cg_to_cgu(u_orien=np.array([np.sqrt(2)/2, 0, np.sqrt(2)/2]))
    atom3d.cg_to_cgu()
    LmpParser.create_data_file(atom3d, "data.tilt4", q=False, improper=True)
