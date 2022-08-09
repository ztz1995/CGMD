import sys, os
project_root = os.path.abspath(__file__ + "../../../../")
sys.path.append(project_root)
from pyMD import ibi
from pyMD import collective_structure_class as csc
from pyMD.file_parser import LmpParser
import numpy as np



if __name__ == '__main__':

    fp1 = project_root + "/data/20200720_lr=0.05_1blk_50_cg_ex_hard/record/"
    r1 = ibi.Record(fp1, restart=True)
    param1 = r1.get_param(52, 1)

    # print(param1["Angle"])
    param1["Angle"][ ('Ph', 'U', 'Ph')] = np.array([63.52, 117.4])
    param1["Angle"][ ('U', 'Ph', 'U')] = np.array([48.50302414, 165.23629667])

    blocks = 8

    H = ["Ph", "U"]
    # pattern =  (["TO(2)"] * 6 + ["Es"] + H * 4 + ["Ph", "Es"] + ["TO(2)"] * 7) * blocks
    # pattern1 =  ["TO(2)"] * 6 + ["Es"] + H * 2 + ["Ph", "Es"] + ["TO(2)"] * 7
    # pattern2 =  ["TO(2)"] * 6 + ["Es"] + H * 6 + ["Ph", "Es"] + ["TO(2)"] * 7
    # pattern = (pattern1 + pattern2) * (blocks//2)

    pattern1 =  ["TO(2)"] * 6 + ["Es"] + H * 2 + ["Ph", "Es"] + ["TO(2)"] * 7
    pattern2 =  ["TO(2)"] * 6 + ["Es"] + H * 10 + ["Ph", "Es"] + ["TO(2)"] * 7
    pattern = (pattern1 * 3 + pattern2) * (blocks//4)

    pattern[-1] = "TO(1)"
    pattern[0] = "TO(1)"

    chains = 400
    atom3d = csc.Atom3D()
    atom3d.create_from_pattern(pattern, chains, param=param1, density=1.05)

    fp = "/home/centos/ztz/dpd/dipersity/AB_22210_8/"
    os.makedirs(fp + "trj/", exist_ok=True)
    LmpParser.create_data_file(atom3d, fp + "data.%dblk_%d" % (blocks, chains), improper=False,dihedral=False)


    # ibi_ob = ibi.IBI("create", atom3d)
    # ibi_ob.improper = False
    # ibi_ob.dihedral = False
    # ibi_ob.coefficient = param1
    # ibi_ob.prepare_simulation_file(fp, "nvt.in")


