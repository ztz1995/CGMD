from pyMD.file_parser import LmpParser
import os
import numpy as np


if __name__ == '__main__':
    parser = LmpParser()

    fp1 = "/home/centos/work/hyper/cg/S7HS7/S3_G5_10/init/"
    cg_atom3d = parser.load_data_file(fp1 + "hyper.data")
    cg_atom3d.renew_coordinate_file(fp1 + "trj/hyper.lammpstrj.000960000")

    # cg_atom3d.cg_to_cgu(u_orien=[1., 1., 1])
    cg_atom3d.cg_to_cgu()
    fp2 = "/home/centos/work/hyper/cgu/S7HS7/S3_G5_10/anneal/"
    os.makedirs(fp2, exist_ok=True)
    # cg_atom3d.box_h[0] = 175.
    # cg_atom3d.box_l[0] = 0.
    # cg_atom3d.lattice_parameter[0] = 175.
    LmpParser.create_data_file(cg_atom3d, fp2 + "hyper.data", q=False, improper=True)