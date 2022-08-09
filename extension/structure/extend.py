from pyMD.file_parser import LmpParser
from pyMD import collective_structure_class as csc
import os
import copy


def load_seg_id(_fp):
    mole_seg = dict()
    with open(_fp, "r") as file:
        seg_lines = file.readlines()
    for line in seg_lines:
        if line:
            if line[0] != "#":
                args = line.split(" ")
                mole_seg[int(args[0])] = int(args[1])
    return mole_seg


def output(_fp, mole_segment):
    file_lines = ["# mole id -> segment id\n"]
    for k, v in mole_segment.items():
        file_lines.append("%d %d\n" % (k, v))
    with open(_fp, "w") as file:
        file.writelines(file_lines)
    pass


if __name__ == '__main__':

    parser = LmpParser()
    # fp1 = "/home/centos/ztz/stress/line_test/parallel2/equi/"
    # atom3d = parser.load_data_file(fp1 + "data.test_cgu")
    # atom3d.renew_coordinate_file(fp1 + "trj/init.test_cgu.lammpstrj.830000")

    atom3d = parser.load_data_file("data.tilt2")

    new_atom3d = csc.Atom3D()
    mole_seg = load_seg_id("segment.para2")
    # e_mole_seg = new_atom3d.extend_atom3d(atom3d, times=[1, 3, 3], mole_seg=mole_seg, _type="x", lines=2)
    new_atom3d.extend_atom3d(atom3d, times=[3, 8, 3], mole_seg=mole_seg,_type="x")

    # fp2 = "/home/centos/ztz/stress/line_test/parallel8/init/"

    # fp3 = "/home/centos/ztz/stress/line_test/para18/init/"
    fp3 = "/home/centos/ztz/stress/line_test/tilt68/init/"
    # output(fp3 + "segment.para18", e_mole_seg)
    os.makedirs(fp3, exist_ok=True)
    # LmpParser.create_data_file(new_atom3d, fp3+"data.para18", q=False, improper=True)
    LmpParser.create_data_file(new_atom3d, fp3+"data.tilt68", q=False, improper=True)
