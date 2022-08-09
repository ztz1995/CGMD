from pyMD.file_parser import LmpParser
from pyMD.collective_structure_class import Atom3D
from pyMD.basic_structure_class import Bond
import os
import copy


def get_end(atom3d: Atom3D, mole_id):
    mole_ids = atom3d.mole_dict[mole_id]
    to0 = 0
    for idx in mole_ids:
        # print(atom3d.Atoms[idx].type)
        if atom3d.Atoms[idx].type == "TO(1)":
            to0 = idx
            break
    if to0 == 0:
        raise KeyError("No TO(1)!")
    return to0


def get_mole_list(atom3d: Atom3D, end_id):
    visited = atom3d.traverse_bfs_visit(end_id)
    return visited


def get_bond(atom3d: Atom3D, a1, a2):
    for b in atom3d.Atoms[a1].Bonds:
        if atom3d.Bonds[b].get_other_atom_id(a1) == a2:
            return b
    raise KeyError("No such bond! %d %d" % (a1, a2))


def split_mole(atom3d: Atom3D, mole_list, mole_max, mole_id):
    def change_atom_bond_type(to1, to2):
        btt_id = get_bond(atom3d, to1, to2)
        atom3d.Atoms[to1].type = "TO(1)"
        atom3d.Bonds[btt_id] = Bond(btt_id, to1, to2, _type=("TO(1)", "TO(2)"))

    # print(len(mole_list))
    l = int(len(mole_list) / 2)
    to0, to1 = mole_list[l - 1], mole_list[l]
    change_atom_bond_type(to0, mole_list[l - 2])
    change_atom_bond_type(to1, mole_list[l + 1])

    b_id = get_bond(atom3d, to0, to1)
    atom3d.delete_bond(b_id)
    atom3d.mole_dict[mole_id] = mole_list[:l]
    atom3d.mole_dict[mole_max + 1] = mole_list[l:]
    for a_id in mole_list[l:]:
        atom3d.Atoms[a_id].mole_id = mole_max + 1


def split_all_mole(atom3d: Atom3D):
    mole_max = max(atom3d.mole_dict.keys())
    all_mole = list(atom3d.mole_dict.keys())
    for mole_id in all_mole:
        end_id = get_end(atom3d, mole_id)
        mole_list = get_mole_list(atom3d, end_id)
        split_mole(atom3d, mole_list, mole_max, mole_id)
        mole_max += 1
    atom3d.clear_all_joints()
    atom3d.sort_atom_bond_id()
    atom3d.cal_all_joints(improper="ignore")


if __name__ == '__main__':
    parser = LmpParser()
    fp1 = "/home/centos/ztz/stress/small_split/8blk_200/init/"
    cgu_atom3d = parser.load_data_file(fp1 + "data.8blk_200")
    cgu_atom3d.renew_coordinate_file(fp1 + "trj/8blk_200.lammpstrj.2000000", v=True)

    # print(len(cgu_atom3d.Atoms))
    # cgu_atom3d.convert_to_default()
    # out_path = "/home/centos/ztz/stress/small_split/8blk_200/1E-6_x/"
    # parser.create_data_file(cgu_atom3d, out_path + "data.8blk_200", q=False, improper=True)
    # parser.write_coordinate_file(cgu_atom3d, out_path + "init.lammpstrj.0")
    #
    # out_path = "/home/centos/ztz/stress/small_split/8blk_200/1E-6_y/"
    # parser.create_data_file(cgu_atom3d, out_path + "data.8blk_200", q=False, improper=True)
    # parser.write_coordinate_file(cgu_atom3d, out_path + "init.lammpstrj.0", perm=[1, 2, 0])
    #
    # out_path = "/home/centos/ztz/stress/small_split/8blk_200/1E-6_z/"
    # parser.create_data_file(cgu_atom3d, out_path + "data.8blk_200", q=False, improper=True)
    # parser.write_coordinate_file(cgu_atom3d, out_path + "init.lammpstrj.0", perm=[2, 0, 1])


    # print(len(cgu_atom3d.Atoms))
    # split_all_mole(cgu_atom3d)
    # cgu_atom3d.convert_to_default()
    # out_path = "/home/centos/ztz/stress/small_split/4blk_400/1E-6_x/"
    # parser.create_data_file(cgu_atom3d, out_path + "data.4blk_400", q=False, improper=True)
    # parser.write_coordinate_file(cgu_atom3d, out_path + "init.lammpstrj.0")
    #
    # out_path = "/home/centos/ztz/stress/small_split/4blk_400/1E-6_y/"
    # parser.create_data_file(cgu_atom3d, out_path + "data.4blk_400", q=False, improper=True)
    # parser.write_coordinate_file(cgu_atom3d, out_path + "init.lammpstrj.0", perm=[1, 2, 0])
    #
    # out_path = "/home/centos/ztz/stress/small_split/4blk_400/1E-6_z/"
    # parser.create_data_file(cgu_atom3d, out_path + "data.4blk_400", q=False, improper=True)
    # parser.write_coordinate_file(cgu_atom3d, out_path + "init.lammpstrj.0", perm=[2, 0, 1])

    # print(len(cgu_atom3d.Atoms))
    # split_all_mole(cgu_atom3d)
    # split_all_mole(cgu_atom3d)
    # cgu_atom3d.convert_to_default()
    # out_path = "/home/centos/ztz/stress/small_split/2blk_800/1E-6_x/"
    # parser.create_data_file(cgu_atom3d, out_path + "data.2blk_800", q=False, improper=True)
    # parser.write_coordinate_file(cgu_atom3d, out_path + "init.lammpstrj.0")
    #
    # out_path = "/home/centos/ztz/stress/small_split/2blk_800/1E-6_y/"
    # parser.create_data_file(cgu_atom3d, out_path + "data.2blk_800", q=False, improper=True)
    # parser.write_coordinate_file(cgu_atom3d, out_path + "init.lammpstrj.0", perm=[1, 2, 0])
    #
    # out_path = "/home/centos/ztz/stress/small_split/2blk_800/1E-6_z/"
    # parser.create_data_file(cgu_atom3d, out_path + "data.2blk_800", q=False, improper=True)
    # parser.write_coordinate_file(cgu_atom3d, out_path + "init.lammpstrj.0", perm=[2, 0, 1])

    print(len(cgu_atom3d.Atoms))
    split_all_mole(cgu_atom3d)
    split_all_mole(cgu_atom3d)
    split_all_mole(cgu_atom3d)
    cgu_atom3d.convert_to_default()
    out_path = "/home/centos/ztz/stress/small_split/1blk_1600/1E-6_x/"
    parser.create_data_file(cgu_atom3d, out_path + "data.1blk_1600", q=False, improper=True)
    parser.write_coordinate_file(cgu_atom3d, out_path + "init.lammpstrj.0")

    out_path = "/home/centos/ztz/stress/small_split/1blk_1600/1E-6_y/"
    parser.create_data_file(cgu_atom3d, out_path + "data.1blk_1600", q=False, improper=True)
    parser.write_coordinate_file(cgu_atom3d, out_path + "init.lammpstrj.0", perm=[1, 2, 0])

    out_path = "/home/centos/ztz/stress/small_split/1blk_1600/1E-6_z/"
    parser.create_data_file(cgu_atom3d, out_path + "data.1blk_1600", q=False, improper=True)
    parser.write_coordinate_file(cgu_atom3d, out_path + "init.lammpstrj.0", perm=[2, 0, 1])