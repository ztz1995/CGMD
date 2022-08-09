import numpy as np
import os
from pyMD import basic_structure_class as bsc
from pyMD import collective_structure_class as csc
from pyMD import ibi
from pyMD.file_parser import LmpParser


def create_hyperAn(joint, end, arm, n_arm=3, generation=3):
    def growth(_gl):
        new_gl = list()
        for g in _gl:
            new_gl.append(add_arm(g))
        return new_gl

    def add_arm(start_id):
        nonlocal a_id, b_id
        for j in range(len(arm)):
            new_a = bsc.Atom(a_id, _type=arm[j], mole_id=1)
            hyper_mole.append_element(new_a)
            l_id = start_id if j == 0 else a_id - 1
            _type, rev = csc.sort_types((hyper_mole.Atoms[l_id].type, arm[j]))
            if rev:
                new_b = bsc.Bond(b_id, a_id, l_id, _type=_type)
            else:
                new_b = bsc.Bond(b_id, l_id, a_id, _type=_type)
            hyper_mole.append_element(new_b)
            a_id += 1
            b_id += 1
        return a_id - 1

    def add_joint(_gl):
        nonlocal a_id, b_id
        new_growth = list()
        for g in _gl:
            new_center = bsc.Atom(a_id, _type=joint, mole_id=1)
            _type, rev = csc.sort_types((hyper_mole.Atoms[g].type, joint))
            if rev:
                new_b = bsc.Bond(b_id, a_id, g, _type=_type)
            else:
                new_b = bsc.Bond(b_id, g, a_id, _type=_type)
            hyper_mole.append_element(new_center)
            hyper_mole.append_element(new_b)
            new_growth += [a_id] * (n_arm - 1)
            a_id += 1
            b_id += 1
        return new_growth

    def add_end(_gl):
        nonlocal a_id, b_id
        for g in _gl:
            new_center = bsc.Atom(a_id, _type=end, mole_id=1)
            _type, rev = csc.sort_types((hyper_mole.Atoms[g].type, end))
            if rev:
                new_b = bsc.Bond(b_id, a_id, g, _type=_type)
            else:
                new_b = bsc.Bond(b_id, g, a_id, _type=_type)
            hyper_mole.append_element(new_center)
            hyper_mole.append_element(new_b)
            a_id += 1
            b_id += 1

    hyper_mole = csc.Molecule()
    a_id = 1
    b_id = 1
    center = bsc.Atom(a_id, _type=joint, mole_id=1)
    hyper_mole.append_element(center)
    hyper_mole.center_id = a_id
    growth_list = [a_id] * n_arm
    a_id += 1
    for i in range(generation):
        growth_list = growth(growth_list)
        if i < generation - 1:
            growth_list = add_joint(growth_list)
        elif i == generation - 1:
            add_end(growth_list)

    hyper_mole.cal_all_joints(improper="ignore")
    return hyper_mole


def main():
    joint = "TO(2)"
    end = "TO(1)"
    s_len = 7
    arm = ["TO(2)"] * s_len + ["Es", "Ph", "U", "Ph", "Me", "Ph", "U", "Ph", "Es"] + ["TO(2)"] * s_len
    n_arm = 3
    generation = 5
    hyper_mole = create_hyperAn(joint, end, arm, n_arm, generation)

    fp1 = "../../data/20200720_lr=0.05_1blk_50_cg_ex_hard/record/"
    r1 = ibi.Record(fp1, restart=True)
    param1 = r1.get_param(52, 1)

    # hyper_mole.random_coordinate(1, np.zeros(3), param1)
    # print(hyper_mole)

    atom3d = csc.Atom3D()
    mole_num = 10
    atom3d.create_from_molecule(hyper_mole, mole_num, param=param1, density=1., seed=1)
    # print(atom3d)
    fp = "/home/centos/work/hyper/S%dHS%d/S%d_G%d_%d/init/" % (s_len, s_len, n_arm, generation, mole_num)
    os.makedirs(fp, exist_ok=True)
    LmpParser.create_data_file(atom3d, fp + "hyper.data", q=False, improper=False)


#


if __name__ == '__main__':
    main()
