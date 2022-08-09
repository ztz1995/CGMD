from pyMD import collective_structure_class as csc
from pyMD import basic_structure_class as bsc
from pyMD import functions as f
from pyMD.file_parser import LmpParser
import numpy as np


if __name__ == '__main__':

    b = 1.4
    b_each_m = 5
    gen_m = 7
    x1 = np.asarray([0., 0., 0.])
    x2 = np.asarray([10., 5., 0.])
    co_list = f.create_fix_chain(x1, x2, b, b_each_m * gen_m)

    atom3d = csc.Atom3D([0.,0.,0.], [10.,10.,10.])
    atom_id = 1
    bond_id = 1
    atom3d.append_element(bsc.Atom(atom_id, "Es", x1))
    atom_id += 1
    for co in co_list:
        atom3d.append_element(bsc.Atom(atom_id, "c1", co))
        atom3d.append_element(bsc.Bond(bond_id, atom_id, atom_id-1))
        atom_id += 1
        bond_id += 1

    atom3d.append_element(bsc.Atom(atom_id, "Es", x2))
    atom3d.append_element(bsc.Bond(bond_id, atom_id, atom_id - 1))
    LmpParser.create_data_file(atom3d, "fix_chain_test.data", q=False, improper=False)