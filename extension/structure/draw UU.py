import pickle as pkl
import copy as cp
from pyMD import collective_structure_class as csc
from pyMD.file_parser import LmpParser
import numpy as np


if __name__ == '__main__':
    with open("../../template/U1.pkl", "rb") as file:
        Uh = pkl.load(file)

    # print(Uh)

    # u_orien = [1., 0., 0.]
    # old_orient = np.asarray([0., 1., 0.])
    # new_orient = np.asarray(u_orien)
    # Uh.rotate_orientation(old_orient, new_orient)
    #
    # print(Uh)


    Uh2 = cp.deepcopy(Uh)

    Uh.set_centroid([5,5,5])
    Uh2.set_centroid([5,0.5,5])
    # Uh2.set_centroid([0, 0., -4.5])

    # print(Uh.__dict__)
    atom3d = csc.Atom3D(box_h=[10,10,10])
    for atom in Uh.Atoms.values():
        atom3d.append_element(atom)
    for atom in Uh2.Atoms.values():
        n_atom = cp.copy(atom)
        n_atom.id+=6
        atom3d.append_element(n_atom)
    for bond in Uh.Bonds.values():
        atom3d.append_element(bond)
    for bond in Uh2.Bonds.values():
        n_bond = cp.copy(bond)
        n_bond.id += 6
        n_bond.atoms = [_+6 for _ in bond.atoms]
        atom3d.append_element(n_bond)
    LmpParser.create_data_file(atom3d, "uu2.data", q=False, improper=False)
