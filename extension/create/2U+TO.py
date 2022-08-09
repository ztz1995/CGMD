import pickle as pkl
import copy as cp
from pyMD import collective_structure_class as csc
from pyMD import basic_structure_class as bsc
from pyMD import ibi
from pyMD.file_parser import LmpParser

if __name__ == '__main__':
    with open("../../template/U1.pkl", "rb") as file:
        U = pkl.load(file)
    print(U)
    Ph_l = bsc.Atom(7, "Ph", _coordinate=[-3., 0., 0.])
    Ph_r = bsc.Atom(8, "Ph", _coordinate=[3., 0., 0.])
    b1 = bsc.Bond(6, 7, 1, _type=("Ph", "n1"))
    b2 = bsc.Bond(7, 8, 2, _type=("Ph", "n1"))
    U.append_element(Ph_l)
    U.append_element(Ph_r)
    U.append_element(b1)
    U.append_element(b2)
    U2 = cp.deepcopy(U)
    U.set_centroid([20, 17.5, 20])
    U2.set_centroid([20, 22.5, 20])

    fp1 = "../../data/20200720_lr=0.05_1blk_50_cg_ex_hard/record/"
    r1 = ibi.Record(fp1, restart=True)
    # param1 = r1.get_param(22, 1)
    param1 = r1.get_param(52, 1)
    pattern = ["TO(2)"] * 8
    pattern[0] = "TO(1)"
    pattern[-1] = "TO(1)"
    chains = 60
    atom3d = csc.Atom3D()
    atom3d.create_from_pattern(pattern, chains, param=param1, density=0.9)
    ma = atom3d.get_max_atom_id()
    mb = atom3d.get_max_bond_id()

    for atom in U.Atoms.values():
        n_atom = cp.copy(atom)
        n_atom.id += ma
        atom3d.append_element(n_atom)
    for atom in U2.Atoms.values():
        n_atom = cp.copy(atom)
        n_atom.id += 8 + ma
        atom3d.append_element(n_atom)
    for bond in U.Bonds.values():
        n_bond = cp.copy(bond)
        n_bond.id += mb
        n_bond.atoms = [_ + ma for _ in bond.atoms]
        atom3d.append_element(n_bond)
    for bond in U2.Bonds.values():
        n_bond = cp.copy(bond)
        n_bond.id += 8 + mb
        n_bond.atoms = [_ + 8 + ma for _ in bond.atoms]
        atom3d.append_element(n_bond)
    atom3d.sort_atom_bond_id()
    # print("2")
    atom3d.cal_all_joints(improper="ignore")
    LmpParser.create_data_file(atom3d, "uuTO.data", q=False, improper=True)
