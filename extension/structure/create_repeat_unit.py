from pyMD import collective_structure_class as csc
from pyMD import basic_structure_class as bsc
import numpy as np
import pickle as pkl

if __name__ == '__main__':
    ##### create Ph
    # Ph = csc.RepeatUnit()
    # a1 = bsc.Atom(1, "c2", [0, 0, 0])
    # a2 = bsc.Atom(2, "c2", [1.39, 0, 0])
    # a3 = bsc.Atom(3, "c2", [2.085, 1.203775, 0])
    # a4 = bsc.Atom(4, "c2", [1.39, 2.40755, 0])
    # a5 = bsc.Atom(5, "c2", [0, 2.40755, 0])
    # a6 = bsc.Atom(6, "c2", [-0.695, 1.203775, 0])
    # Ph.append_element(a1)
    # Ph.append_element(a2)
    # Ph.append_element(a3)
    # Ph.append_element(a4)
    # Ph.append_element(a5)
    # Ph.append_element(a6)
    # b1 = bsc.Bond(1, 1, 2, _type=("c2", "c2"))
    # b2 = bsc.Bond(2, 2, 3, _type=("c2", "c2"))
    # b3 = bsc.Bond(3, 3, 4, _type=("c2", "c2"))
    # b4 = bsc.Bond(4, 4, 5, _type=("c2", "c2"))
    # b5 = bsc.Bond(5, 5, 6, _type=("c2", "c2"))
    # b6 = bsc.Bond(6, 1, 6, _type=("c2", "c2"))
    # Ph.append_element(b1)
    # Ph.append_element(b2)
    # Ph.append_element(b3)
    # Ph.append_element(b4)
    # Ph.append_element(b5)
    # Ph.append_element(b6)
    # Ph.head = 1
    # Ph.end = 4
    # Ph.set_centroid(np.array([0., 0., 0.]))
    #
    # with open("../template/Ph.pkl", "wb") as f:
    #     pkl.dump(Ph, f)

    ## create U
    # U = csc.RepeatUnit()
    # a1 = bsc.Atom(1, "n1", [-1.17779, -0.68, 0.])
    # a2 = bsc.Atom(2, "n1", [1.17779, -0.68, 0.])
    # a3 = bsc.Atom(3, "c3", [0., 0., 0.])
    # a4 = bsc.Atom(4, "o1", [0., 1.22, 0.])
    # U.append_element(a1)
    # U.append_element(a2)
    # U.append_element(a3)
    # U.append_element(a4)
    # b1 = bsc.Bond(1, 3, 1, _type=("c3", "n1"))
    # b2 = bsc.Bond(2, 3, 2, _type=("c3", "n1"))
    # b3 = bsc.Bond(3, 3, 4, _type=("c3", "o1"))
    # U.append_element(b1)
    # U.append_element(b2)
    # U.append_element(b3)
    # U.head = 1
    # U.end = 2
    # U.set_centroid(np.array([0., 0., 0.]))
    #
    # with open("../template/U.pkl", "wb") as f:
    #     pkl.dump(U, f)

    # Uh = csc.RepeatUnit()
    # a1 = bsc.Atom(1, "n1", [-1.17779, -0.68, 0.])
    # a2 = bsc.Atom(2, "n1", [1.17779, -0.68, 0.])
    # a3 = bsc.Atom(3, "c3", [0., 0., 0.])
    # a4 = bsc.Atom(4, "o1", [0., 1.22, 0.])
    # a5 = bsc.Atom(5, "h1", [-1.17779, -1.68, 0.])
    # a6 = bsc.Atom(6, "h1", [1.17779, -1.68, 0.])
    # Uh.append_element(a1)
    # Uh.append_element(a2)
    # Uh.append_element(a3)
    # Uh.append_element(a4)
    # Uh.append_element(a5)
    # Uh.append_element(a6)
    # b1 = bsc.Bond(1, 3, 1, _type=("c3", "n1"))
    # b2 = bsc.Bond(2, 3, 2, _type=("c3", "n1"))
    # b3 = bsc.Bond(3, 3, 4, _type=("c3", "o1"))
    # b4 = bsc.Bond(4, 5, 1, _type=("h1", "n1"))
    # b5 = bsc.Bond(5, 6, 2, _type=("h1", "n1"))
    #
    # Uh.append_element(b1)
    # Uh.append_element(b2)
    # Uh.append_element(b3)
    # Uh.append_element(b4)
    # Uh.append_element(b5)
    # Uh.head = 1
    # Uh.end = 2
    # Uh.set_centroid(np.array([0., 0., 0.]))
    #
    # with open("../template/Uh.pkl", "wb") as f:
    #     pkl.dump(Uh, f)


    ### create B

    # with open("../template/Ph.pkl", "rb") as f:
    #     Ph = pkl.load(f)
    #
    # atom3d = csc.Atom3D([0,0,0], [50,50,50])
    # Ph.set_centroid([25,25,25])
    # Ph.rotate_to_vector([0, 0, 10])
    #
    # for atom in Ph.Atoms.values():
    #     atom3d.append_element(atom)
    # for bond in Ph.Bonds.values():
    #     atom3d.append_element(bond)
    # from pyMD.file_parser import LmpParser
    #
    # LmpParser.create_data_file(atom3d, "ph.data")

    Uh = csc.RepeatUnit()
    a1 = bsc.Atom(1, "n1", [-1.17779, -0.68, 0.])
    a2 = bsc.Atom(2, "n1", [1.17779, -0.68, 0.])
    a3 = bsc.Atom(3, "c1", [0., 0., 0.])
    a4 = bsc.Atom(4, "o1", [0., 1.22, 0.])
    a5 = bsc.Atom(5, "h1", [-1.17779, -1.68, 0.])
    a6 = bsc.Atom(6, "h1", [1.17779, -1.68, 0.])
    Uh.append_element(a1)
    Uh.append_element(a2)
    Uh.append_element(a3)
    Uh.append_element(a4)
    Uh.append_element(a5)
    Uh.append_element(a6)
    b1 = bsc.Bond(1, 3, 1, _type=("c1", "n1"))
    b2 = bsc.Bond(2, 3, 2, _type=("c1", "n1"))
    b3 = bsc.Bond(3, 3, 4, _type=("c1", "o1"))
    b4 = bsc.Bond(4, 5, 1, _type=("h1", "n1"))
    b5 = bsc.Bond(5, 6, 2, _type=("h1", "n1"))

    Uh.append_element(b1)
    Uh.append_element(b2)
    Uh.append_element(b3)
    Uh.append_element(b4)
    Uh.append_element(b5)
    Uh.head = 1
    Uh.end = 2
    Uh.set_centroid(np.array([0., 0., 0.]))

    with open("../../template/U1.pkl", "wb") as f:
        pkl.dump(Uh, f)