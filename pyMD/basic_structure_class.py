import numpy as np
import copy as cp


class Atom:
    def __init__(self, _id, _type=None, _coordinate=None, _image=None, aa_ids=None, mole_id=0, charge=0., _mv2_aa=None, _mv2_cg=None, v=None):
        self.id = _id
        self.aa_ids = [] if aa_ids is None else aa_ids
        self.coordinate = np.zeros(3) if _coordinate is None else np.asarray(_coordinate, dtype=float)
        self.mv2_aa = 0 if _mv2_aa is None else _mv2_aa
        self.mv2_cg = 0 if _mv2_cg is None else _mv2_cg
        self.type = _type
        if _image is None:
            self.image = np.zeros(3, dtype=int)
        else:
            assert len(_image) == 3
            self.image = np.asarray(_image, dtype=int)
        self.Bonds = list()
        self.Angles = list()
        self.Dihedrals = list()
        self.Impropers = list()
        self.mole_id = mole_id
        self.charge = charge
        self.v = np.zeros(3) if v is None else v
        self.f = np.zeros(3)

    def __str__(self):
        np.set_printoptions(formatter={'float': lambda x: ('%.3f' % x).rjust(8)})
        return "Atom  ID: " + str(
            self.id) + ",\tType: " + self.type + ",\tCoordinate: " + str(
            self.coordinate) + ",\t\tPeriod: " + str(self.image) + ",\t\tMole_id: " + str(self.mole_id)

    def get_image_coordinate(self, image, lattice_parameter):
        co = cp.copy(self.coordinate)
        co -= (np.asarray(image) - np.asarray(self.image)) * lattice_parameter
        return co

    def set_image(self, image, lattice_parameter):
        new_co = self.get_image_coordinate(image, lattice_parameter)
        self.coordinate = new_co
        self.image = image

    def set_in_cell_coordinate(self, lattice_parameter, box_l):
        co = self.coordinate
        self.image = np.asarray((co - box_l) // lattice_parameter, dtype=int) + self.image
        self.coordinate = (co - box_l) % lattice_parameter + box_l

    def set_default_coordinate(self, lattice_parameter):
        self.coordinate += lattice_parameter * self.image
        self.image = np.zeros(3, dtype=int)

    # no bond information
    def get_copy(self):
        return Atom(self.id, self.type, self.coordinate.copy(), self.image.copy(), mole_id=self.mole_id, charge=self.charge, v=self.v.copy())

    def move_along_vector(self, vec):
        self.coordinate += vec

    @property
    def mass(self):
        type_mass_dict = {"TO(1)": 73.1142, "TO(2)": 72.1062, "TO(3)": 73.1142, "Es": 72.0593, "Me": 14.0267,
                          "Ph": 76.093640, "U": 58.039890,
                          "c1": 12.0107, "o1": 15.9994, "n1": 14.0067, "h1": 1.008, "sv1": 16, "sv2": 56, "sv": 75}
        a_dict = {"c": 12.0107, "o": 15.9994, "n": 14.0067, "h": 1.008}
        if self.type in type_mass_dict:
            return type_mass_dict[self.type]
        else:
            assert self.type[0].islower()
            return a_dict[self.type[0]]


class Joint:
    def __init__(self, _id, *atoms, _type=None):
        self.id = _id
        self.type = _type
        self.atoms = list(atoms)

    def __str__(self):
        return self.__class__.__name__ + " ID: " + str(
            self.id) + "\tType: " + self.type.__repr__() + "\tAtoms: " + self.atoms.__repr__()


class Bond(Joint):

    def get_other_atom_id(self, atom_id):
        if atom_id == self.atoms[0]:
            return self.atoms[1]
        elif atom_id == self.atoms[1]:
            return self.atoms[0]
        else:
            raise KeyError("Atom id not in bond!")


class Angle(Joint):
    pass


class Dihedral(Joint):
    pass


class Improper(Joint):
    pass


if __name__ == '__main__':
    atom1 = Atom(12, _type="C")
    co = atom1.get_image_coordinate([1, 1, 1], 10)
    print(atom1)
