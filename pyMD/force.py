from pyMD import basic_structure_class as bsc
from pyMD import functions as func
from pyMD import cartesian_operation as cso
from pyMD.file_parser import LmpParser, get_date
from pyMD import collective_structure_class as csc


class PairCalculator:
    def __init__(self, cut_off=12.):
        self.cut_off = cut_off

    def cal_force(self, dist):
        return 0

    def cal_energy(self, dist):
        return 0


class PairCalAA(PairCalculator):
    def __init__(self, cut_off=12., *args, **kwargs):
        super().__init__(cut_off)
        self.type_tuple = kwargs["tt"]

    def cal_energy(self, dist):
        eng = func.cal_aa_force(dist, self.type_tuple, cut_off=self.cut_off, q=True, dielectric=1.0, dsf=True)
        return eng

    def cal_force(self, dist):
        eng2 = self.cal_energy(dist + 0.0001)
        eng1 = self.cal_energy(dist - 0.0001)
        return (eng1 - eng2) / 0.0002


class Pair:
    style_map = {"aa": PairCalAA}

    def __init__(self):
        self.type_calculator = dict()

    def new_type(self, style, type_tuple, cut_off=12., *args, **kwargs):
        new_cal = self.style_map[style](cut_off, *args, **kwargs)
        self.type_calculator[type_tuple] = new_cal

    def cal_force(self, atom1: bsc.Atom, atom2: bsc.Atom):
        type_tuple = tuple(sorted([atom1.type, atom2.type]))
        if type_tuple not in self.type_calculator:
            raise KeyError("Pair type not defined!")
        dist = cso.calculate_distance(atom1.coordinate, atom2.coordinate)
        return self.type_calculator[type_tuple].cal_force(dist)

    def cal_force_tt(self, dist, tt):
        if tt not in self.type_calculator:
            raise KeyError("Pair type not defined!")
        return self.type_calculator[tt].cal_force(dist)

    def cal_energy(self, atom1: bsc.Atom, atom2: bsc.Atom):
        type_tuple = tuple(sorted([atom1.type, atom2.type]))
        if type_tuple not in self.type_calculator:
            raise KeyError("Pair type not defined!")
        dist = cso.calculate_distance(atom1.coordinate, atom2.coordinate)
        return self.type_calculator[type_tuple].cal_energy(dist)

    def cal_energy_tt(self, dist, tt):
        if tt not in self.type_calculator:
            raise KeyError("Pair type not defined!")
        return self.type_calculator[tt].cal_energy(dist)


def urea_pair():
    pair = Pair()
    aa_types = ["c1", "h1", "o1", "n1"]
    for i in range(len(aa_types)):
        for j in range(i, len(aa_types)):
            type_tuple = tuple(sorted([aa_types[i], aa_types[j]]))
            pair.new_type("aa", type_tuple, tt=type_tuple)
    return pair


class DFF:
    def __init__(self):
        self.types = list()
        self.type_mass = dict()
        self.vdw = dict()
        self.Bond = dict()
        self.Angle = dict()
        self.BondBond = dict()
        self.BondAngle = dict()
        self.Dihedral = dict()
        self.Improper = dict()
        self.bond_inc = dict()

    def load_from_data(self, data_path):
        with open(data_path, "r") as f:
            lines = f.readlines()
        # get line index of different part
        index_dict = dict()
        last_flag = None
        for idx, line in enumerate(lines):
            args = line.split()
            if args:
                if args[0][0].isupper():
                    index_dict[line.strip("\n")] = [idx]
                    if last_flag is not None:
                        index_dict[last_flag].append(idx - 1)
                    last_flag = line.strip("\n")
        index_dict[last_flag].append(len(lines) - 1)

        # types and type_mass
        num_type = dict()
        for line in lines[index_dict["Masses"][0] + 1: index_dict["Masses"][1] + 1]:
            args = line.split()
            if args:
                num = int(args[0])
                mass = float(args[1])
                t_name = args[-1]
                num_type[num] = t_name
                self.types.append(t_name)
                self.type_mass[t_name] = mass

        # Atoms
        atoms = dict()
        for line in lines[index_dict["Atoms"][0] + 1:index_dict["Atoms"][1] + 1]:
            args = line.split()
            if args:
                a_id = int(args[0])
                m_id = int(args[1])
                a_type = num_type[int(args[2])]
                q = float(args[3])
                atoms[a_id] = [m_id, a_type, q]

        # Bonds
        bset = set()
        num_bond = dict()
        for line in lines[index_dict["Bonds"][0] + 1:index_dict["Bonds"][1] + 1]:
            args = line.split()
            if args:
                if args[1] in bset:
                    # t1 = atoms[int(args[2])][1]
                    # t2 = atoms[int(args[3])][1]
                    # if num_bond[int(args[1])] != (t1, t2):
                    #     print(num_bond[int(args[1])], (t1, t2))
                    continue
                else:
                    t1 = atoms[int(args[2])][1]
                    t2 = atoms[int(args[3])][1]
                    num_bond[int(args[1])] = (t1, t2)
                    bset.add(args[1])

        # Angles
        aset = set()
        num_angle = dict()
        for line in lines[index_dict["Angles"][0] + 1:index_dict["Angles"][1] + 1]:
            args = line.split()
            if args:
                if args[1] in aset:
                    # t1 = atoms[int(args[2])][1]
                    # t2 = atoms[int(args[3])][1]
                    # t3 = atoms[int(args[4])][1]
                    #
                    # if num_angle[int(args[1])] != (t1, t2, t3):
                    #     print(num_angle[int(args[1])], (t1, t2, t3))
                    continue
                else:
                    t1 = atoms[int(args[2])][1]
                    t2 = atoms[int(args[3])][1]
                    t3 = atoms[int(args[4])][1]
                    num_angle[int(args[1])] = (t1, t2, t3)
                    aset.add(args[1])
        # print(num_angle)

        # Dihedrals
        dset = set()
        num_dihedral = dict()
        for line in lines[index_dict["Dihedrals"][0] + 1:index_dict["Dihedrals"][1] + 1]:
            args = line.split()
            if args:
                if args[1] in dset:
                    continue
                else:
                    t1 = atoms[int(args[2])][1]
                    t2 = atoms[int(args[3])][1]
                    t3 = atoms[int(args[4])][1]
                    t4 = atoms[int(args[5])][1]
                    num_dihedral[int(args[1])] = (t1, t2, t3, t4)
                    dset.add(args[1])

        # Impropers
        iset = set()
        num_improper = dict()
        for line in lines[index_dict["Impropers"][0] + 1:index_dict["Impropers"][1] + 1]:
            args = line.split()
            if args:
                if args[1] in iset:
                    continue
                else:
                    t1 = atoms[int(args[2])][1]
                    t2 = atoms[int(args[3])][1]
                    t3 = atoms[int(args[4])][1]
                    t4 = atoms[int(args[5])][1]
                    num_improper[int(args[1])] = (t1, t2, t3, t4)
                    iset.add(args[1])

        # Pair Coeffs
        for line in lines[index_dict["Pair Coeffs"][0] + 1:index_dict["Pair Coeffs"][1] + 1]:
            args = line.split()
            if args:
                self.vdw[num_type[int(args[0])]] = [float(args[1]), float(args[2])]

        # Bond Coeffs
        for line in lines[index_dict["Bond Coeffs"][0] + 1:index_dict["Bond Coeffs"][1] + 1]:
            args = line.split()
            if args:
                self.Bond[num_bond[int(args[0])]] = [float(args[1]), float(args[2]), float(args[3]), float(args[4])]

        # Angle Coeffs
        for line in lines[index_dict["Angle Coeffs"][0] + 1:index_dict["Angle Coeffs"][1] + 1]:
            args = line.split()
            if args:
                self.Angle[num_angle[int(args[0])]] = [float(args[1]), float(args[2]), float(args[3]), float(args[4])]

        # BondBond Coeffs
        for line in lines[index_dict["BondBond Coeffs"][0] + 1:index_dict["BondBond Coeffs"][1] + 1]:
            args = line.split()
            if args:
                self.BondBond[num_angle[int(args[0])]] = [float(args[1]), float(args[2]), float(args[3])]

        # BondAngle Coeffs
        for line in lines[index_dict["BondAngle Coeffs"][0] + 1:index_dict["BondAngle Coeffs"][1] + 1]:
            args = line.split()
            if args:
                self.BondAngle[num_angle[int(args[0])]] = [float(args[1]), float(args[2]), float(args[3]),
                                                           float(args[4])]

        # Dihedral Coeffs
        for line in lines[index_dict["Dihedral Coeffs"][0] + 1:index_dict["Dihedral Coeffs"][1] + 1]:
            args = line.split()
            if args:
                self.Dihedral[num_dihedral[int(args[0])]] = [float(args[1]), float(args[2]), float(args[3]),
                                                             float(args[4])]

        # Improper Coeffs
        for line in lines[index_dict["Improper Coeffs"][0] + 1:index_dict["Improper Coeffs"][1] + 1]:
            args = line.split()
            if args:
                self.Improper[num_improper[int(args[0])]] = [float(args[1]), float(args[2])]

    def cal_charge(self, atom3d: csc.Atom3D):
        for atom in atom3d.Atoms.values():
            atom.charge = 0
            for bond_id in atom.Bonds:
                bond = atom3d.Bonds[bond_id]
                s = 1 if atom.id == bond.atoms[0] else -1
                if bond.type in self.bond_inc:
                    charge = self.bond_inc[bond.type]
                elif bond.type[::-1] in self.bond_inc:
                    charge = - self.bond_inc[bond.type[::-1]]
                else:
                    print(bond.type)
                    raise KeyError("no bond inc defined")
                atom.charge += charge * s

    def output_data(self, atom3d: csc.Atom3D, outpath):
        self.cal_charge(atom3d)
        atom3d.clear_all_joints()
        print(self.Improper)
        atom3d.cal_all_joints(improper="DFF", improper_dict=self.Improper)

        atom_list = list(atom3d.Atom_type_dict.keys())
        atom_list.sort()
        atom_type_dict = {_type: idx + 1 for idx, _type in enumerate(atom_list)}
        joint_type_dict = dict()
        for joint_n in ["Bond", "Angle", "Dihedral", "Improper"]:
            joint_list = list(atom3d.__dict__[joint_n + "_type_dict"].keys())
            joint_list.sort()
            joint_type_dict[joint_n] = {_type: idx + 1 for idx, _type in enumerate(joint_list)}

        string_lines = "LAMMPS data file. %s\n\n" % get_date()
        string_lines += "\t%d atoms\n" % len(atom3d.Atoms)
        string_lines += "\t%d bonds\n" % len(atom3d.Bonds)
        string_lines += "\t%d angles\n" % len(atom3d.Angles)
        string_lines += "\t%d dihedrals\n" % len(atom3d.Dihedrals)
        string_lines += "\t%d impropers\n\n" % len(atom3d.Impropers)
        string_lines += "\t%d atom types\n" % len(atom_type_dict)
        string_lines += "\t%d bond types\n" % len(joint_type_dict["Bond"])
        string_lines += "\t%d angle types\n" % len(joint_type_dict["Angle"])
        string_lines += "\t%d dihedral types\n" % len(joint_type_dict["Dihedral"])
        string_lines += "\t%d improper types\n\n" % len(joint_type_dict["Improper"])
        string_lines += "\t%.9f    %.9f xlo xhi\n" % (atom3d.box_l[0], atom3d.box_h[0])
        string_lines += "\t%.9f    %.9f ylo yhi\n" % (atom3d.box_l[1], atom3d.box_h[1])
        string_lines += "\t%.9f    %.9f zlo zhi\n\n" % (atom3d.box_l[2], atom3d.box_h[2])
        string_lines += "Masses\n\n"
        for idx, _type in enumerate(atom_list):
            string_lines += "{:>7}    {:.6f}".format(idx + 1, self.type_mass[_type])
            string_lines += "\t\t# %s\n" % _type
        string_lines += "\nAtoms\n\n"
        for atom_id, atom in atom3d.Atoms.items():
            atom_line = str(atom_id).rjust(12) + str(atom.mole_id).rjust(12) + str(
                atom_type_dict[atom.type]).rjust(12)
            atom_line += ("%.4f" % atom.charge).rjust(12)
            atom_line += "  " + "".join([("%.8f" % _).rjust(15) for _ in atom.coordinate]) + "\n"
            string_lines += atom_line
        for joint_n in ["Bond", "Angle", "Dihedral", "Improper"]:
            string_lines += "\n%ss\n\n" % joint_n
            for joint_id, joint in atom3d.__dict__[joint_n + "s"].items():
                if joint.type in self.__dict__[joint_n]:
                    jatoms = joint.atoms
                elif joint.type[::-1] in self.__dict__[joint_n]:
                    jatoms = joint.atoms[::-1]
                else:
                    print(joint_n, joint.type)
                    raise KeyError("No joint type defined")
                string_lines += str(joint_id).rjust(12) + str(joint_type_dict[joint_n][joint.type]).rjust(12) \
                                + "".join([str(_).rjust(12) for _ in jatoms]) + "\n"

        string_lines += "\nPair Coeffs\n\n"
        for at in atom_type_dict:
            string_lines += str(atom_type_dict[at]).rjust(12) + "".join(
                [("%.4f" % _).rjust(12) for _ in self.vdw[at]]) + "\n"

        for joint_n in ["Bond", "Angle", "Dihedral", "Improper"]:
            string_lines += "\n%s Coeffs\n\n" % joint_n
            for jt in joint_type_dict[joint_n]:
                jjt = jt if jt in self.__dict__[joint_n] else jt[::-1]
                coeff = self.__dict__[joint_n][jjt]
                string_lines += str(joint_type_dict[joint_n][jt]).rjust(12) + "".join(
                    [("%.4f" % _).rjust(12) for _ in coeff]) + "\n"
        # print(self.BondAngle)
        # print(self.BondBond)

        for co_n in ["BondBond", "BondAngle"]:
            string_lines += "\n%s Coeffs\n\n" % co_n
            for jt in joint_type_dict["Angle"]:
                jjt = jt if jt in self.__dict__[co_n] else jt[::-1]
                coeff = self.__dict__[co_n][jjt]
                string_lines += str(joint_type_dict["Angle"][jt]).rjust(12) + "".join(
                    [("%.4f" % _).rjust(12) for _ in coeff]) + "\n"

        with open(outpath, "w") as f:
            f.write(string_lines)
        pass


if __name__ == '__main__':
    # upair = urea_pair()
    # a1 = bsc.Atom(1, "h1", [0,0,1])
    # a2 = bsc.Atom(2, "o1", [0,0,3])
    # f = upair.cal_force(a1, a2)
    # print(f)

    dpath = "../UPy/UPy.data"
    ff = DFF()
    ff.load_from_data(dpath)
