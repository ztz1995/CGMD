import xml.dom.minidom as minidom
from pyMD import basic_structure_class as bsc
from pyMD import collective_structure_class as csc
import os
from collections import OrderedDict
import time
import logging
import sys
import matplotlib.pyplot as plt
import numpy as np


def create_log(file_path, restart):
    # def __init__(self, file_path, restart):
    logger = logging.getLogger('mylogger')
    logger.setLevel(logging.DEBUG)
    sc_handler = logging.StreamHandler(stream=sys.stdout)
    sc_handler.setLevel(logging.DEBUG)
    sc_handler.setFormatter(
        logging.Formatter("%(asctime)s - %(filename)s[:%(lineno)d] - %(message)s"))

    if (not restart) and os.path.exists(file_path + 'run.log'):
        os.remove(file_path + 'run.log')
    f_handler = logging.FileHandler(file_path + 'run.log')
    f_handler.setLevel(logging.DEBUG)
    f_handler.setFormatter(
        logging.Formatter("%(asctime)s - %(levelname)s - %(filename)s[:%(lineno)d] - %(message)s"))

    logger.addHandler(sc_handler)
    logger.addHandler(f_handler)
    return logger

    # def debug(self, message):
    #     self.logger.debug(message)
    #
    # def info(self, message):
    #     self.logger.info(message)


def sort_key(old_dict, reverse=False):
    keys = sorted(old_dict.keys(), reverse=reverse)
    new_dict = OrderedDict()
    for key in keys:
        new_dict[key] = old_dict[key]
    return new_dict


def get_date():
    return time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))


class MsParser:
    def __init__(self):
        pass

    @staticmethod
    def load_file(file_path, type_dict=None):
        if type_dict is None:
            type_dict = {}

        def input_xml():
            with open(file_path, 'r', encoding='latin1') as my_file:
                file_string = my_file.read()
            new_string = file_string.replace("\n", "").replace("\t", "")
            xml_minidom_document = minidom.parseString(new_string)
            return xml_minidom_document

        xml = input_xml()
        root_node = xml.documentElement
        lattice_parameter = None
        try:
            space_group = root_node.getElementsByTagName("SpaceGroup")
            for space in space_group:
                a_vector = space.getAttribute("AVector")
                lattice_parameter = a_vector.split(",")[0]
                lattice_parameter = np.asarray([float(lattice_parameter)] * 3)
        finally:
            if lattice_parameter is not None:
                extract_structure = csc.Atom3D(box_h=lattice_parameter)
            else:
                extract_structure = csc.Collection()
        all_atom3d = root_node.getElementsByTagName("Atom3d")
        for aa in all_atom3d:
            if aa.hasAttribute("XYZ"):
                id_num = int(aa.getAttribute("ID"))
                xyz = aa.getAttribute("XYZ")
                if lattice_parameter is not None:
                    split_xyz = np.asarray([float(i) for i in xyz.split(",")]) * lattice_parameter
                else:
                    split_xyz = np.asarray([float(i) for i in xyz.split(",")])
                frc_type = aa.getAttribute("ForcefieldType")
                if frc_type in type_dict:
                    frc_type = type_dict[frc_type]
                extract_structure.append_element(bsc.Atom(id_num, frc_type, split_xyz))
        # print(extract_structure.Atoms)
        all_bond = root_node.getElementsByTagName("Bond")
        for bond in all_bond:
            if bond.hasAttribute("Connects"):
                id_num = int(bond.getAttribute("ID"))
                connects = bond.getAttribute("Connects").split(",")
                a1 = extract_structure.Atoms[int(connects[0])]
                a2 = extract_structure.Atoms[int(connects[1])]
                bond_type = tuple(sorted([a1.type, a2.type]))
                if bond_type == (a1.type, a2.type):
                    atoms = [a1.id, a2.id]
                else:
                    atoms = [a2.id, a1.id]
                extract_structure.append_element(bsc.Bond(id_num, *atoms, _type=bond_type))
        extract_structure.clear_all_joints()
        extract_structure.sort_atom_bond_id()
        return extract_structure


class LmpParser:
    def __init__(self):
        pass

    @staticmethod
    def create_data_file(atom3d: csc.Atom3D, file_path, q=False, improper=True, dihedral=True, velocity=False):
        if q:
            atom3d.cal_charge()
        # atom_list = list(atom3d.type_mass_dict.keys())
        atom_list = list(atom3d.Atom_type_dict.keys())
        atom_list.sort()
        atom_type_dict = {_type: idx + 1 for idx, _type in enumerate(atom_list)}
        joint_type_dict = dict()
        for joint_n in ["Bond", "Angle", "Dihedral", "Improper"]:
            joint_list = list(atom3d.__dict__[joint_n + "_type_dict"].keys())
            joint_list.sort()
            joint_type_dict[joint_n] = {_type: idx + 1 for idx, _type in enumerate(joint_list)}
            # print(joint_n)
            # print(joint_type_dict)

        string_lines = "LAMMPS data file. %s\n\n" % get_date()
        string_lines += "\t%d atoms\n" % len(atom3d.Atoms)
        string_lines += "\t%d bonds\n" % len(atom3d.Bonds)
        string_lines += "\t%d angles\n" % len(atom3d.Angles)
        if dihedral:
            string_lines += "\t%d dihedrals\n" % len(atom3d.Dihedrals)
        if improper:
            string_lines += "\t%d impropers\n\n" % len(atom3d.Impropers)
        string_lines += "\t%d atom types\n" % len(atom_type_dict)
        string_lines += "\t%d bond types\n" % len(joint_type_dict["Bond"])
        string_lines += "\t%d angle types\n" % len(joint_type_dict["Angle"])
        if dihedral:
            string_lines += "\t%d dihedral types\n" % len(joint_type_dict["Dihedral"])
        if improper:
            string_lines += "\t%d improper types\n\n" % len(joint_type_dict["Improper"])
        string_lines += "\t%.9f    %.9f xlo xhi\n" % (atom3d.box_l[0], atom3d.box_h[0])
        string_lines += "\t%.9f    %.9f ylo yhi\n" % (atom3d.box_l[1], atom3d.box_h[1])
        string_lines += "\t%.9f    %.9f zlo zhi\n\n" % (atom3d.box_l[2], atom3d.box_h[2])
        string_lines += "Masses\n\n"
        for idx, _type in enumerate(atom_list):
            string_lines += "{:>7}    {:.6f}".format(idx + 1, atom3d.get_atom_type_mass(_type))
            string_lines += "\t\t# %s\n" % _type
        string_lines += "\nAtoms\n\n"
        for atom_id, atom in atom3d.Atoms.items():
            atom_line = str(atom_id).rjust(15) + str(atom.mole_id).rjust(15) + str(
                atom_type_dict[atom.type]).rjust(15)
            if q:
                atom_line += str(atom.charge).rjust(15)
            atom_line += "  " + " ".join([("%.8f" % _).rjust(15) for _ in atom.coordinate]) + "\n"
            string_lines += atom_line
        for joint_n in ["Bond", "Angle", "Dihedral", "Improper"]:
            if joint_n == "Dihedral" and not dihedral:
                continue
            if joint_n == "Improper" and not improper:
                continue
            string_lines += "\n%ss\n\n" % joint_n
            for joint_id, joint in atom3d.__dict__[joint_n + "s"].items():
                string_lines += str(joint_id).rjust(15) + str(joint_type_dict[joint_n][joint.type]).rjust(15) \
                                + "".join([str(_).rjust(15) for _ in joint.atoms]) + "\n"
        
        if velocity:
            string_lines += "\nVelocities\n\n"
            for atom_id, atom in atom3d.Atoms.items():
                atom_line = str(atom_id).rjust(15)
                atom_line += "  " + " ".join([("%.8f" % _).rjust(15) for _ in atom.v]) + "\n"
                string_lines += atom_line
        
        with open(file_path, "w") as f:
            f.write(string_lines)
        pass

    @staticmethod
    def create_data_file_aa(atom3d: csc.Atom3D, file_path, q=True, cal_q=True, improper=True, dihedral=True, velocity=False):
        atom3d.set_type_define("aa")
        if q and cal_q:
            atom3d.cal_charge_aa(cal_q)
        # atom_list = list(atom3d.type_mass_dict.keys())
        atom_list = list(atom3d.Atom_type_dict.keys())
        atom_list.sort()
        atom_type_dict = OrderedDict([(_type, idx + 1) for idx, _type in enumerate(atom_list)])
        # atom_type_dict = {_type: idx + 1 for idx, _type in enumerate(atom_list)}
        joint_type_dict = dict()
        for joint_n in ["Bond", "Angle", "Dihedral", "Improper"]:
            joint_list = list(atom3d.__dict__[joint_n + "_type_dict"].keys())
            joint_list.sort()
            joint_type_dict[joint_n] = OrderedDict([(_type, idx + 1) for idx, _type in enumerate(joint_list)])
            # joint_type_dict[joint_n] = {_type: idx + 1 for idx, _type in enumerate(joint_list)}
            # print(joint_n)
            # print(joint_type_dict)

        string_lines = "LAMMPS data file. %s\n\n" % get_date()
        string_lines += "\t%d atoms\n" % len(atom3d.Atoms)
        string_lines += "\t%d bonds\n" % len(atom3d.Bonds)
        string_lines += "\t%d angles\n" % len(atom3d.Angles)
        if dihedral:
            string_lines += "\t%d dihedrals\n" % len(atom3d.Dihedrals)
        if improper:
            string_lines += "\t%d impropers\n\n" % len(atom3d.Impropers)
        string_lines += "\t%d atom types\n" % len(atom_type_dict)
        string_lines += "\t%d bond types\n" % len(joint_type_dict["Bond"])
        string_lines += "\t%d angle types\n" % len(joint_type_dict["Angle"])
        if dihedral:
            string_lines += "\t%d dihedral types\n" % len(joint_type_dict["Dihedral"])
        if improper:
            string_lines += "\t%d improper types\n\n" % len(joint_type_dict["Improper"])
        string_lines += "\t%.9f    %.9f xlo xhi\n" % (atom3d.box_l[0], atom3d.box_h[0])
        string_lines += "\t%.9f    %.9f ylo yhi\n" % (atom3d.box_l[1], atom3d.box_h[1])
        string_lines += "\t%.9f    %.9f zlo zhi\n\n" % (atom3d.box_l[2], atom3d.box_h[2])
        string_lines += "Masses\n\n"
        for idx, _type in enumerate(atom_list):
            string_lines += "{:>7}    {:.6f}".format(idx + 1, atom3d.get_atom_type_mass(_type))
            string_lines += "\t\t# %s\n" % _type
        from pyMD.class2 import COMPASS
        string_lines += COMPASS.create_potential_lines(atom_type_dict, joint_type_dict)
        string_lines += "\nAtoms\n\n"
        for atom_id, atom in atom3d.Atoms.items():
            atom_line = str(atom_id).rjust(15) + str(atom.mole_id).rjust(15) + str(
                atom_type_dict[atom.type]).rjust(15)
            if q:
                atom_line += str(round(atom.charge, 5)).rjust(15)
            atom_line += "  " + " ".join([("%.8f" % _).rjust(15) for _ in atom.coordinate]) + "\n"
            string_lines += atom_line
        for joint_n in ["Bond", "Angle", "Dihedral", "Improper"]:
            if joint_n == "Dihedral" and not dihedral:
                continue
            if joint_n == "Improper" and not improper:
                continue
            string_lines += "\n%ss\n\n" % joint_n
            for joint_id, joint in atom3d.__dict__[joint_n + "s"].items():
                string_lines += str(joint_id).rjust(15) + str(joint_type_dict[joint_n][joint.type]).rjust(15) \
                                + "".join([str(_).rjust(15) for _ in joint.atoms]) + "\n"
        
        if velocity:
            string_lines += "\nVelocities\n\n"
            for atom_id, atom in atom3d.Atoms.items():
                atom_line = str(atom_id).rjust(15)
                atom_line += "  " + " ".join([("%.8f" % _).rjust(15) for _ in atom.v]) + "\n"
                string_lines += atom_line
        
        with open(file_path, "w") as f:
            f.write(string_lines)
        return atom_type_dict, joint_type_dict

    @staticmethod
    def load_data_file(file_path, _type="fg", q=False, map_dict=None, param_dict=None):
        with open(file_path, "r") as f:
            lines = f.readlines()
        index_dict = dict()
        last_flag = None
        for idx, line in enumerate(lines):
            args = line.split()
            if args:
                if (last_flag is None and args[0] == "#") or args[0] == "LAMMPS":
                    index_dict["LAMMPS"] = [idx]
                    last_flag = "LAMMPS"

                elif args[0][0].isupper():
                    index_dict[args[0]] = [idx]
                    if last_flag is not None:
                        index_dict[last_flag].append(idx - 1)
                    last_flag = args[0]
        index_dict[last_flag].append(len(lines) - 1)
        # print(index_dict)
        if param_dict is not None:
            assert isinstance(param_dict, dict)
            param_list = ['Bond', 'Angle', 'Dihedral', 'Improper', 'BondBond', 'BondAngle', 'AngleAngle', 'AngleAngleTorsion', 'EndBondTorsion',
                          'MiddleBondTorsion', 'BondBond13', 'AngleTorsion']

            for param in param_list:
                param_dict[param] = dict()
                for line in lines[index_dict[param][0] + 1: index_dict[param][1] + 1]:
                    args = line.split()
                    if args:
                        param_dict[param][int(args[0])] = ",".join(args[1:])

        if _type == "cg":
            type_list = ["TO(0)", "TO(1)", "TO(2)", "Es(0)", "Ph(0)", "U(0)", "B(0)"]
        elif _type == "dpd":
            type_list = ["Es", "Me", "Ph", "TO(1)", "TO(2)", "U"]
        elif _type == "dpd_fg":
            type_list = ["Es", "c3", "Ph", "TO(1)", "TO(2)", "U"]
        elif _type == "dpd_fg_old":
            type_list = ["Es", "c1", "Ph", "TO(1)", "TO(2)", "U"]
        elif _type == "fg_old2new":
            type_list = ["Es", "TO(1)", "TO(2)", "c3", "c2", "c1", "h1", "n1", "o1"]
        elif _type == "cg_fg":
            type_list = ["TO(0)", "TO(1)", "TO(2)", "Es", "Ph(0)", "U(0)", "B(0)"]
        elif _type == "cg_fgh":
            type_list = ["TO(0)", "TO(1)", "TO(1)", "Es", "Ph(0)", "U(0)", "B(0)"]
        elif _type == "fg":
            type_list = list()
            for line in lines[index_dict["Masses"][0] + 1: index_dict["Masses"][1] + 1]:
                args = line.split("#")
                a_type = args[-1].split()
                if a_type:
                    if map_dict:
                        if a_type[0] in map_dict:
                            a_type[0] = map_dict[a_type[0]]
                    type_list.append(a_type[0])
        else:
            type_list = None
        print(type_list)
        box_l = list()
        box_h = list()
        for line in lines[index_dict["LAMMPS"][0] + 1:index_dict["LAMMPS"][1] + 1]:
            args = line.split()
            if args:
                if args[-1][1:] == "hi":
                    box_l.append(float(args[0]))
                    box_h.append(float(args[1]))
        # print(box_h)
        atom3d = csc.Atom3D(box_l, box_h)

        for line in lines[index_dict["Atoms"][0] + 1:index_dict["Atoms"][1] + 1]:
            args = line.split()
            if args:
                a_id = int(args[0])
                m_id = int(args[1])
                a_type = args[2] if type_list is None else type_list[int(args[2]) - 1]
                if q:
                    charge = float(args[3])
                    co = [float(_) for _ in args[4:7]]
                    atom = bsc.Atom(a_id, _type=a_type, _coordinate=co, mole_id=m_id, charge=charge)
                else:
                    co = [float(_) for _ in args[3:6]]
                    atom = bsc.Atom(a_id, _type=a_type, _coordinate=co, mole_id=m_id)
                atom3d.append_element(atom)

        for line in lines[index_dict["Bonds"][0] + 1:index_dict["Bonds"][1] + 1]:
            args = line.split()
            if args:
                b_id = int(args[0])
                a1 = int(args[2])
                a2 = int(args[3])
                b_type = (atom3d.Atoms[a1].type, atom3d.Atoms[a2].type)
                bond = bsc.Bond(b_id, a1, a2, _type=b_type)
                atom3d.append_element(bond)
        for line in lines[index_dict["Angles"][0] + 1:index_dict["Angles"][1] + 1]:
            args = line.split()
            if args:
                b_id = int(args[0])
                a1 = int(args[2])
                a2 = int(args[3])
                a3 = int(args[4])
                b_type = (atom3d.Atoms[a1].type, atom3d.Atoms[a2].type, atom3d.Atoms[a3].type)
                bond = bsc.Angle(b_id, a1, a2, a3, _type=b_type)
                atom3d.append_element(bond)
        if "Dihedrals" in index_dict:
            for line in lines[index_dict["Dihedrals"][0] + 1:index_dict["Dihedrals"][1] + 1]:
                args = line.split()
                if args:
                    b_id = int(args[0])
                    a1 = int(args[2])
                    a2 = int(args[3])
                    a3 = int(args[4])
                    a4 = int(args[5])
                    b_type = (atom3d.Atoms[a1].type, atom3d.Atoms[a2].type, atom3d.Atoms[a3].type, atom3d.Atoms[a4].type)
                    bond = bsc.Dihedral(b_id, a1, a2, a3, a4, _type=b_type)
                    atom3d.append_element(bond)
        if "Impropers" in index_dict:
            for line in lines[index_dict["Impropers"][0] + 1:index_dict["Impropers"][1] + 1]:
                args = line.split()
                if args:
                    b_id = int(args[0])
                    a1 = int(args[2])
                    a2 = int(args[3])
                    a3 = int(args[4])
                    a4 = int(args[5])
                    b_type = (
                        atom3d.Atoms[a1].type, atom3d.Atoms[a2].type, atom3d.Atoms[a3].type, atom3d.Atoms[a4].type)
                    bond = bsc.Improper(b_id, a1, a2, a3, a4, _type=b_type)
                    atom3d.append_element(bond)
        return atom3d

    @staticmethod
    def renew_coordinate_file(atom3d, file_path, step, one_file=False, only_aa=False, only_u=False, only_type=None, pad=0, renew_v=False,
                              renew_f=False):
        if one_file:
            with open(file_path, "r") as f:
                lines = f.readlines()
            begin = 0
            if step is not None:
                for idx, line in enumerate(lines):
                    if line == "ITEM: TIMESTEP\n":
                        if int(lines[idx + 1].split()[0]) == step:
                            begin = idx
                            break
        else:
            with open(file_path + str(int(step)).zfill(pad), "r") as f:
                lines = f.readlines()
            begin = 0
        box_l = list()
        box_h = list()
        for line in lines[begin + 5: begin + 8]:
            arg = line.split()
            box_l.append(float(arg[0]))
            box_h.append(float(arg[1]))
        atom3d.box_l = np.array(box_l)
        atom3d.box_h = np.array(box_h)
        atom3d.lattice_parameter = atom3d.box_h - atom3d.box_l

        for line in lines[begin + 9:]:
            arg = line.split()
            if arg:
                if arg[0] == "ITEM:":
                    break
                atom_id = int(arg[0])
                if only_aa and not atom3d.Atoms[atom_id].type[0].islower():
                    continue
                if only_u and atom3d.Atoms[atom_id].type not in ["c1", "U"]:
                    continue
                if only_type and atom3d.Atoms[atom_id].type not in only_type:
                    continue
                co_s = np.asarray([float(arg[1]), float(arg[2]), float(arg[3])])
                image = np.asarray([int(arg[4]), int(arg[5]), int(arg[6])], dtype=int)
                if renew_v:
                    v = np.asarray([float(arg[7]), float(arg[8]), float(arg[9])])
                    atom3d.Atoms[atom_id].v = v
                if renew_f:
                    f = np.asarray([float(arg[10]), float(arg[11]), float(arg[12])])
                    atom3d.Atoms[atom_id].f = f
                atom3d.Atoms[atom_id].coordinate = co_s * atom3d.lattice_parameter + box_l
                atom3d.Atoms[atom_id].image = image
        atom3d.default = False
        atom3d.in_cell = False
        return atom3d

    @staticmethod
    def write_coordinate_file(atom3d, file_path, v=True, timestep=0, perm=None, pro_dict=None):
        atom3d.convert_to_in_cell()
        lines = ["ITEM: TIMESTEP\n",
                 "%d\n" % timestep,
                 "ITEM: NUMBER OF ATOMS\n",
                 "%d\n" % len(atom3d.Atoms),
                 "ITEM: BOX BOUNDS pp pp pp\n"]
        if perm is None:
            perm = [0, 1, 2]
        for i in perm:
            lines.append("%.8e %.8e\n" % (atom3d.box_l[i], atom3d.box_h[i]))
        if pro_dict is not None:
            keys = list(pro_dict.keys())
        else:
            keys = []
        lines.append("ITEM: ATOMS id xs ys zs ix iy iz vx vy vz %s\n" % " ".join(keys))
        for a_id, atom in atom3d.Atoms.items():
            rs = (atom.coordinate - atom3d.box_l) / atom3d.lattice_parameter
            image = atom.image
            v = atom.v
            l = "%d" % a_id
            for j in perm:
                l += " %.8f" % rs[j]
            for j in perm:
                l += " %d" % image[j]
            for j in perm:
                l += " %.8f" % v[j]
            for key in keys:
                l += " %.8f" % pro_dict[key][a_id-1]
            l += "\n"
            lines.append(l)
        lines.append("\n")
        with open(file_path, "w") as file:
            file.writelines(lines)
        pass


def get_project_path():
    cur_path = os.path.abspath(os.path.dirname(__file__))
    root_path = cur_path[:cur_path.find("CGMD") + len("CGMD")] + "/"
    return root_path


def invoke_lmp(in_path, in_name, num_parallel=5, screen=False):
    if screen:
        os.system("cd %s && /usr/lib64/openmpi/bin/mpirun -np %d /home/centos/bin/lmp_daily -in %s" % (in_path, num_parallel, in_name))
    else:
        os.system(
            "cd %s && /usr/lib64/openmpi/bin/mpirun -np %d /home/centos/bin/lmp_daily -in %s -screen ./screen.out" % (in_path, num_parallel, in_name))
    pass


def plot_figs(x_list, y_list, x_label, y_label, name, save_dictionary=None, label=None, axh=None, axv=None, show=False,
              xlim=None, ylim=None):
    assert len(x_list) == len(y_list)
    plt.figure(0, figsize=(10, 6), dpi=100)
    plt.rcParams['savefig.dpi'] = 300
    if label:
        for (x, y, l) in zip(x_list, y_list, label):
            plt.plot(x, y, label=l, linewidth=2)
        plt.legend()
    else:
        for (x, y) in zip(x_list, y_list):
            plt.plot(x, y)
    if axh:
        plt.axhline(axh, linestyle='--', color='#000000')
    if axv:
        plt.axvline(axv, linestyle='--', color='#000000')
    plt.title(name)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)

    if save_dictionary:
        # if not os.path.exists(save_dictionary):
        #     os.makedirs(save_dictionary)
        plt.savefig(save_dictionary + name + '.png')
    if show:
        plt.show()
    plt.close(0)


if __name__ == '__main__':
    # fp = "/home/centos/model/aa/1blk_40/0/struct_info.pkl"
    # m1 = csc.Molecule()
    # m1.create_from_info(fp)
    # atom3d = csc.Atom3D([0,0,0], [50,50,50])
    # atom3d.add_molecules(m1, 40)
    #
    # fp1 = "/home/centos/model/aa/1blk_40/0/1blk_40.data"
    # fp2 = "/home/centos/model/aa/1blk_40/0/1blk_40.lammpstrj"
    # aa = csc.AA3D()
    # aa.get_mass(fp1)
    # aa.renew_coordinate(fp2, step=int(15000000))
    # atom3d.renew_coordinate(aa)
    #
    # LmpParser.create_data_file(atom3d, "temp.data")

    fp = "../data/20200329_lr=0.1_bak/iter_00/pre/1blk_40.data"
    a = LmpParser.load_data_file(fp)
    LmpParser.create_data_file(a, "test.data")
    pass
