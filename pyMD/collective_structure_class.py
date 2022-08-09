from typing import List, Dict, Set
from pyMD import cartesian_operation as cso
from pyMD import basic_structure_class as bsc
from pyMD import functions as f
import numpy as np
import pickle as pkl
import copy as cp
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree



def sort_types(type_tuple):
    max_iter = int(len(type_tuple) / 2)
    for iter_ in range(max_iter):
        if type_tuple[iter_] < type_tuple[-iter_ - 1]:
            return type_tuple, False
        elif type_tuple[iter_] > type_tuple[-iter_ - 1]:
            return type_tuple[::-1], True
        elif iter_ == max_iter - 1:
            return type_tuple[::-1], True


class AA3D:
    mass_dict = {1: 12.011150, 2: 1.007970, 3: 12.011150, 4: 15.999400, 5: 15.999400, 6: 12.011150, 7: 12.011150,
                 8: 14.006700, 9: 1.007970, 10: 12.011150, 11: 15.999400}
    aa_type = {8, 9, 10, 11}

    def __init__(self):
        self.box_l = np.zeros(3)
        self.box_h = np.zeros(3)
        self.lattice_parameter = np.zeros(3)
        self.id_mass_dict = dict()
        self.id_type_dict = dict()
        self.id_coordinate_dict = dict()
        self.id_velocity_dict = dict()
        self.id_v2_dict = dict()

    def get_mass(self, file_path):
        with open(file_path, "r") as f:
            lines = f.readlines()
        flag = 0
        for line in lines:
            arg = line.split()
            if arg:
                if not flag:
                    if arg[0] == "Atoms":
                        flag = 1
                else:
                    if arg[0] == "Bonds":
                        break
                    self.id_mass_dict[int(arg[0])] = self.mass_dict[int(arg[2])]
                    self.id_type_dict[int(arg[0])] = int(arg[2])
        pass

    def cal_coordinate(self, co_s, image):
        return self.lattice_parameter * co_s + self.box_l + image * self.lattice_parameter

    def renew_coordinate(self, file_path, step=None, begin=0, lines=None, one_file=True, renew_v=False, pad=0,
                         only_aa=False, only_u=False):
        if one_file:
            if lines is None:
                with open(file_path, "r") as f:
                    lines = f.readlines()

            if step is not None:
                while True:
                    line = lines[begin]
                    if line == "ITEM: TIMESTEP\n":
                        # print(lines[begin + 1])
                        if int(lines[begin + 1].split()[0]) == step:
                            break
                    begin += 1
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
        self.box_l = np.array(box_l)
        self.box_h = np.array(box_h)
        self.lattice_parameter = self.box_h - self.box_l

        end = begin
        for line_idx, line in enumerate(lines[begin + 9:]):
            arg = line.split()
            if arg:
                if arg[0] == "ITEM:":
                    end += line_idx + 8
                    break
                atom_id = int(arg[0])
                if only_aa and self.id_type_dict[atom_id] not in self.aa_type:
                    continue
                if only_u and self.id_type_dict[atom_id] != 10:
                    continue
                co_s = np.asarray([float(arg[1]), float(arg[2]), float(arg[3])])
                image = np.asarray([int(arg[4]), int(arg[5]), int(arg[6])], dtype=int)
                self.id_coordinate_dict[atom_id] = self.cal_coordinate(co_s, image)
                if renew_v:
                    self.id_velocity_dict[atom_id] = np.asarray([float(_) for _ in arg[7:10]])
                    self.id_v2_dict[atom_id] = np.sum(np.square(self.id_velocity_dict[atom_id]))

        if one_file:
            return end, lines

    def cal_centroid(self, id_set):
        coordinate_list = list()
        mass_list = list()
        for idx in id_set:
            coordinate_list.append(self.id_coordinate_dict[idx])
            mass_list.append(self.id_mass_dict[idx])
        centroid = np.dot(np.array(mass_list), np.array(coordinate_list)) / np.sum(mass_list)
        return centroid

    def cal_mv2(self, id_set):
        mv2_aa = 0.
        mv2_cg = 0.
        sum_m = 0.
        for idx in id_set:
            mv2_aa += self.id_mass_dict[idx] * self.id_v2_dict[idx]
            mv2_cg += self.id_mass_dict[idx] * self.id_velocity_dict[idx]
            sum_m += self.id_mass_dict[idx]
        mv2_cg = np.sum(np.square(mv2_cg)) / sum_m
        return mv2_aa, mv2_cg


class Collection:
    # type_mass_dict = {"TO(0)": 73.1142, "TO(1)": 72.1062, "TO(2)": 72.1062, "TO(3)": 72.1062, "Es": 72.0593, "c1": 14.0267, "c2": 13.0187,
    #                   "c3": 12.0107, "o1": 15.9994, "n1": 15.0147, "Es(0)": 72.0593, "Ph(0)": 76.093640,
    #                   "U(0)": 58.039890, "B(0)": 83.112325}

    # old type define
    # type_mass_dict = {"TO(0)": 73.1142, "TO(1)": 72.1062, "TO(2)": 72.1062, "TO(3)": 72.1062, "Es": 72.0593, "c1": 14.0267, "c2": 13.0187,
    #                   "c3": 12.0107, "o1": 15.9994, "n1": 14.0067, "Es(0)": 72.0593, "Ph(0)": 76.093640,
    #                   "U(0)": 58.039890, "B(0)": 83.112325, "h1": 1.008}

    # # fg type define
    # print("use new fg")
    # type_mass_dict = {"TO(1)": 73.1142, "TO(2)": 72.1062, "Es": 72.0593, "c3": 14.0267, "c2": 13.0187,
    #                   "c1": 12.0107, "o1": 15.9994, "n1": 14.0067, "h1": 1.008, "Ph": 76.093640, "U": 58.039890,"Me": 14.0267}
    # improper_types_dict = {("n1", "c1", "n1", "o1"): [59.374, 0.], ("c2", "c2", "c2", "c3"): [7.8153, 0.],
    #                         ("c2", "c2", "c2", "n1"): [17.0526, 0.],
    #                        ("c1", "n1", "c2", "h1"): [4.4181, 0.0000],
    #                         # ("Es", "c2", "c2", "c2"): [17.0526, 0.],
    #                         }

    # # fgh type define
    # type_mass_dict = {"TO(1)": 73.1142, "TO(2)": 72.1062, "Es": 72.0593, "c1": 14.0267, "c2": 13.0187,
    #                   "c3": 12.0107, "o1": 15.9994, "n1": 14.0067, "h1": 1.008, "Ph": 76.093640, "U": 58.039890,"Me": 14.0267}
    # improper_types_dict = {("n1", "c3", "n1", "o1"): [59.374, 0.], ("c2", "c2", "c2", "c1"): [7.8153, 0.],
    #                         ("c2", "c2", "c2", "n1"): [17.0526, 0.],
    #                        ("c3", "n1", "c2", "h1"): [4.4181, 0.0000]}
    #                     #    ("Es", "c2", "c2", "c2"): [17.0526, 0.]}

    # cgu type define
    print("using CGU")
    type_mass_dict = {"TO(1)": 73.1142, "TO(2)": 72.1062, "TO(3)": 73.1142, "Es": 72.0593, "Me": 14.0267,
                      "Ph": 76.093640, "U": 58.039890,
                      "c1": 12.0107, "o1": 15.9994, "n1": 14.0067, "h1": 1.008, "sv1": 16, "sv2": 56, "sv": 75}
    improper_types_dict = {("n1", "c1", "n1", "o1"): [59.374, 0.], ("c1", "n1", "Ph", "h1"): [4.4181, 0.0000]}

    hard_type = {"U", "Ph", "Me", "c1", "n1", "h1", "o1"}

    def set_type_define(self, _type):
        assert _type in ["fg", "cgu", "aa"], "type not defined"
        if _type == "cgu":
            print("using CGU")
            self.__class__.type_mass_dict = {"TO(1)": 73.1142, "TO(2)": 72.1062, "TO(3)": 73.1142, "Es": 72.0593, "Me": 14.0267,
                      "Ph": 76.093640, "U": 58.039890,
                      "c1": 12.0107, "o1": 15.9994, "n1": 14.0067, "h1": 1.008, "sv1": 16, "sv2": 56, "sv": 75}
            self.__class__.improper_types_dict = {("n1", "c1", "n1", "o1"): [59.374, 0.], ("c1", "n1", "Ph", "h1"): [4.4181, 0.0000]}
        elif _type == "fg":
            print("use new fg")
            self.__class__.type_mass_dict = {"TO(1)": 73.1142, "TO(2)": 72.1062, "Es": 72.0593, "c3": 14.0267, "c2": 13.0187,
                      "c1": 12.0107, "o1": 15.9994, "n1": 14.0067, "h1": 1.008, "Ph": 76.093640, "U": 58.039890,"Me": 14.0267}
            self.__class__.improper_types_dict = {("n1", "c1", "n1", "o1"): [59.374, 0.], ("c2", "c2", "c2", "c3"): [7.8153, 0.],
                            ("c2", "c2", "c2", "n1"): [17.0526, 0.],
                           ("c1", "n1", "c2", "h1"): [4.4181, 0.0000],
                            # ("Es", "c2", "c2", "c2"): [17.0526, 0.],
                            }
        elif _type == "aa":
            print("Using aa type define")
            self.__class__.type_mass_dict = {"TO(1)": 73.1142, "TO(2)": 72.1062, "Es": 72.0593, "c3": 14.0267, "c2": 13.0187,
                      "c1": 12.0107, "o1": 15.9994, "n1": 14.0067, "h1": 1.008, "Ph": 76.093640, "U": 58.039890,"Me": 14.0267, 
                      "c": 12.0107, "o": 15.9994, "n": 14.0067, "h": 1.008, "s":1.}
            # self.__class__.improper_types_dict = {("n1", "c1", "n1", "o1"): [59.374, 0.], ("c2", "c2", "c2", "c3"): [7.8153, 0.],
            #                 ("c2", "c2", "c2", "n1"): [17.0526, 0.],
            #                ("c1", "n1", "c2", "h1"): [4.4181, 0.0000],
            #                 # ("Es", "c2", "c2", "c2"): [17.0526, 0.],
            #                 }

    def __init__(self):
        self.Atoms = dict()
        self.Bonds = dict()
        self.Angles = dict()
        self.Dihedrals = dict()
        self.Impropers = dict()
        self.Atom_type_dict = dict()
        self.Bond_type_dict = dict()
        self.Angle_type_dict = dict()
        self.Dihedral_type_dict = dict()
        self.Improper_type_dict = dict()
        self.mole_dict = dict()
        self.atom_connect_dict = dict()
        self.atom_connect_2_dict = dict()
        self.changed = False
        pass

    def __str__(self):
        string = ""
        for atom in self.Atoms.values():
            string += atom.__str__()
            string += "\n"
        return string

    def recreate_atom_type_dict(self):
        new_dict = dict()
        for atom_id, atom in self.Atoms.items():
            if atom.type not in new_dict:
                new_dict[atom.type] = [atom_id]
            else:
                new_dict[atom.type].append(atom_id)
        self.Atom_type_dict = new_dict

    def get_atom_type_mass(self, atom_type):
        if atom_type in self.type_mass_dict:
            return self.type_mass_dict[atom_type]
        elif atom_type[0] in self.type_mass_dict:
            return self.type_mass_dict[atom_type[0]]
        else:
            raise KeyError("atom type (%s) not defined!" % atom_type)

    def append_element(self, element):
        class_n = element.__class__.__name__
        self.__dict__[class_n + "s"][element.id] = element
        if element.type not in self.__dict__[class_n + "_type_dict"]:
            self.__dict__[class_n + "_type_dict"][element.type] = [element.id]
        else:
            self.__dict__[class_n + "_type_dict"][element.type].append(element.id)
        if not isinstance(element, bsc.Atom):
            for atom_id in element.atoms:
                self.Atoms[atom_id].__dict__[class_n + "s"].append(element.id)
        else:
            if element.mole_id in self.mole_dict:
                self.mole_dict[element.mole_id].append(element.id)
            else:
                self.mole_dict[element.mole_id] = [element.id]
        self.changed = True
        self.in_cell = False
        self.default = False

    # create a bond between two atoms, specified by their id
    def create_bond(self, atom_id1, atom_id2, mb=None):
        if mb is None:
            mb = self.get_max_bond_id()
        bond_type, reversed = sort_types((self.Atoms[atom_id1].type, self.Atoms[atom_id2].type))
        if reversed:
            bond = bsc.Bond(mb+1, atom_id2, atom_id1, _type=bond_type)
        else:
            bond = bsc.Bond(mb+1, atom_id1, atom_id2, _type=bond_type)
        self.append_element(bond)
        return mb+1

    # must recal angles, dihedrals and impropers after change topology
    def delete_atom(self, atom_id):
        connected = self.get_connected_atom_id(atom_id)
        atom = self.Atoms.pop(atom_id)
        self.Atom_type_dict[atom.type].remove(atom_id)
        if len(self.Atom_type_dict[atom.type]) == 0:
            self.Atom_type_dict.pop(atom.type)
        for bond_id in atom.Bonds:
            if bond_id in self.Bonds:
                bond = self.Bonds.pop(bond_id)
                self.Bond_type_dict[bond.type].remove(bond_id)
                if len(self.Bond_type_dict[bond.type]) == 0:
                    self.Bond_type_dict.pop(bond.type)
        # for joint_type in ["Bond"]:
        #     for joint_id in atom.__dict__[joint_type + "s"]:
        #         if joint_id in self.__dict__[joint_type + "s"]:
        #             joint = self.__dict__[joint_type + "s"].pop(joint_id)
        #             self.__dict__[joint_type + "_type_dict"][joint.type].remove(joint_id)
        #             if len(self.__dict__[joint_type + "_type_dict"][joint.type]) == 0:
        #                 self.__dict__[joint_type + "_type_dict"].pop(joint.type)
        for c_atom_id in connected:
            c_atom = self.Atoms[c_atom_id]
            for bond_id in atom.Bonds:
                if bond_id in c_atom.Bonds:
                    c_atom.Bonds.remove(bond_id)
        self.mole_dict[atom.mole_id].remove(atom_id)
        self.changed = True
        return atom

    def delete_bond(self, bond_id):
        bond = self.Bonds.pop(bond_id)
        self.Bond_type_dict[bond.type].remove(bond_id)
        if len(self.Bond_type_dict[bond.type]) == 0:
            self.Bond_type_dict.pop(bond.type)
        for atom_id in bond.atoms:
            self.Atoms[atom_id].Bonds.remove(bond_id)
        return bond
    
    def delete_bond_type(self, bond_type):
        bond_ids = list(self.Bond_type_dict[bond_type])
        for bond_id in bond_ids:
            self.delete_bond(bond_id)
        return len(bond_ids)

    def get_connected_atom_id(self, atom_id):
        if self.changed:
            self.atom_connect_dict = dict()
        if atom_id not in self.atom_connect_dict:
            connected = set()
            for bond_id in self.Atoms[atom_id].Bonds:
                connected.add(self.Bonds[bond_id].get_other_atom_id(atom_id))
            self.atom_connect_dict[atom_id] = connected
        self.changed = False
        return cp.copy(self.atom_connect_dict[atom_id])

    # direct connected atoms and not direct
    def get_connected_2(self, atom_id):
        if self.changed:
            self.atom_connect_2_dict = dict()
        if atom_id not in self.atom_connect_2_dict:
            connected = self.get_connected_atom_id(atom_id)
            new_c = set()
            for idx in connected:
                new_c = new_c | self.get_connected_atom_id(idx)
            self.atom_connect_2_dict[atom_id] = new_c | connected
        self.changed = False
        return cp.copy(self.atom_connect_2_dict[atom_id])

    # get all paths
    def traverse_bfs(self, begin_id, depth=3):
        def bfs_iterative(path, current_depth):
            if current_depth == depth:
                all_path.add(path)
            else:
                connected = self.get_connected_atom_id(path[-1])
                for atom_id in connected:
                    if atom_id not in path:
                        bfs_iterative(tuple(list(path) + [atom_id]), current_depth + 1)

        all_path = set()
        bfs_iterative((begin_id,), 1)
        return all_path

    # visit all atoms
    def tranverse_bfs_visit(self, begin_id, visit=None):
        if visit is not None:
            visit(begin_id, 0)
        visited = [begin_id]
        unvisited = list()
        unvisited_0 = set()
        for c_id in self.get_connected_atom_id(begin_id):
            unvisited.append((c_id, begin_id))
            unvisited_0.add(c_id)
        while unvisited:
            # print(visited)
            # print(unvisited)
            b_id, f_id = unvisited.pop(0)
            if visit is not None:
                visit(b_id, f_id)
            visited.append(b_id)
            for c_id in self.get_connected_atom_id(b_id):
                if c_id not in visited and c_id not in unvisited_0:
                    unvisited.append((c_id, b_id))
                    unvisited_0.add(c_id)
        return visited

    def set_block_id_mole(self, begin_block):
        def visit(a_id, f_id):
            this_a = self.Atoms[a_id]
            if f_id == 0:
                this_a.block = begin_block
            else:
                father_a = self.Atoms[f_id]
                if (father_a.type == "Es" and this_a.type in ["Ph", "c2"]) or (father_a.type in ["Ph", "c2"] and this_a.type == "Es"):
                    this_a.block = father_a.block + 1
                else:
                    this_a.block = father_a.block

        return visit

    def set_block_id(self, begin_block=1):
        visited = set()
        for a_id, atom in self.Atoms.items():
            if atom.type == "TO(1)":
                if a_id not in visited:
                    print(a_id)
                    visit = self.set_block_id_mole(begin_block=begin_block)
                    new_visit = self.tranverse_bfs_visit(a_id, visit)
                    for v in new_visit:
                        if self.Atoms[v].block > begin_block:
                            begin_block = self.Atoms[v].block
                        visited.add(v)
                    begin_block += 1
        return visited

    def renew_coordinate(self, aa3d: AA3D, only_aa=False, only_u=False):
        for atom in self.Atoms.values():
            if only_aa:
                if atom.type[0].islower():
                    atom.coordinate = aa3d.cal_centroid(atom.aa_ids)
            elif only_u:
                if atom.type[0] == "c":
                    atom.coordinate = aa3d.cal_centroid(atom.aa_ids)
            else:
                atom.coordinate = aa3d.cal_centroid(atom.aa_ids)
        pass

    def renew_mv2(self, aa3d: AA3D):
        for atom in self.Atoms.values():
            atom.mv2_aa, atom.mv2_cg = aa3d.cal_mv2(atom.aa_ids)
        pass

    def cal_all_angles(self):
        self.Angles = dict()
        # mole_angle = dict()
        # for mole, atom_ids in self.mole_dict.values():
        #     mole_angle[mole] = set()
        #     for atom_id in atom_ids:
        #         angles = self.traverse_bfs(atom_id, 3)
        #         for angle in angles:
        #             if angle not in mole_angle[mole] and angle[::-1] not in mole_angle[mole]:
        #                 mole_angle[mole].add(angle)

        all_angle_id = set()
        for atom_id in self.Atoms:
            angles = self.traverse_bfs(atom_id, 3)
            for angle in angles:
                if angle not in all_angle_id and angle[::-1] not in all_angle_id:
                    all_angle_id.add(angle)
        for idx, angle in enumerate(all_angle_id):
            _type, rev = sort_types(
                (self.Atoms[angle[0]].type, self.Atoms[angle[1]].type, self.Atoms[angle[2]].type))
            if rev:
                angle_ = bsc.Angle(idx + 1, *angle[::-1], _type=_type)
            else:
                angle_ = bsc.Angle(idx + 1, *angle, _type=_type)
            self.append_element(angle_)
            # self.Atoms[angle[0]].Angles.append(idx + 1)
            # self.Atoms[angle[1]].Angles.append(idx + 1)
            # self.Atoms[angle[2]].Angles.append(idx + 1)

    def cal_all_dihedrals(self):
        self.Dihedrals = dict()
        all_dihedral_id = set()
        for atom_id in self.Atoms:
            dihedrals = self.traverse_bfs(atom_id, 4)
            for dihedral in dihedrals:
                if dihedral not in all_dihedral_id and dihedral[::-1] not in all_dihedral_id:
                    all_dihedral_id.add(dihedral)
        for idx, dihedral in enumerate(all_dihedral_id):
            _type, rev = sort_types(
                (self.Atoms[dihedral[0]].type, self.Atoms[dihedral[1]].type, self.Atoms[dihedral[2]].type,
                 self.Atoms[dihedral[3]].type))
            if rev:
                dihedral_ = bsc.Dihedral(idx + 1, *dihedral[::-1], _type=_type)
            else:
                dihedral_ = bsc.Dihedral(idx + 1, *dihedral, _type=_type)
            self.append_element(dihedral_)
            # self.Atoms[dihedral[0]].Dihedrals.append(idx + 1)
            # self.Atoms[dihedral[1]].Dihedrals.append(idx + 1)
            # self.Atoms[dihedral[2]].Dihedrals.append(idx + 1)
            # self.Atoms[dihedral[3]].Dihedrals.append(idx + 1)

    def cal_all_impropers(self, improper="class2", improper_dict=None):
        def get_tuples(type_list):
            tuple_list = list()
            tuple1 = tuple(type_list)
            tuple_list.append(tuple1)
            tuple2 = (type_list[0], type_list[2], type_list[1])
            if tuple2 not in tuple_list:
                tuple_list.append(tuple2)
            tuple3 = (type_list[1], type_list[0], type_list[2])
            if tuple3 not in tuple_list:
                tuple_list.append(tuple3)
            tuple4 = (type_list[1], type_list[2], type_list[0])
            if tuple4 not in tuple_list:
                tuple_list.append(tuple4)
            tuple5 = (type_list[2], type_list[1], type_list[0])
            if tuple5 not in tuple_list:
                tuple_list.append(tuple5)
            tuple6 = (type_list[2], type_list[0], type_list[1])
            if tuple6 not in tuple_list:
                tuple_list.append(tuple6)
            return tuple_list

        def find_type(_aid, _it, _iid):
            improper_type = None
            for tail in get_tuples(_it):
                if improper == "DFF":
                    improper_type = (self.Atoms[_aid].type, self.Atoms[tail[0]].type, self.Atoms[tail[1]].type,
                                        self.Atoms[tail[2]].type,)
                else:
                    improper_type = (self.Atoms[tail[0]].type, self.Atoms[_aid].type, self.Atoms[tail[1]].type,
                                        self.Atoms[tail[2]].type,)
                    # print(improper_type)
                    
                if improper_type in improper_dict:
                    if improper == "DFF":
                        n_improper = bsc.Improper(_iid, _aid, tail[0], tail[1], tail[2],
                                                    _type=improper_type)
                    else:
                        n_improper = bsc.Improper(_iid, tail[0], _aid, tail[1], tail[2],
                                                    _type=improper_type)
                    self.append_element(n_improper)
                    _iid += 1
                    break
                improper_type = None
            if improper_type is None:
                if improper in ["ignore", "DFF"]:
                    print((self.Atoms[_aid].type, self.Atoms[_it[0]].type, self.Atoms[_it[1]].type,
                            self.Atoms[_it[2]].type,))
                    print("Improper type not found!")
                else:
                    print((self.Atoms[_aid].type, self.Atoms[_it[0]].type, self.Atoms[_it[1]].type,
                            self.Atoms[_it[2]].type,))
                    raise KeyError("Improper type not defined!")
            return _iid

        self.Impropers = dict()
        print(improper)
        if improper == "COMPASS":
            from pyMD.class2 import COMPASS
            improper_dict = COMPASS.type_map['Impropers']
            improper = "ignore"
        if improper_dict is None:
            improper_dict = self.improper_types_dict
        improper_id = 1
        for atom_id in self.Atoms:
            connected = self.get_connected_atom_id(atom_id)
            if len(connected) == 3:
                improper_tail = list(connected)
                improper_id = find_type(atom_id, improper_tail, improper_id)
            elif len(connected) == 4:
                for i in range(4):
                    improper_tail = list(connected)
                    improper_tail.pop(i)
                    improper_id = find_type(atom_id, improper_tail, improper_id)

    def clear_all_joints(self):
        self.Angles = dict()
        self.Angle_type_dict = dict()
        self.Dihedrals = dict()
        self.Dihedral_type_dict = dict()
        self.Impropers = dict()
        self.Improper_type_dict = dict()
        for atom in self.Atoms.values():
            atom.Angles = list()
            atom.Dihedrals = list()
            atom.Impropers = list()

    def cal_all_joints(self, **kwargs):
        self.cal_all_angles()
        self.cal_all_dihedrals()
        self.cal_all_impropers(**kwargs)

    def get_max_atom_id(self):
        return max(self.Atoms.keys())

    def get_max_bond_id(self):
        if not self.Bonds:
            return 0
        else:
            return max(self.Bonds.keys())

    def get_max_angle_id(self):
        return max(self.Angles.keys())

    def get_max_mole_id(self):
        m = 0
        for a in self.Atoms.values():
            m = a.mole_id if a.mole_id > m else m
        return m

    def sort_atom_bond_id(self):
        self.changed = True
        atom_reflect = dict()
        bond_reflect = dict()
        Atoms = self.Atoms
        self.Atoms = dict()
        self.Atom_type_dict = dict()
        Bonds = self.Bonds
        self.Bonds = dict()
        self.Bond_type_dict = dict()
        atom_id_list = list(Atoms.keys())
        atom_id_list.sort()
        for idx, atom_id in enumerate(atom_id_list):
            atom = Atoms.pop(atom_id)
            atom.id = idx + 1
            atom_reflect[atom_id] = idx + 1
            atom.Bonds = list()
            self.append_element(atom)
        for idx, bond_id in enumerate(list(Bonds.keys())):
            bond = Bonds.pop(bond_id)
            bond.id = idx + 1
            bond_atoms = list()
            for atom_id in bond.atoms:
                bond_atoms.append(atom_reflect[atom_id])
            bond.atoms = bond_atoms
            self.append_element(bond)
        return atom_reflect

    def get_mass(self):
        mass = 0.
        for atom in self.Atoms.values():
            mass += self.get_atom_type_mass(atom.type)
        return mass

    def pattern_mass(self, mole_pattern):
        mass_sum = 0.
        for pa in mole_pattern:
            mass_sum += self.get_atom_type_mass(pa)
        return mass_sum


class Molecule(Collection):

    def create_from_info(self, file_path):
        with open(file_path, "rb") as f:
            info_dict = pkl.load(f)
        # create atoms
        for idx, aa_ids in info_dict["id_ids"].items():
            atom = bsc.Atom(idx, _type=info_dict["id_type"][idx], aa_ids=aa_ids)
            self.append_element(atom)
        # create bonds
        for idx, bond in enumerate(info_dict["bonds"]):
            a1 = self.Atoms[bond[0]]
            a2 = self.Atoms[bond[1]]
            _type, rev = sort_types((a1.type, a2.type))
            if rev:
                bond_ = bsc.Bond(idx + 1, *bond[::-1], _type=_type)
            else:
                bond_ = bsc.Bond(idx + 1, *bond, _type=_type)
            self.append_element(bond_)
        self.cal_all_joints()

    def __str__(self):
        return "Molecule\n" + super().__str__()

    def cal_centroid(self):
        coordinate_list = list()
        mass_list = list()
        for atom in self.Atoms.values():
            coordinate_list.append(atom.coordinate)
            mass_list.append(self.get_atom_type_mass(atom.type))
        centroid = np.dot(np.array(mass_list), np.array(coordinate_list)) / np.sum(mass_list)
        return centroid

    def set_centroid(self, co):
        centroid = self.cal_centroid()
        vec = co - centroid
        for atom in self.Atoms.values():
            atom.move_along_vector(vec)

    def random_coordinate(self, start_id, start_co, param):
        len_dict = dict()
        for type_tuple, energy in param["Bond"].items():
            len_dict[type_tuple] = energy[1]
        for atom in self.Atoms.values():
            atom.create_co = False
        self.Atoms[start_id].coordinate = start_co.copy()
        self.Atoms[start_id].create_co = True
        expanding = list()
        for connect in self.get_connected_atom_id(start_id):
            expanding.append([start_id, connect])
        while expanding:
            pair = expanding.pop(0)
            _type, rev = sort_types((self.Atoms[pair[0]].type, self.Atoms[pair[1]].type))
            n_co = self.Atoms[pair[0]].coordinate + len_dict[_type] * cso.random_unit()
            self.Atoms[pair[1]].coordinate = n_co.copy()
            self.Atoms[pair[1]].create_co = True
            for connect in self.get_connected_atom_id(pair[1]):
                if not self.Atoms[connect].create_co:
                    expanding.append([pair[1], connect])


class RepeatUnit(Molecule):
    def __init__(self):
        super().__init__()
        self.head = 0
        self.end = 0


    def rotate_to_vector(self, vector):
        vector = np.array(vector)
        old_vector = self.Atoms[self.end].coordinate - self.Atoms[self.head].coordinate
        n_vector, theta = cso.cal_n_vector_theta(old_vector, vector)
        if theta < 1E-10:
            return
        centroid = self.cal_centroid()
        for atom in self.Atoms.values():
            atom.coordinate = cso.rotate_n_theta(atom.coordinate, centroid, n_vector, theta)

    def rotate_orientation(self, old_orient, new_orient):
        n_vector, theta = cso.cal_n_vector_theta(old_orient, new_orient)
        # ee_vector = self.Atoms[self.end].coordinate - self.Atoms[self.head].coordinate
        centroid = self.cal_centroid()
        for atom in self.Atoms.values():
            atom.coordinate = cso.rotate_n_theta(atom.coordinate, centroid, n_vector, theta)
            # atom.coordinate = cso.rotate_n_theta(atom.coordinate, centroid, ee_vector, theta)

    def rotate_n_random(self):
        centroid = self.cal_centroid()
        n_vector = self.Atoms[self.end].coordinate - self.Atoms[self.head].coordinate
        n_vector /= cso.calculate_distance(n_vector)
        theta = np.random.rand() * 2 * np.pi
        for atom in self.Atoms.values():
            atom.coordinate = cso.rotate_n_theta(atom.coordinate, centroid, n_vector, theta)


class Atom3D(Collection):
    def __init__(self, box_l=None, box_h=None):
        super().__init__()
        self.box_l = np.array(box_l) if box_l is not None else np.zeros(3)
        self.box_h = np.array(box_h) if box_h is not None else np.zeros(3)
        self.lattice_parameter = self.box_h - self.box_l
        self.in_cell = False
        self.default = False
        self.ce_dict = dict()

    def __str__(self):
        return "Atom3d, lattice parameter: " + self.lattice_parameter.__str__() + "\n" + super().__str__()

    def cal_charge(self):
        # charge_dict = {"n1": -0.462, "h1": 0.351, "o1": -0.5, "c3": 0.532, "c2": 0.095}
        charge_dict = {"n1": -0.595, "h1": 0.388, "c3": 0.884, "o1": -0.66, "c2": 0.095}

        for atom_id in self.Atoms:
            atom = self.Atoms[atom_id]
            if atom.type in charge_dict:
                if atom.type != "c2":
                    atom.charge = charge_dict[atom.type]
                else:
                    connected = self.get_connected_atom_id(atom_id)
                    connected_type = {self.Atoms[_].type for _ in connected}
                    if "n1" in connected_type:
                        atom.charge = charge_dict[atom.type]

    def cal_charge_aa(self, _type):
        from pyMD.class2 import COMPASS
        if _type == True:
            bond_inc = COMPASS.bond_increment
        elif _type == "no_h":
            bond_inc = COMPASS.bond_increment_no_h
        else:
            raise KeyError("bond inc not defined")
        for atom in self.Atoms.values():
            atom.charge = 0
            for connected in self.get_connected_atom_id(atom.id):
                type_tuple = tuple([atom.type, self.Atoms[connected].type])
                if type_tuple in bond_inc:
                    atom.charge += bond_inc[type_tuple]
                elif type_tuple[::-1] in bond_inc:
                    atom.charge -= bond_inc[type_tuple[::-1]]
                else:
                    print(type_tuple)
                    raise KeyError("Bond Increment Parameter Not Found!")
                # atom.charge = round(atom.charge, 3)


    def cal_hard_end(self, hard_end_list: List):
        def cal_hard_end_mole():
            def visit(a_id, f_id):
                this_a = self.Atoms[a_id]
                if this_a.type == "Es":
                    hard_end_list.append(a_id)

            return visit

        visited = set()
        visit_func = cal_hard_end_mole()
        for a_id, atom in self.Atoms.items():
            if atom.type == "TO(1)":
                if a_id not in visited:
                    new_visit = self.tranverse_bfs_visit(a_id, visit_func)
                    for v in new_visit:
                        visited.add(v)

    def cal_seg_end(self):
        def cal_end_mole():
            nonlocal hard_end_list, soft_end_list
            def visit(a_id, f_id):
                this_a = self.Atoms[a_id]
                if this_a.type == "Es":
                    hard_end_list.append(a_id)
                    soft_end_list.append(a_id)
            return visit

        visited = set()
        hard_end_list, soft_end_list = list(), list()
        visit_func = cal_end_mole()
        mole_num = 0
        for a_id, atom in self.Atoms.items():
            if atom.type == "TO(1)":
                if a_id not in visited:
                    mole_num += 1
                    new_visit = self.tranverse_bfs_visit(a_id, visit_func)
                    for v in new_visit:
                        visited.add(v)
        hard_end_list = np.asarray(hard_end_list)
        hard_end_list = np.reshape(hard_end_list, (-1, 2))

        # print(mole_num)
        soft_end_list = np.asarray(soft_end_list)
        soft_end_list = np.reshape(soft_end_list, (mole_num, -1))
        # print(soft_end_list.shape)

        soft_end_list = soft_end_list[:, 1:-1]
        # print(soft_end_list.shape)1111
        soft_end_list = np.reshape(soft_end_list, (-1, 2))
        return hard_end_list, soft_end_list

    def set_block_id(self, begin_block=1):
        visited = set()
        for a_id, atom in self.Atoms.items():
            if atom.type == "TO(1)":
                if a_id not in visited:
                    # print(a_id)
                    visit = self.set_block_id_mole(begin_block=begin_block)
                    new_visit = self.tranverse_bfs_visit(a_id, visit)
                    for v in new_visit:
                        if self.Atoms[v].block > begin_block:
                            begin_block = self.Atoms[v].block
                        visited.add(v)
                    begin_block += 1
        return visited

    def set_seg_id_mole(self, seg_dict:Dict, head_tail_set:Set, begin_block=1):
        def visit(a_id, f_id):
            this_a = self.Atoms[a_id]
            if f_id == 0:
                this_a.seg_id = begin_block
                seg_dict[this_a.seg_id] = [a_id]
            else:
                father_a = self.Atoms[f_id]
                if (father_a.type == "Es" and this_a.type in ["Ph", "c2"]) or (father_a.type in ["Ph", "c2"] and this_a.type == "Es"):
                    this_a.seg_id = father_a.seg_id + 1
                    seg_dict[this_a.seg_id] = [a_id]
                else:
                    this_a.seg_id = father_a.seg_id
                    seg_dict[this_a.seg_id].append(a_id)
            if this_a.type == "TO(1)":
                head_tail_set.add(this_a.seg_id)
        return visit

    # set_seg_id, get seg dict
    def set_seg_id(self, begin_block=1):
        visited = set()
        seg_dict = dict()
        head_tail_set = set()
        for a_id, atom in self.Atoms.items():
            if atom.type == "TO(1)":
                if a_id not in visited:
                    # print(a_id)
                    visit = self.set_seg_id_mole(seg_dict, head_tail_set, begin_block=begin_block)
                    new_visit = self.tranverse_bfs_visit(a_id, visit)
                    for v in new_visit:
                        if self.Atoms[v].seg_id > begin_block:
                            begin_block = self.Atoms[v].seg_id
                        visited.add(v)
                    begin_block += 1
        return seg_dict, head_tail_set



    def set_block_id_c1(self, begin_block=1):
        visited = set()
        for a_id, atom in self.Atoms.items():
            if atom.type == "TO(1)":
                if a_id not in visited:
                    print(a_id)
                    visit = self.set_block_id_mole(begin_block=begin_block)
                    new_visit = self.tranverse_bfs_visit(a_id, visit)
                    for v in new_visit:
                        if self.Atoms[v].block > begin_block:
                            begin_block = self.Atoms[v].block
                        visited.add(v)
                    begin_block += 1
        return visited

    def add_molecule(self, molecule: Molecule, aa_num=244):
        atom_num = len(self.Atoms)
        bond_num = len(self.Bonds)
        angle_num = len(self.Angles)
        dihedral_num = len(self.Dihedrals)
        improper_num = len(self.Impropers)
        mole_num = len(self.mole_dict)
        self.mole_dict[mole_num + 1] = list()
        for atom_id, atom in molecule.Atoms.items():
            self.mole_dict[mole_num + 1].append(atom_id + atom_num)
            n_atom = bsc.Atom(atom_id + atom_num, atom.type, _coordinate=atom.coordinate, mole_id=mole_num + 1)
            if aa_num > 0:
                for aa_id in atom.aa_ids:
                    n_atom.aa_ids.append(aa_id + mole_num * aa_num)
            n_atom.Bonds = [_ + bond_num for _ in atom.Bonds]
            n_atom.Angles = [_ + bond_num for _ in atom.Angles]
            n_atom.Dihedrals = [_ + bond_num for _ in atom.Dihedrals]
            n_atom.Impropers = [_ + bond_num for _ in atom.Impropers]
            self.append_element(n_atom)
        for bond_id, bond in molecule.Bonds.items():
            n_bond = bsc.Bond(bond_id + bond_num, *[_ + atom_num for _ in bond.atoms], _type=bond.type)
            self.append_element(n_bond)
        for angle_id, angle in molecule.Angles.items():
            n_angle = bsc.Angle(angle_id + angle_num, *[_ + atom_num for _ in angle.atoms], _type=angle.type)
            self.append_element(n_angle)
        for dihedral_id, dihedral in molecule.Dihedrals.items():
            n_dihedral = bsc.Dihedral(dihedral_id + dihedral_num, *[_ + atom_num for _ in dihedral.atoms],
                                      _type=dihedral.type)
            self.append_element(n_dihedral)
        for improper_id, improper in molecule.Impropers.items():
            n_improper = bsc.Improper(improper_id + improper_num, *[_ + atom_num for _ in improper.atoms],
                                      _type=improper.type)
            self.append_element(n_improper)

    def add_molecules(self, molecule: Molecule, num=1):
        for _ in range(num):
            self.add_molecule(molecule)

    def create_from_info(self, file_path, improper="class2"):
        with open(file_path, "rb") as f:
            info_dict = pkl.load(f)
        # create atoms
        # print(info_dict.keys())
        for idx, aa_ids in info_dict["id_ids"].items():
            atom = bsc.Atom(idx, _type=info_dict["id_type"][idx], aa_ids=aa_ids, mole_id=info_dict["mole_id"][idx])
            self.append_element(atom)
        # create bonds
        for idx, bond in enumerate(info_dict["bonds"]):
            a1 = self.Atoms[bond[0]]
            a2 = self.Atoms[bond[1]]
            _type, rev = sort_types((a1.type, a2.type))
            if rev:
                bond_ = bsc.Bond(idx + 1, *bond[::-1], _type=_type)
            else:
                bond_ = bsc.Bond(idx + 1, *bond, _type=_type)
            self.append_element(bond_)
        self.cal_all_joints(improper=improper)

    def create_from_pattern_box(self, mole_pattern, param, density=1.1, seed=None, anglecheck=True, box=(10,10,10)):
        def angle_check(t_angle, co1, co2, co3):
            c_angle = cso.calculate_angle(co1, co2, co3)
            if abs(c_angle - t_angle) <= 20:
                return True
            else:
                return False

        if seed is not None:
            np.random.seed(seed)
        assert len(box) == 3
        lattice_parameter = np.asarray(box)
        mole_num = box[0] * box[1] * box[2] / self.pattern_mass(mole_pattern) * density * 0.602
        self.lattice_parameter = lattice_parameter * int(mole_num) / mole_num
        mole_num = int(mole_num)
        self.box_l = np.zeros(3)
        self.box_h = self.lattice_parameter

        len_dict = dict()
        angle_dict = dict()
        for type_tuple, energy in param["Bond"].items():
            len_dict[type_tuple] = energy[1]
        for type_tuple, energy in param["Angle"].items():
            angle_dict[type_tuple] = energy[1]
        atom_id = 1
        bond_id = 1
        has_bond = False
        has_angle = 0
        last_type = ""
        ll_type = ""
        ll_co = None
        for mole_id in range(1, mole_num + 1):
            l_co = lattice_parameter * np.random.rand(3)
            for a_type in mole_pattern:
                if has_bond:
                    _type, rev = sort_types((last_type, a_type))
                    n_co = l_co + len_dict[_type] * cso.random_unit()
                    if has_angle >= 2 and anglecheck:
                        angle = angle_dict[sort_types((ll_type, last_type, a_type))[0]]
                        while not angle_check(angle, ll_co, l_co, n_co):
                            n_co = l_co + len_dict[_type] * cso.random_unit()
                    atom = bsc.Atom(atom_id, _type=a_type, mole_id=mole_id, _coordinate=n_co.copy())
                    self.append_element(atom)
                    ll_co = l_co
                    l_co = n_co
                    if rev:
                        bond = bsc.Bond(bond_id, atom_id, atom_id - 1, _type=_type)
                    else:
                        bond = bsc.Bond(bond_id, atom_id - 1, atom_id, _type=_type)
                    self.append_element(bond)
                    bond_id += 1
                else:
                    atom = bsc.Atom(atom_id, _type=a_type, mole_id=mole_id, _coordinate=l_co.copy())
                    self.append_element(atom)
                atom_id += 1
                has_bond = True
                has_angle += 1
                ll_type = last_type
                last_type = a_type
            has_bond = False
            has_angle = 0
        self.cal_all_joints()
        return mole_num

    def create_linear(self, mole_pattern:list, mole_num:int, param:dict, density=1.0, mass_dict=None, seed=None):
        """
        create linear polymers, placed in a cubic box randomly

        Parameters:
            mole_pattern: list of atom type in a single chain, e.g., ["A"] * 40
            mole_num: number of polymer chains
            param: dict of structure parameter, e.g, {"Bond": {("A", "A"): 1.0}}
            density: density of simulation box
            mass_dict: dict of atom mass, e.g., {"A": 50.}
            seed: random number seed

        Returns: None
        """

        # this function calculate mass of a single chain
        def cal_mass(_mole_pattern, _mass_dict):
            m = 0.
            for p in _mole_pattern:
                m += _mass_dict[p]
            return m
        
        # set random number seed
        if seed is not None:
            np.random.seed(seed)
        
        # calculate and set lattice parameter according to the mole pattern and mass dict. unit conversion is required
        lattice_parameter = pow(mole_num * cal_mass(mole_pattern, mass_dict) / density / 0.602, 1 / 3)
        self.lattice_parameter = np.asarray([lattice_parameter] * 3)
        self.box_l = np.zeros(3)
        self.box_h = self.lattice_parameter

        # initialize parameters
        atom_id = 1
        bond_id = 1

        # create all chains
        for mole_id in range(1, mole_num + 1):
            # reset parameters
            last_type = ""
            has_bond = False
            # create single chain
            # initialize coord of one atom randomly
            l_co = lattice_parameter * np.random.rand(3)
            for a_type in mole_pattern:
                # if not the first atom
                if has_bond:
                    _type, rev = sort_types((last_type, a_type))
                    # create new atom
                    n_co = l_co + param["Bond"][_type] * cso.random_unit()
                    atom = bsc.Atom(atom_id, _type=a_type, mole_id=mole_id, _coordinate=n_co.copy())
                    self.append_element(atom)
                    # record current coord
                    l_co = n_co
                    # connect the current atom to last atom
                    if rev:
                        bond = bsc.Bond(bond_id, atom_id, atom_id - 1, _type=_type)
                    else:
                        bond = bsc.Bond(bond_id, atom_id - 1, atom_id, _type=_type)
                    self.append_element(bond)
                    bond_id += 1
                # first atom
                else:
                    # create new atom
                    atom = bsc.Atom(atom_id, _type=a_type, mole_id=mole_id, _coordinate=l_co.copy())
                    self.append_element(atom)
                atom_id += 1
                has_bond = True
                has_angle += 1
                last_type = a_type
        
        # this function create angles and dihedrals automatically
        self.cal_all_joints()
        pass


    def create_from_pattern(self, mole_pattern, mole_num, param, density=1.1, seed=None, anglecheck=True):
        def angle_check(t_angle, co1, co2, co3):
            c_angle = cso.calculate_angle(co1, co2, co3)
            if abs(c_angle - t_angle) <= 20:
                return True
            else:
                return False

        if seed is not None:
            np.random.seed(seed)
        lattice_parameter = pow(mole_num * self.pattern_mass(mole_pattern) / density / 0.602, 1 / 3)
        self.lattice_parameter = np.asarray([lattice_parameter] * 3)
        self.box_l = np.zeros(3)
        self.box_h = self.lattice_parameter

        len_dict = dict()
        angle_dict = dict()
        for type_tuple, energy in param["Bond"].items():
            len_dict[type_tuple] = energy[1]
        for type_tuple, energy in param["Angle"].items():
            angle_dict[type_tuple] = energy[1]
        atom_id = 1
        bond_id = 1
        has_bond = False
        has_angle = 0
        last_type = ""
        ll_type = ""
        ll_co = None
        for mole_id in range(1, mole_num + 1):
            l_co = lattice_parameter * np.random.rand(3)
            for a_type in mole_pattern:
                if has_bond:
                    _type, rev = sort_types((last_type, a_type))
                    n_co = l_co + len_dict[_type] * cso.random_unit()
                    if has_angle >= 2 and anglecheck:
                        angle = angle_dict[sort_types((ll_type, last_type, a_type))[0]]
                        while not angle_check(angle, ll_co, l_co, n_co):
                            n_co = l_co + len_dict[_type] * cso.random_unit()
                    atom = bsc.Atom(atom_id, _type=a_type, mole_id=mole_id, _coordinate=n_co.copy())
                    self.append_element(atom)
                    ll_co = l_co
                    l_co = n_co
                    if rev:
                        bond = bsc.Bond(bond_id, atom_id, atom_id - 1, _type=_type)
                    else:
                        bond = bsc.Bond(bond_id, atom_id - 1, atom_id, _type=_type)
                    self.append_element(bond)
                    bond_id += 1
                else:
                    atom = bsc.Atom(atom_id, _type=a_type, mole_id=mole_id, _coordinate=l_co.copy())
                    self.append_element(atom)
                atom_id += 1
                has_bond = True
                has_angle += 1
                ll_type = last_type
                last_type = a_type
            has_bond = False
            has_angle = 0
        self.cal_all_joints()
        pass

    def create_from_pattern_wall(self, mole_pattern, mole_num, param, density=1.1, seed=None, wall=(1, 1, 1)):

        def angle_check(t_angle, co1, co2, co3):
            c_angle = cso.calculate_angle(co1, co2, co3)
            p = 180 * 180 - (c_angle - t_angle) * (c_angle - t_angle)
            if np.exp(-p) <= np.random.rand():
                return True
            else:
                return False

        def wall_check(_co):

            check = (_co > (self.lattice_parameter - 5.)) | (_co < 5.)
            # print(_co, self.lattice_parameter)
            # print(check)
            # print(not (check & wall).any())
            return not (check & wall).any()

        if seed is not None:
            np.random.seed(seed)
        lattice_parameter = pow(mole_num * self.pattern_mass(mole_pattern) / density / 0.602, 1 / 3)
        assert len(wall) == 3
        wall = np.asarray(wall, dtype=bool)
        self.lattice_parameter = np.asarray([lattice_parameter] * 3)
        self.box_l = np.zeros(3)
        self.box_h = self.lattice_parameter

        len_dict = dict()
        angle_dict = dict()
        for type_tuple, energy in param["Bond"].items():
            len_dict[type_tuple] = energy[1]
        for type_tuple, energy in param["Angle"].items():
            angle_dict[type_tuple] = energy[1]
        atom_id = 1
        bond_id = 1
        has_bond = False
        has_angle = 0
        last_type = ""
        ll_type = ""
        ll_co = None
        for mole_id in range(1, mole_num + 1):
            l_co = lattice_parameter * np.random.rand(3)
            for a_type in mole_pattern:
                if has_bond:
                    _type, rev = sort_types((last_type, a_type))
                    n_co = l_co + len_dict[_type] * cso.random_unit()
                    # print(n_co)
                    # print("call 1")
                    while not wall_check(n_co):
                        n_co = l_co + len_dict[_type] * cso.random_unit()
                    # print("call 2")
                    if has_angle >= 2:
                        angle = angle_dict[sort_types((ll_type, last_type, a_type))[0]]
                        while not (angle_check(angle, ll_co, l_co, n_co) and wall_check(n_co)):
                            n_co = l_co + len_dict[_type] * cso.random_unit()
                    # print(n_co)
                    atom = bsc.Atom(atom_id, _type=a_type, mole_id=mole_id, _coordinate=n_co.copy())
                    self.append_element(atom)
                    ll_co = l_co
                    l_co = n_co
                    if rev:
                        bond = bsc.Bond(bond_id, atom_id, atom_id - 1, _type=_type)
                    else:
                        bond = bsc.Bond(bond_id, atom_id - 1, atom_id, _type=_type)
                    self.append_element(bond)
                    bond_id += 1
                else:
                    atom = bsc.Atom(atom_id, _type=a_type, mole_id=mole_id, _coordinate=l_co.copy())
                    self.append_element(atom)
                atom_id += 1
                has_bond = True
                has_angle += 1
                ll_type = last_type
                last_type = a_type
            has_bond = False
            has_angle = 0
        self.cal_all_joints()
        pass

    def create_line(self, mole_pattern, mole_num, param, density=1.0, lines=2):
        len_dict = dict()
        angle_dict = dict()
        for type_tuple, energy in param["Bond"].items():
            len_dict[type_tuple] = energy[1]
        for type_tuple, energy in param["Angle"].items():
            angle_dict[type_tuple] = energy[1]
        V = mole_num * self.pattern_mass(mole_pattern) / density / 0.602
        l1 = mole_num / lines * 4.6
        l2 = np.sqrt(V / l1)
        self.lattice_parameter = np.asarray([l1, l2, l2])
        self.box_l = np.zeros(3)
        self.box_h = self.lattice_parameter
        center = mole_pattern.index("Me")
        hard_type = {"Me", "Ph", "U", "Es"}
        atom_id = 1
        bond_id = 1
        mole_id = 1
        for i in range(lines):
            for j in range(int(mole_num / lines)):
                center_co = np.asarray([2.5 + j * 4.6, l2 / 2, l2 / lines / 2 + i * (l2 / lines)])
                atom = bsc.Atom(atom_id, _type="Me", mole_id=mole_id, _coordinate=center_co.copy())
                self.append_element(atom)
                c_id = atom_id
                atom_id += 1
                l_co = center_co
                l_id = c_id
                last_type = "Me"
                for k in range(center + 1, len(mole_pattern)):
                    a_type = mole_pattern[k]
                    _type, rev = sort_types((last_type, a_type))
                    if a_type in hard_type:
                        n_co = l_co + len_dict[_type] * np.array([0., -1., 0.])
                    else:
                        n_co = l_co + len_dict[_type] * cso.random_unit()
                        i_co = n_co % self.lattice_parameter
                        z = (i_co[2] - l2 / lines / 2) % (l2 / lines)
                        # print("0", z)
                        while z < 2. or z > l2 / lines - 2:
                            n_co = l_co + len_dict[_type] * cso.random_unit()
                            i_co = n_co % self.lattice_parameter
                            z = (i_co[2] - l2 / lines / 2) % (l2 / lines)
                            # print("1", z)

                    atom = bsc.Atom(atom_id, _type=a_type, mole_id=mole_id, _coordinate=n_co.copy())
                    self.append_element(atom)
                    if rev:
                        bond = bsc.Bond(bond_id, l_id, atom_id, _type=_type)
                    else:
                        bond = bsc.Bond(bond_id, atom_id, l_id, _type=_type)
                    self.append_element(bond)
                    l_id = atom_id
                    l_co = n_co
                    last_type = a_type
                    bond_id += 1
                    atom_id += 1
                l_co = center_co
                l_id = c_id
                last_type = "Me"
                for k in range(center - 1, -1, -1):
                    a_type = mole_pattern[k]
                    _type, rev = sort_types((last_type, a_type))
                    if a_type in hard_type:
                        n_co = l_co + len_dict[_type] * np.array([0., 1., 0.])
                    else:
                        n_co = l_co + len_dict[_type] * cso.random_unit()
                    atom = bsc.Atom(atom_id, _type=a_type, mole_id=mole_id, _coordinate=n_co.copy())
                    self.append_element(atom)
                    if rev:
                        bond = bsc.Bond(bond_id, l_id, atom_id, _type=_type)
                    else:
                        bond = bsc.Bond(bond_id, atom_id, l_id, _type=_type)
                    self.append_element(bond)
                    l_id = atom_id
                    l_co = n_co
                    last_type = a_type
                    bond_id += 1
                    atom_id += 1
                mole_id += 1
        self.cal_all_joints()

    def create_line_tilt(self, mole_pattern, mole_num, param, density=1.0):
        len_dict = dict()
        angle_dict = dict()
        for type_tuple, energy in param["Bond"].items():
            len_dict[type_tuple] = energy[1]
        for type_tuple, energy in param["Angle"].items():
            angle_dict[type_tuple] = energy[1]
        V = mole_num * self.pattern_mass(mole_pattern) / density / 0.602
        root2 = 1.414213562
        duu = 4.6
        l1 = mole_num * duu / root2
        l2 = V / l1 / l1
        print(l2)
        self.lattice_parameter = np.asarray([l1, l2, l1])
        self.box_l = np.zeros(3)
        self.box_h = self.lattice_parameter
        center = mole_pattern.index("Me")
        hard_type = {"Me", "Ph", "U", "Es"}
        atom_id = 1
        bond_id = 1
        mole_id = 1
        c1 = duu * mole_num / 4
        c2 = c1 * 3
        for i in range(2):
            for j in range(int(mole_num / 2)):
                if i == 0:
                    x = (duu / 2 + j * duu) / root2
                    y = l2 / 2
                    z = (duu / 2 + j * duu) / root2 + l1 / 2
                else:
                    x = (duu / 2 + j * duu) / root2 + l1 / 2
                    y = l2 / 2
                    z = (duu / 2 + j * duu) / root2
                center_co = np.asarray([x, y, z])
                atom = bsc.Atom(atom_id, _type="Me", mole_id=mole_id, _coordinate=center_co.copy())
                self.append_element(atom)
                c_id = atom_id
                atom_id += 1
                l_co = center_co
                l_id = c_id
                last_type = "Me"
                for k in range(center + 1, len(mole_pattern)):
                    a_type = mole_pattern[k]
                    _type, rev = sort_types((last_type, a_type))
                    if a_type in hard_type:
                        n_co = l_co + len_dict[_type] * np.array([0., -1., 0.])
                    else:
                        n_co = l_co + len_dict[_type] * cso.random_unit()
                        i_co = n_co % self.lattice_parameter
                        z = np.sqrt(np.square(l1 - i_co[0]) + i_co[2] * i_co[2])
                        # print("0", z)
                        while c1 - 3.0 < z < c1 + 3.0 or c2 - 3.0 < z < c2 + 3.0:
                            n_co = l_co + len_dict[_type] * cso.random_unit()
                            i_co = n_co % self.lattice_parameter
                            z = np.sqrt(np.square(l1 - i_co[0]) + i_co[2] * i_co[2])
                            # print("1", z)

                    atom = bsc.Atom(atom_id, _type=a_type, mole_id=mole_id, _coordinate=n_co.copy())
                    self.append_element(atom)
                    if rev:
                        bond = bsc.Bond(bond_id, l_id, atom_id, _type=_type)
                    else:
                        bond = bsc.Bond(bond_id, atom_id, l_id, _type=_type)
                    self.append_element(bond)
                    l_id = atom_id
                    l_co = n_co
                    last_type = a_type
                    bond_id += 1
                    atom_id += 1
                l_co = center_co
                l_id = c_id
                last_type = "Me"
                for k in range(center - 1, -1, -1):
                    a_type = mole_pattern[k]
                    _type, rev = sort_types((last_type, a_type))
                    if a_type in hard_type:
                        n_co = l_co + len_dict[_type] * np.array([0., 1., 0.])
                    else:
                        n_co = l_co + len_dict[_type] * cso.random_unit()
                    atom = bsc.Atom(atom_id, _type=a_type, mole_id=mole_id, _coordinate=n_co.copy())
                    self.append_element(atom)
                    if rev:
                        bond = bsc.Bond(bond_id, l_id, atom_id, _type=_type)
                    else:
                        bond = bsc.Bond(bond_id, atom_id, l_id, _type=_type)
                    self.append_element(bond)
                    l_id = atom_id
                    l_co = n_co
                    last_type = a_type
                    bond_id += 1
                    atom_id += 1
                mole_id += 1
        self.cal_all_joints()

    def create_line_tilt2(self, mole_pattern, mole_num, param, density=1.0):
        len_dict = dict()
        angle_dict = dict()
        for type_tuple, energy in param["Bond"].items():
            len_dict[type_tuple] = energy[1]
        for type_tuple, energy in param["Angle"].items():
            angle_dict[type_tuple] = energy[1]
        V = mole_num * self.pattern_mass(mole_pattern) / density / 0.602
        root2 = 1.414213562
        duu = 4.6
        l1 = mole_num * duu / root2
        l2 = V / l1 / l1
        print(l2)
        self.lattice_parameter = np.asarray([l1, l2, l1])
        self.box_l = np.zeros(3)
        self.box_h = self.lattice_parameter
        center = mole_pattern.index("Me")
        hard_type = {"Me", "Ph", "U", "Es"}
        atom_id = 1
        bond_id = 1
        mole_id = 1
        c1 = duu * mole_num / 4
        c2 = c1 * 3
        for i in range(2):
            for j in range(int(mole_num / 2)):
                if i == 0:
                    x = (duu / 2 + j * duu) / root2
                    y = l2 / 2
                    z = (duu / 2 + j * duu) / root2 + l1 / 2
                else:
                    x = (duu / 2 + j * duu) / root2 + l1 / 2
                    y = l2 / 2
                    z = (duu / 2 + j * duu) / root2
                center_co = np.asarray([x, y, z])
                atom = bsc.Atom(atom_id, _type="Me", mole_id=mole_id, _coordinate=center_co.copy())
                self.append_element(atom)
                c_id = atom_id
                atom_id += 1
                l_co = center_co
                l_id = c_id
                last_type = "Me"
                for k in range(center + 1, len(mole_pattern)):
                    a_type = mole_pattern[k]
                    _type, rev = sort_types((last_type, a_type))
                    if a_type in hard_type:
                        n_co = l_co + len_dict[_type] * np.array([root2 / 2, 0., -root2 / 2])
                    else:
                        n_co = l_co + len_dict[_type] * cso.random_unit()
                        i_co = n_co % self.lattice_parameter
                        y = i_co[1]
                        while -3. < y - l2 / 2 < 3.0:
                            n_co = l_co + len_dict[_type] * cso.random_unit()
                            i_co = n_co % self.lattice_parameter
                            y = i_co[1]

                    atom = bsc.Atom(atom_id, _type=a_type, mole_id=mole_id, _coordinate=n_co.copy())
                    self.append_element(atom)
                    if rev:
                        bond = bsc.Bond(bond_id, l_id, atom_id, _type=_type)
                    else:
                        bond = bsc.Bond(bond_id, atom_id, l_id, _type=_type)
                    self.append_element(bond)
                    l_id = atom_id
                    l_co = n_co
                    last_type = a_type
                    bond_id += 1
                    atom_id += 1
                l_co = center_co
                l_id = c_id
                last_type = "Me"
                for k in range(center - 1, -1, -1):
                    a_type = mole_pattern[k]
                    _type, rev = sort_types((last_type, a_type))
                    if a_type in hard_type:
                        n_co = l_co + len_dict[_type] * np.array([-root2 / 2, 0., root2 / 2])
                    else:
                        n_co = l_co + len_dict[_type] * cso.random_unit()
                        i_co = n_co % self.lattice_parameter
                        y = i_co[1]
                        while -3. < y - l2 / 2 < 3.0:
                            n_co = l_co + len_dict[_type] * cso.random_unit()
                            i_co = n_co % self.lattice_parameter
                            y = i_co[1]
                    atom = bsc.Atom(atom_id, _type=a_type, mole_id=mole_id, _coordinate=n_co.copy())
                    self.append_element(atom)
                    if rev:
                        bond = bsc.Bond(bond_id, l_id, atom_id, _type=_type)
                    else:
                        bond = bsc.Bond(bond_id, atom_id, l_id, _type=_type)
                    self.append_element(bond)
                    l_id = atom_id
                    l_co = n_co
                    last_type = a_type
                    bond_id += 1
                    atom_id += 1
                mole_id += 1
        self.cal_all_joints()

    def create_from_pattern_center_p(self, mole_pattern, mole_num, param, density=1.1, center=0, p=None, seed=None):

        def angle_check(t_angle, co1, co2, co3):
            c_angle = cso.calculate_angle(co1, co2, co3)
            if abs(c_angle - t_angle) <= 20:
                return True
            else:
                return False

        if seed is not None:
            np.random.seed(seed)
        lattice_parameter = pow(mole_num * self.pattern_mass(mole_pattern) / density / 0.602, 1 / 3)
        self.lattice_parameter = np.asarray([lattice_parameter] * 3)
        self.box_l = np.zeros(3)
        self.box_h = self.lattice_parameter
        if p is None:
            p = [f.pdf_uniform, f.pdf_uniform, f.pdf_uniform]

        len_dict = dict()
        angle_dict = dict()
        for type_tuple, energy in param["Bond"].items():
            len_dict[type_tuple] = energy[1]
        for type_tuple, energy in param["Angle"].items():
            angle_dict[type_tuple] = energy[1]
        bond_id = 1
        has_bond = False
        has_angle = 0
        last_type = ""
        ll_type = ""
        ll_co = None

        for mole_id in range(1, mole_num + 1):
            # l_co = lattice_parameter * np.random.rand(3)
            r_co_list = list()
            for i in range(3):
                r_co = f.random_gen(p[i], 0, 1)
                r_co_list.append(r_co)
            l_co = lattice_parameter * np.asarray(r_co_list)
            # before center
            for a_idx, a_type in enumerate(mole_pattern[center::-1]):
                atom_id = (mole_id - 1) * len(mole_pattern) + center - a_idx + 1
                if has_bond:
                    _type, rev = sort_types((last_type, a_type))
                    n_co = l_co + len_dict[_type] * cso.random_unit()
                    if has_angle >= 2:
                        angle = angle_dict[sort_types((ll_type, last_type, a_type))[0]]
                        while not angle_check(angle, ll_co, l_co, n_co):
                            n_co = l_co + len_dict[_type] * cso.random_unit()
                    atom = bsc.Atom(atom_id, _type=a_type, mole_id=mole_id, _coordinate=n_co.copy())
                    self.append_element(atom)
                    ll_co = l_co
                    l_co = n_co
                    if rev:
                        bond = bsc.Bond(bond_id, atom_id, atom_id + 1, _type=_type)
                    else:
                        bond = bsc.Bond(bond_id, atom_id + 1, atom_id, _type=_type)
                    self.append_element(bond)
                    bond_id += 1
                else:
                    atom = bsc.Atom(atom_id, _type=a_type, mole_id=mole_id, _coordinate=l_co.copy())
                    self.append_element(atom)
                has_bond = True
                has_angle += 1
                ll_type = last_type
                last_type = a_type

            # after center
            has_bond = True
            has_angle = 2
            l_a = self.Atoms[(mole_id - 1) * len(mole_pattern) + center + 1]
            ll_a = self.Atoms[(mole_id - 1) * len(mole_pattern) + center]
            last_type = l_a.type
            ll_type = ll_a.type
            l_co = l_a.coordinate
            ll_co = ll_a.coordinate
            for a_idx, a_type in enumerate(mole_pattern[center + 1:]):
                atom_id = (mole_id - 1) * len(mole_pattern) + center + a_idx + 2
                if has_bond:
                    _type, rev = sort_types((last_type, a_type))
                    n_co = l_co + len_dict[_type] * cso.random_unit()
                    if has_angle >= 2:
                        angle = angle_dict[sort_types((ll_type, last_type, a_type))[0]]
                        while not angle_check(angle, ll_co, l_co, n_co):
                            n_co = l_co + len_dict[_type] * cso.random_unit()
                    atom = bsc.Atom(atom_id, _type=a_type, mole_id=mole_id, _coordinate=n_co.copy())
                    self.append_element(atom)
                    ll_co = l_co
                    l_co = n_co
                    if rev:
                        bond = bsc.Bond(bond_id, atom_id, atom_id - 1, _type=_type)
                    else:
                        bond = bsc.Bond(bond_id, atom_id - 1, atom_id, _type=_type)
                    self.append_element(bond)
                    bond_id += 1
                else:
                    atom = bsc.Atom(atom_id, _type=a_type, mole_id=mole_id, _coordinate=l_co.copy())
                    self.append_element(atom)
                has_bond = True
                ll_type = last_type
                last_type = a_type
            has_bond = False
            has_angle = 0
        self.cal_all_joints()
        pass

    def create_from_pattern_fix(self, mole_pattern, mole_num, param, density=1.1, center=0):

        def angle_check(t_angle, co1, co2, co3):
            c_angle = cso.calculate_angle(co1, co2, co3)
            if abs(c_angle - t_angle) <= 20:
                return True
            else:
                return False

        lattice_parameter = pow(mole_num * self.pattern_mass(mole_pattern) / density / 0.602, 1 / 3)
        self.lattice_parameter = np.asarray([lattice_parameter] * 3)
        self.box_l = np.zeros(3)
        self.box_h = self.lattice_parameter

        len_dict = dict()
        angle_dict = dict()
        for type_tuple, energy in param["Bond"].items():
            len_dict[type_tuple] = energy[1]
        for type_tuple, energy in param["Angle"].items():
            angle_dict[type_tuple] = energy[1]
        bond_id = 1
        has_bond = False
        has_angle = 0
        last_type = ""
        ll_type = ""
        ll_co = None

        for mole_id in range(1, mole_num + 1):
            chain_idx = (mole_id - 1) // (mole_num // 4)
            in_idx = ((mole_id - 1) % (mole_num // 4))
            if in_idx % 3 == 0:
                delta = -0.05
            elif in_idx % 3 == 1:
                delta = 0.
            else:
                delta = 0.05
            y = chain_idx % 2
            z = chain_idx // 2
            l_co = lattice_parameter * np.asarray([in_idx // 3 * 3 * 4 / mole_num, y / 2 + 0.25, z / 2 + 0.25 + delta])
            # before center
            for a_idx, a_type in enumerate(mole_pattern[center::-1]):
                atom_id = (mole_id - 1) * len(mole_pattern) + center - a_idx + 1
                if has_bond:
                    _type, rev = sort_types((last_type, a_type))
                    n_co = l_co + len_dict[_type] * cso.random_unit()
                    if has_angle >= 2:
                        angle = angle_dict[sort_types((ll_type, last_type, a_type))[0]]
                        while not angle_check(angle, ll_co, l_co, n_co):
                            n_co = l_co + len_dict[_type] * cso.random_unit()
                    atom = bsc.Atom(atom_id, _type=a_type, mole_id=mole_id, _coordinate=n_co.copy())
                    self.append_element(atom)
                    ll_co = l_co
                    l_co = n_co
                    if rev:
                        bond = bsc.Bond(bond_id, atom_id, atom_id + 1, _type=_type)
                    else:
                        bond = bsc.Bond(bond_id, atom_id + 1, atom_id, _type=_type)
                    self.append_element(bond)
                    bond_id += 1
                else:
                    atom = bsc.Atom(atom_id, _type=a_type, mole_id=mole_id, _coordinate=l_co.copy())
                    self.append_element(atom)
                has_bond = True
                has_angle += 1
                ll_type = last_type
                last_type = a_type

            # after center
            has_bond = True
            has_angle = 2
            l_a = self.Atoms[(mole_id - 1) * len(mole_pattern) + center + 1]
            ll_a = self.Atoms[(mole_id - 1) * len(mole_pattern) + center]
            last_type = l_a.type
            ll_type = ll_a.type
            l_co = l_a.coordinate
            ll_co = ll_a.coordinate
            for a_idx, a_type in enumerate(mole_pattern[center + 1:]):
                atom_id = (mole_id - 1) * len(mole_pattern) + center + a_idx + 2
                if has_bond:
                    _type, rev = sort_types((last_type, a_type))
                    n_co = l_co + len_dict[_type] * cso.random_unit()
                    if has_angle >= 2:
                        angle = angle_dict[sort_types((ll_type, last_type, a_type))[0]]
                        while not angle_check(angle, ll_co, l_co, n_co):
                            n_co = l_co + len_dict[_type] * cso.random_unit()
                    atom = bsc.Atom(atom_id, _type=a_type, mole_id=mole_id, _coordinate=n_co.copy())
                    self.append_element(atom)
                    ll_co = l_co
                    l_co = n_co
                    if rev:
                        bond = bsc.Bond(bond_id, atom_id, atom_id - 1, _type=_type)
                    else:
                        bond = bsc.Bond(bond_id, atom_id - 1, atom_id, _type=_type)
                    self.append_element(bond)
                    bond_id += 1
                else:
                    atom = bsc.Atom(atom_id, _type=a_type, mole_id=mole_id, _coordinate=l_co.copy())
                    self.append_element(atom)
                has_bond = True
                ll_type = last_type
                last_type = a_type
            has_bond = False
            has_angle = 0
        self.cal_all_joints()
        pass

    def create_from_pattern_fix2(self, mole_pattern, mole_num, param, density=1.1, center=0, soft=True):

        def angle_check(t_angle, co1, co2, co3):
            c_angle = cso.calculate_angle(co1, co2, co3)
            if abs(c_angle - t_angle) <= 20:
                return True
            else:
                return False

        def random_new_co(_l_co, _ll_co, _bond_type, _has_angle, _ll_type, _last_type, _a_type):
            _n_co = _l_co + len_dict[_bond_type] * cso.random_unit()
            if _has_angle >= 2:
                _angle = angle_dict[sort_types((_ll_type, _last_type, _a_type))[0]]
                while not angle_check(_angle, _ll_co, _l_co, _n_co):
                    _n_co = _l_co + len_dict[_bond_type] * cso.random_unit()
            return _n_co

        lattice_parameter = pow(mole_num * self.pattern_mass(mole_pattern) / density / 0.602, 1 / 3)
        self.lattice_parameter = np.asarray([lattice_parameter] * 3)
        self.box_l = np.zeros(3)
        self.box_h = self.lattice_parameter

        if not soft:
            hard_pattern = []
            for a_type in mole_pattern:
                if a_type[:2] != "TO":
                    hard_pattern.append(a_type)
            mole_pattern = hard_pattern
            for a_type in mole_pattern:
                if a_type[:2] != "TO":
                    center -= 1
                else:
                    break

        len_dict = dict()
        angle_dict = dict()
        for type_tuple, energy in param["Bond"].items():
            len_dict[type_tuple] = energy[1]
        for type_tuple, energy in param["Angle"].items():
            angle_dict[type_tuple] = energy[1]
        bond_id = 1
        has_bond = False
        has_angle = 0
        last_type = ""
        ll_type = ""
        ll_co = None

        for mole_id in range(1, mole_num + 1):
            chain_idx = (mole_id - 1) // (mole_num // 4)
            in_idx = ((mole_id - 1) % (mole_num // 4))
            if in_idx % 3 == 0:
                delta = -0.05
            elif in_idx % 3 == 1:
                delta = 0.
            else:
                delta = 0.05
            y = chain_idx % 2
            z = chain_idx // 2
            l_co = lattice_parameter * np.asarray([in_idx // 3 * 3 * 4 / mole_num, y / 2 + 0.25, z / 2 + 0.25 + delta])
            # before center
            for a_idx, a_type in enumerate(mole_pattern[center::-1]):
                atom_id = (mole_id - 1) * len(mole_pattern) + center - a_idx + 1
                if has_bond:
                    bond_type, rev = sort_types((last_type, a_type))
                    if a_type in self.hard_type:
                        n_co = l_co - np.asarray([0., len_dict[bond_type], 0.])
                    else:
                        n_co = random_new_co(l_co, ll_co, bond_type, has_angle, ll_type, last_type, a_type)
                    atom = bsc.Atom(atom_id, _type=a_type, mole_id=mole_id, _coordinate=n_co.copy())
                    self.append_element(atom)
                    ll_co = l_co
                    l_co = n_co
                    if rev:
                        bond = bsc.Bond(bond_id, atom_id, atom_id + 1, _type=bond_type)
                    else:
                        bond = bsc.Bond(bond_id, atom_id + 1, atom_id, _type=bond_type)
                    self.append_element(bond)
                    bond_id += 1
                else:
                    atom = bsc.Atom(atom_id, _type=a_type, mole_id=mole_id, _coordinate=l_co.copy())
                    self.append_element(atom)
                has_bond = True
                has_angle += 1
                ll_type = last_type
                last_type = a_type

            # after center
            has_bond = True
            has_angle = 2
            l_a = self.Atoms[(mole_id - 1) * len(mole_pattern) + center + 1]
            ll_a = self.Atoms[(mole_id - 1) * len(mole_pattern) + center]
            last_type = l_a.type
            ll_type = ll_a.type
            l_co = l_a.coordinate
            ll_co = ll_a.coordinate
            for a_idx, a_type in enumerate(mole_pattern[center + 1:]):
                atom_id = (mole_id - 1) * len(mole_pattern) + center + a_idx + 2
                if has_bond:
                    bond_type, rev = sort_types((last_type, a_type))
                    if a_type in self.hard_type:
                        n_co = l_co + np.asarray([0., len_dict[bond_type], 0.])
                    else:
                        n_co = random_new_co(l_co, ll_co, bond_type, has_angle, ll_type, last_type, a_type)
                    atom = bsc.Atom(atom_id, _type=a_type, mole_id=mole_id, _coordinate=n_co.copy())
                    self.append_element(atom)
                    ll_co = l_co
                    l_co = n_co
                    if rev:
                        bond = bsc.Bond(bond_id, atom_id, atom_id - 1, _type=bond_type)
                    else:
                        bond = bsc.Bond(bond_id, atom_id - 1, atom_id, _type=bond_type)
                    self.append_element(bond)
                    bond_id += 1
                else:
                    atom = bsc.Atom(atom_id, _type=a_type, mole_id=mole_id, _coordinate=l_co.copy())
                    self.append_element(atom)
                has_bond = True
                ll_type = last_type
                last_type = a_type
            has_bond = False
            has_angle = 0
        self.cal_all_joints()
        pass

    def create_from_patterns(self, mole_patterns, mole_nums, param, density=1.1, seed=None):

        def angle_check(t_angle, co1, co2, co3):
            c_angle = cso.calculate_angle(co1, co2, co3)
            if abs(c_angle - t_angle) <= 20:
                return True
            else:
                return False

        def add_pattern(mole_pattern, begin_a_id, begin_b_id, begin_m_id):
            a_id = begin_a_id
            b_id = begin_b_id
            m_id = begin_m_id
            has_bond = False
            has_angle = 0
            last_type = ""
            ll_type = ""
            ll_co = None
            l_co = lattice_parameter * np.random.rand(3)
            for a_type in mole_pattern:
                if has_bond:
                    _type, rev = sort_types((last_type, a_type))
                    n_co = l_co + len_dict[_type] * cso.random_unit()
                    if has_angle >= 2:
                        angle = angle_dict[sort_types((ll_type, last_type, a_type))[0]]
                        while not angle_check(angle, ll_co, l_co, n_co):
                            n_co = l_co + len_dict[_type] * cso.random_unit()
                    atom = bsc.Atom(a_id, _type=a_type, mole_id=m_id, _coordinate=n_co.copy())
                    self.append_element(atom)
                    ll_co = l_co
                    l_co = n_co
                    if rev:
                        bond = bsc.Bond(b_id, a_id, a_id - 1, _type=_type)
                    else:
                        bond = bsc.Bond(b_id, a_id - 1, a_id, _type=_type)
                    self.append_element(bond)
                    b_id += 1
                else:
                    atom = bsc.Atom(a_id, _type=a_type, mole_id=m_id, _coordinate=l_co.copy())
                    self.append_element(atom)
                a_id += 1
                has_bond = True
                has_angle += 1
                ll_type = last_type
                last_type = a_type
            m_id += 1
            return a_id, b_id, m_id

        if seed is not None:
            np.random.seed(seed)
        mass_sum = 0.
        for pat_idx, mole_pat in enumerate(mole_patterns):
            mass_sum += mole_nums[pat_idx] * self.pattern_mass(mole_pat)
        lattice_parameter = pow(mass_sum / density / 0.602, 1 / 3)
        self.lattice_parameter = np.asarray([lattice_parameter] * 3)
        self.box_l = np.zeros(3)
        self.box_h = self.lattice_parameter
        len_dict = dict()
        angle_dict = dict()
        for type_tuple, energy in param["Bond"].items():
            len_dict[type_tuple] = energy[1]
        for type_tuple, energy in param["Angle"].items():
            angle_dict[type_tuple] = energy[1]
        atom_id, bond_id, mole_id = (1, 1, 1)
        for pat_idx, mole_pat in enumerate(mole_patterns):
            for _ in range(mole_nums[pat_idx]):
                atom_id, bond_id, mole_id = add_pattern(mole_pat, atom_id, bond_id, mole_id)
        self.cal_all_joints()
        pass

    def create_from_patterns_exclude(self, mole_patterns, mole_nums, param, density=1.1, seed=None):

        def angle_check(t_angle, co1, co2, co3):
            c_angle = cso.calculate_angle(co1, co2, co3)
            if abs(c_angle - t_angle) <= 40:
                return True
            else:
                return False

        def exclude_check(co, mole_co_list, _tree, rc=2.):
            r2 = rc * rc
            for m_co in mole_co_list[:-1]:
                if np.square(m_co - co).sum() < r2:
                    print("1 error")
                    return False
            if tree is not None:
                if len(_tree.query_ball_point(co, rc)) > 0:
                    print("2 error")
                    return False
            return True

        def add_pattern(mole_pattern, begin_a_id, begin_b_id, begin_m_id):
            a_id = begin_a_id
            b_id = begin_b_id
            m_id = begin_m_id
            has_bond = False
            has_angle = 0
            last_type = ""
            ll_type = ""
            ll_co = None
            mole_co_list = list()
            l_co = lattice_parameter * np.random.rand(3)
            while not exclude_check(l_co, mole_co_list, tree, rc=3.):
                l_co = lattice_parameter * np.random.rand(3)
            mole_co_list.append(l_co.copy())
            for a_type in mole_pattern:
                if has_bond:
                    _type, rev = sort_types((last_type, a_type))
                    n_co = l_co + len_dict[_type] * cso.random_unit() * np.random.uniform(0.9, 1.1)
                    n_co_in = n_co % lattice_parameter
                    while not exclude_check(n_co_in, mole_co_list, tree):
                        n_co = l_co + len_dict[_type] * cso.random_unit() * np.random.uniform(0.9, 1.1)
                        n_co_in = n_co % lattice_parameter
                    if has_angle >= 2:
                        angle = angle_dict[sort_types((ll_type, last_type, a_type))[0]]
                        while not angle_check(angle, ll_co, l_co, n_co):
                            n_co = l_co + len_dict[_type] * cso.random_unit() * np.random.uniform(0.9, 1.1)
                            n_co_in = n_co % lattice_parameter
                            while not exclude_check(n_co_in, mole_co_list, tree):
                                n_co = l_co + len_dict[_type] * cso.random_unit()
                                n_co_in = n_co % lattice_parameter
                    mole_co_list.append(n_co_in.copy())
                    atom = bsc.Atom(a_id, _type=a_type, mole_id=m_id, _coordinate=n_co.copy())
                    self.append_element(atom)
                    ll_co = l_co
                    l_co = n_co
                    if rev:
                        bond = bsc.Bond(b_id, a_id, a_id - 1, _type=_type)
                    else:
                        bond = bsc.Bond(b_id, a_id - 1, a_id, _type=_type)
                    self.append_element(bond)
                    b_id += 1
                else:
                    atom = bsc.Atom(a_id, _type=a_type, mole_id=m_id, _coordinate=l_co.copy())
                    self.append_element(atom)
                a_id += 1
                has_bond = True
                has_angle += 1
                ll_type = last_type
                last_type = a_type
                print(a_id)
            m_id += 1
            return a_id, b_id, m_id, mole_co_list

        if seed is not None:
            np.random.seed(seed)
        mass_sum = 0.
        for pat_idx, mole_pat in enumerate(mole_patterns):
            mass_sum += mole_nums[pat_idx] * self.pattern_mass(mole_pat)
        lattice_parameter = pow(mass_sum / density / 0.602, 1 / 3)
        self.lattice_parameter = np.asarray([lattice_parameter] * 3)
        self.box_l = np.zeros(3)
        self.box_h = self.lattice_parameter
        len_dict = dict()
        angle_dict = dict()
        for type_tuple, energy in param["Bond"].items():
            len_dict[type_tuple] = energy[1]
        for type_tuple, energy in param["Angle"].items():
            angle_dict[type_tuple] = energy[1]
        atom_id, bond_id, mole_id = (1, 1, 1)
        atom_co_list = list()
        tree = None
        for pat_idx, mole_pat in enumerate(mole_patterns):
            for _ in range(mole_nums[pat_idx]):
                atom_id, bond_id, mole_id, m_co_list = add_pattern(mole_pat, atom_id, bond_id, mole_id)
                atom_co_list += m_co_list
                tree = cKDTree(atom_co_list)
        self.cal_all_joints()
        pass

    def create_from_molecule(self, mole: Molecule, mole_num, param, density=1.1, seed=None):
        if seed is not None:
            np.random.seed(seed)
        lattice_parameter = pow(mole_num * mole.get_mass() / density / 0.602, 1 / 3)
        self.lattice_parameter = np.asarray([lattice_parameter] * 3)
        self.box_l = np.zeros(3)
        self.box_h = self.lattice_parameter
        for i in range(mole_num):
            mole.random_coordinate(1, lattice_parameter * np.random.rand(3), param)
            self.add_molecule(mole, aa_num=0)

    def extend_atom3d(self, o_a3d, times, mole_seg, _type, lines=2):
        assert len(times) == 3
        times = np.array(times)
        self.lattice_parameter = o_a3d.lattice_parameter * times
        self.box_l = o_a3d.box_l.copy()
        self.box_h = self.box_l + self.lattice_parameter
        o_a3d.convert_to_default()
        count = 0
        n_atom = o_a3d.get_max_atom_id()
        n_bond = o_a3d.get_max_bond_id()
        n_mole = o_a3d.get_max_mole_id()
        e_mole_seg = dict()
        for i in range(int(times[0])):
            for j in range(int(times[1])):
                for k in range(int(times[2])):
                    multi = np.asarray([i, j, k])
                    atom_reflect = dict()
                    for atom in o_a3d.Atoms.values():
                        n_id = atom.id + count * n_atom
                        n_co = atom.coordinate + o_a3d.lattice_parameter * multi
                        n_m_id = atom.mole_id + count * n_mole
                        na = bsc.Atom(n_id, atom.type, n_co, atom.image, mole_id=n_m_id)
                        atom_reflect[atom.id] = na.id
                        self.append_element(na)
                        if n_m_id not in e_mole_seg:
                            if _type == "x":
                                e_mole_seg[n_m_id] = (j + k * times[1]) * lines + mole_seg[atom.mole_id]
                            else:
                                raise KeyError("type not defined")

                    for bond in o_a3d.Bonds.values():
                        nb = cp.deepcopy(bond)
                        nb.id += count * n_bond
                        nb.atoms = [atom_reflect[_] for _ in bond.atoms]
                        self.append_element(nb)
                    count += 1
        self.cal_all_joints()
        return e_mole_seg

    def renew_coordinate(self, aa3d: AA3D, only_aa=False, only_u=False):
        self.box_l = aa3d.box_l
        self.box_h = aa3d.box_h
        self.lattice_parameter = aa3d.lattice_parameter
        super().renew_coordinate(aa3d, only_aa=only_aa, only_u=only_u)
        self.ce_dict = dict()
        self.default = False
        self.in_cell = False

    def renew_coordinate_file(self, file_path, v=False, extra=None, scaled=True):
        with open(file_path, "r") as f:
            lines = f.readlines()
        box_l = list()
        box_h = list()
        for line in lines[5: 8]:
            arg = line.split()
            box_l.append(float(arg[0]))
            box_h.append(float(arg[1]))
        self.box_l = np.array(box_l)
        self.box_h = np.array(box_h)
        self.lattice_parameter = self.box_h - self.box_l

        extra_list = dict()
        if extra is not None:
            for key in extra:
                extra_list[key] = np.zeros(len(self.Atoms.keys()))

        for line in lines[9:]:
            arg = line.split()
            if arg:
                atom_id = int(arg[0])
                co_s = np.asarray([float(arg[1]), float(arg[2]), float(arg[3])])
                image = np.asarray([int(arg[4]), int(arg[5]), int(arg[6])], dtype=int)
                if v:
                    vs = np.asarray([float(arg[7]), float(arg[8]), float(arg[9])])
                    self.Atoms[atom_id].v = vs
                if extra is not None:
                    for key in extra:
                        extra_list[key][atom_id - 1] = float(arg[key])
                if scaled:
                    self.Atoms[atom_id].coordinate = co_s * self.lattice_parameter + box_l
                else:
                    self.Atoms[atom_id].coordinate = co_s
                self.Atoms[atom_id].image = image
        self.default = False
        self.in_cell = False
        return extra_list

    def convert_to_in_cell(self, only_aa=False):
        if only_aa:
            if not self.in_cell:
                for atom in self.Atoms.values():
                    if atom.type[0].islower():
                        atom.set_in_cell_coordinate(self.lattice_parameter, self.box_l)
                self.in_cell = False
                self.default = False
        else:
            if not self.in_cell:
                for atom in self.Atoms.values():
                    atom.set_in_cell_coordinate(self.lattice_parameter, self.box_l)
                self.in_cell = True
                self.default = False
        pass

    def convert_to_default(self, rebuild=False):
        def default_atom(b_id, f_id):
            # print("called 0")
            if f_id == 0: return
            # print("called")
            b_co = self.Atoms[b_id].coordinate
            f_co = self.Atoms[f_id].coordinate
            d = b_co - f_co
            print("1", np.sum(np.square(d)))
            d -= np.round(d / self.lattice_parameter) * self.lattice_parameter
            print("2" , np.sum(np.square(d)))
            self.Atoms[b_id].coordinate = f_co + d

        if not rebuild:
            if not self.default:
                for atom in self.Atoms.values():
                    atom.set_default_coordinate(self.lattice_parameter)
                self.in_cell = False
                self.default = True
        else:
            visited = set()
            for a_id, atom in self.Atoms.items():
                if atom.type == "TO(1)" or atom.type == "c4":
                    if a_id not in visited:
                        new_visit = self.tranverse_bfs_visit(a_id, default_atom)
                        for v in new_visit:
                            visited.add(v)

    def cal_intra_para(self, para, type_tuple, cal_points):
        self.convert_to_default()
        each_r = np.append(cal_points - (cal_points[1] - cal_points[0]) / 2,
                           cal_points[-1] + (cal_points[1] - cal_points[0]) / 2)
        if para == "Bond":
            cal_func = cso.calculate_distance
        elif para == "Angle":
            cal_func = cso.calculate_angle
        elif para == "Dihedral":
            cal_func = cso.calculate_dihedral_no_sign
        # elif para == "dihedral":
        #     cal_func = cso.calculate_dihedral
        else:
            raise KeyError("No para defined!")
        data_list = list()
        for joint_id in self.__dict__[para + "_type_dict"][type_tuple]:
            joint = self.__dict__[para + "s"][joint_id]
            d = cal_func(*map(lambda atom_id: self.Atoms[atom_id].coordinate, joint.atoms))
            # if para=="Angle" and  d < 10:
            #     print(type_tuple ,d)
            #     for atom_id in joint.atoms:
            #         print(self.Atoms[atom_id])
            data_list.append(d)

        hist, bin_edges = np.histogram(data_list, each_r)
        gr = hist / len(self.__dict__[para + "_type_dict"][type_tuple])
        return gr

    # only extend atoms
    def extend_self(self, extend_range, change_id=True):
        def get_moved_atom(_atom, x, y, z, new_id):
            new_atom = _atom.get_copy()
            d = np.asarray([x, y, z], dtype=int)
            new_atom.coordinate -= d * self.lattice_parameter
            new_atom.image += d
            if change_id:
                new_atom.id = new_id
            return new_atom

        self.convert_to_in_cell()
        upper = self.box_h - extend_range
        lower = self.box_l + extend_range
        m_id = self.get_max_atom_id() + 1
        new_atom_dict = dict()
        for atom in self.Atoms.values():
            record = (atom.coordinate > upper).astype(int) + (atom.coordinate <= lower).astype(int) * -1
            # print(record)
            if record[0] != 0:
                new_atom_dict[m_id] = get_moved_atom(atom, record[0], 0, 0, m_id)
                m_id += 1
                if record[1] != 0:
                    new_atom_dict[m_id] = get_moved_atom(atom, record[0], record[1], 0, m_id)
                    m_id += 1
                    if record[2] != 0:
                        new_atom_dict[m_id] = get_moved_atom(atom, record[0], record[1], record[2], m_id)
                        m_id += 1
                if record[2] != 0:
                    new_atom_dict[m_id] = get_moved_atom(atom, record[0], 0, record[2], m_id)
                    m_id += 1
            if record[1] != 0:
                new_atom_dict[m_id] = get_moved_atom(atom, 0, record[1], 0, m_id)
                m_id += 1
                if record[2] != 0:
                    new_atom_dict[m_id] = get_moved_atom(atom, 0, record[1], record[2], m_id)
                    m_id += 1
            if record[2] != 0:
                new_atom_dict[m_id] = get_moved_atom(atom, 0, 0, record[2], m_id)
                m_id += 1
        self.Atoms.update(new_atom_dict)

    def get_central_and_extended_in_range(self, atom_type, extend_range):
        def get_moved_atom(_atom, x, y, z):
            new_atom = _atom.get_copy()
            d = np.asarray([x, y, z], dtype=int)
            new_atom.coordinate -= d * self.lattice_parameter
            new_atom.image += d
            return new_atom

        if (atom_type, extend_range) not in self.ce_dict:
            self.convert_to_in_cell()
            if atom_type in self.Atom_type_dict:
                central_list = [self.Atoms[_] for _ in self.Atom_type_dict[atom_type]]
            else:
                central_list = list(filter(lambda atom: atom.type == atom_type, self.Atoms.values()))
            extended_list = list()
            # Here used to be an error which influence a lot
            # upper = self.lattice_parameter + self.box_h - extend_range
            upper = self.box_h - extend_range
            lower = self.box_l + extend_range
            for atom in central_list:
                record = (atom.coordinate > upper).astype(int) + (atom.coordinate <= lower).astype(int) * -1
                if record[0] != 0:
                    extended_list.append(get_moved_atom(atom, record[0], 0, 0))
                    if record[1] != 0:
                        extended_list.append(get_moved_atom(atom, record[0], record[1], 0))
                        if record[2] != 0:
                            extended_list.append(get_moved_atom(atom, record[0], record[1], record[2]))
                    if record[2] != 0:
                        extended_list.append(get_moved_atom(atom, record[0], 0, record[2]))
                if record[1] != 0:
                    extended_list.append(get_moved_atom(atom, 0, record[1], 0))
                    if record[2] != 0:
                        extended_list.append(get_moved_atom(atom, 0, record[1], record[2]))
                if record[2] != 0:
                    extended_list.append(get_moved_atom(atom, 0, 0, record[2]))
            self.ce_dict[(atom_type, extend_range)] = (central_list, extended_list)
        return self.ce_dict[(atom_type, extend_range)]

    def cal_inter_para(self, type_tuple, cal_points, exclude=None):
        self.convert_to_in_cell()
        central_a_list, extend_a_list = self.get_central_and_extended_in_range(type_tuple[0], max(cal_points) + 2.)
        central_b_list, extend_b_list = self.get_central_and_extended_in_range(type_tuple[1], max(cal_points) + 2.)
        central_a_co_list = [_.coordinate for _ in central_a_list]
        # print(central_a_co_list)
        central_b_co_list = [_.coordinate for _ in central_b_list]
        extend_b_co_list = [_.coordinate for _ in extend_b_list]

        num_items = len(central_a_list)
        d = cal_points[-1] - cal_points[-2]
        points = np.append(cal_points, np.linspace(cal_points[-1] + d, cal_points[-1] + 200 * d, 200))
        each_r = np.append(points - (points[1] - points[0]) / 2, points[-1] + (points[1] - points[0]) / 2)
        volume = 4 / 3 * np.pi * np.array(
            [np.power(each_r[i + 1], 3) - np.power(each_r[i], 3) for i in range(len(each_r) - 1)])
        each_r_square = np.square(each_r)
        black_list = list()
        hard_types = {"Es", "Ph", "U", "Me", "c1", "n1", "o1", "h1", "c2", "c3"}
        for index_a, central_a in enumerate(central_a_list):
            for index_b, central_b in enumerate(central_b_list + extend_b_list):
                if exclude == 1:
                    if central_b.id in self.get_connected_atom_id(central_a.id):
                        black_list.append((index_a, index_b))
                elif exclude == 2:
                    if central_b.id in self.get_connected_2(central_a.id):
                        black_list.append((index_a, index_b))
                else:
                    if central_a.mole_id == central_b.mole_id:
                        if type_tuple[0] in hard_types and type_tuple[1] in hard_types:
                            black_list.append((index_a, index_b))
                        elif central_b.id in self.get_connected_atom_id(central_a.id):
                            # print(self.get_connected_atom_id(central_a.id))
                            black_list.append((index_a, index_b))
        sum_dif1 = cso.matrix_distance_range(central_a_co_list, central_b_co_list + extend_b_co_list,
                                             (each_r[0], each_r[-1]), black_list)
        hist1, bin_edges1 = np.histogram(sum_dif1, each_r_square)
        # black_list = list()
        # for index_A, central_A in enumerate(central_A_list):
        #     for index_B, extended_B in enumerate(extended_B_list):
        #         if central_A.mole_id[0] == extended_B.mole_id[0] and abs(
        #                 central_A.mole_id[1] - extended_B.mole_id[
        #                     1]) <= excluded_range and central_A.period == extended_B.period:
        #             black_list.append((index_A, index_B))
        # sum_dif2 = cso.matrix_distance_range(central_A_centroid_list, extended_B_centroid_list, (each_r[0], each_r[-1]),
        #                                      black_list)
        # hist2, bin_edges2 = np.histogram(sum_dif2, each_r_square)
        # bandwidth = 0.2
        # A_B_count = np.array(hist1) + np.array(hist2)
        gr = hist1 / float(num_items) / volume
        return gr

    # only atom
    def get_fake_copy(self, only_aa=False, only_u=False, only_type=None):
        new_atom3d = Atom3D(self.box_l, self.box_h)
        for atom in self.Atoms.values():
            if only_aa and not atom.type[0].islower():
                continue
            if only_u and atom.type not in ["c1", "U"]:
                continue
            if only_type and atom.type not in only_type:
                continue
            new_atom3d.append_element(atom.get_copy())
        return new_atom3d

    # # there is something wrong
    # def cg_to_fgh(self, u_orien=None):
    #     def append_fg_atom(_ru, _cent, end1, end2, _c_head, _c_end, _m_id, _type):
    #         _ru.set_centroid(_cent)
    #         _ru.rotate_to_vector(end1 - end2)
    #         if _type == "U":
    #             old_orient = _ru.Atoms[4].coordinate - _ru.Atoms[3].coordinate
    #         else:
    #             old_orient = _ru.Atoms[2].coordinate + _ru.Atoms[3].coordinate - _ru.Atoms[5].coordinate - _ru.Atoms[6].coordinate
    #         old_orient /= cso.calculate_distance(old_orient)

    #         if u_orien is not None:
    #             new_orient = np.asarray(u_orien)
    #         else:
    #             new_orient = (end1 + end2) / 2 - _cent
    #         _ru.rotate_orientation(old_orient, new_orient)

    #         for r_atom in _ru.Atoms.values():
    #             new_atom = bsc.Atom(r_atom.id + ma, r_atom.type, r_atom.coordinate.copy(), mole_id=_m_id)
    #             self.append_element(new_atom)
    #         for r_bond in _ru.Bonds.values():
    #             new_bond = bsc.Bond(r_bond.id + mb, r_bond.atoms[0] + ma, r_bond.atoms[1] + ma, _type=r_bond.type)
    #             self.append_element(new_bond)
    #         _type, re = sort_types((self.Atoms[_c_head].type, _ru.Atoms[_ru.head].type))
    #         if re:
    #             head_bond = bsc.Bond(len(_ru.Bonds) + mb + 1, _ru.head + ma, _c_head, _type=_type)
    #         else:
    #             head_bond = bsc.Bond(len(_ru.Bonds) + mb + 1, _c_head, _ru.head + ma, _type=_type)
    #         self.append_element(head_bond)
    #         _type, re = sort_types((self.Atoms[_c_end].type, _ru.Atoms[_ru.end].type))
    #         if re:
    #             end_bond = bsc.Bond(len(_ru.Bonds) + mb + 2, _ru.end + ma, _c_end, _type=_type)
    #         else:
    #             end_bond = bsc.Bond(len(_ru.Bonds) + mb + 2, _c_end, _ru.end + ma, _type=_type)
    #         self.append_element(end_bond)

    #     self.convert_to_default()
    #     with open("/home/centos/CGMD/template/Uh.pkl", "rb") as f:
    #         U = pkl.load(f)
    #     with open("/home/centos/CGMD/template/Ph.pkl", "rb") as f:
    #         Ph = pkl.load(f)

    #     # print(U)
    #     U.rotate_n_random()
    #     Ph.rotate_n_random()
    #     self.clear_all_joints()

    #     bb_list = list()
    #     atom_ids = list(self.Atoms.keys())
    #     atom_ids.sort()
    #     for atom_id in atom_ids:
    #         cg_atom = self.Atoms[atom_id]
    #         if cg_atom.type in ["U", "Ph"]:
    #             head_end = self.get_connected_atom_id(cg_atom.id)
    #             chain_head = head_end.pop()
    #             chain_end = head_end.pop()
    #             if chain_head < chain_end:
    #                 chain_head, chain_end = chain_end, chain_head
    #             centroid = cg_atom.coordinate
    #             self.delete_atom(atom_id)
    #             # vec =
    #             if cg_atom.type == "U":
    #                 repeat_unit = cp.deepcopy(U)
    #             else:
    #                 repeat_unit = cp.deepcopy(Ph)
    #             ma = self.get_max_atom_id()
    #             mb = self.get_max_bond_id()
    #             append_fg_atom(repeat_unit, centroid, self.Atoms[chain_end].coordinate,
    #                            self.Atoms[chain_head].coordinate, chain_head, chain_end, cg_atom.mole_id, cg_atom.type)
    #     # print("1")
    #     self.sort_atom_bond_id()
    #     # print("2")
    #     self.cal_all_joints(improper="ignore")

    def cgu_to_fgh(self, u_orien=None):
        def append_fg_atom(_ru, _cent, end1, end2, _c_head, _c_end, _m_id, _type):
            _ru.set_centroid(_cent)
            _ru.rotate_to_vector(end1 - end2)
            if _type == "U":
                old_orient = _ru.Atoms[4].coordinate - _ru.Atoms[3].coordinate
            else:
                old_orient = _ru.Atoms[2].coordinate + _ru.Atoms[3].coordinate - _ru.Atoms[5].coordinate - _ru.Atoms[6].coordinate
            old_orient /= cso.calculate_distance(old_orient)

            if u_orien is not None:
                new_orient = np.asarray(u_orien)
            else:
                new_orient = (end1 + end2) / 2 - _cent
            # _ru.rotate_orientation(old_orient, new_orient)

            for r_atom in _ru.Atoms.values():
                new_atom = bsc.Atom(r_atom.id + ma, r_atom.type, r_atom.coordinate.copy(), mole_id=_m_id)
                self.append_element(new_atom)
            for r_bond in _ru.Bonds.values():
                new_bond = bsc.Bond(r_bond.id + mb, r_bond.atoms[0] + ma, r_bond.atoms[1] + ma, _type=r_bond.type)
                self.append_element(new_bond)
            _type, re = sort_types((self.Atoms[_c_head].type, _ru.Atoms[_ru.head].type))
            if re:
                head_bond = bsc.Bond(len(_ru.Bonds) + mb + 1, _ru.head + ma, _c_head, _type=_type)
            else:
                head_bond = bsc.Bond(len(_ru.Bonds) + mb + 1, _c_head, _ru.head + ma, _type=_type)
            self.append_element(head_bond)
            _type, re = sort_types((self.Atoms[_c_end].type, _ru.Atoms[_ru.end].type))
            if re:
                end_bond = bsc.Bond(len(_ru.Bonds) + mb + 2, _ru.end + ma, _c_end, _type=_type)
            else:
                end_bond = bsc.Bond(len(_ru.Bonds) + mb + 2, _c_end, _ru.end + ma, _type=_type)
            self.append_element(end_bond)
        
        def append_fg_atom_tail(_ru, _cent, end1, end2, _c_head, _m_id, _type):
            _ru.set_centroid(_cent)
            _ru.rotate_to_vector(end1 - end2)
            if _type == "U":
                old_orient = _ru.Atoms[4].coordinate - _ru.Atoms[3].coordinate
            else:
                old_orient = _ru.Atoms[2].coordinate + _ru.Atoms[3].coordinate - _ru.Atoms[5].coordinate - _ru.Atoms[6].coordinate
            old_orient /= cso.calculate_distance(old_orient)

            for r_atom in _ru.Atoms.values():
                new_atom = bsc.Atom(r_atom.id + ma, r_atom.type, r_atom.coordinate.copy(), mole_id=_m_id)
                self.append_element(new_atom)
            for r_bond in _ru.Bonds.values():
                new_bond = bsc.Bond(r_bond.id + mb, r_bond.atoms[0] + ma, r_bond.atoms[1] + ma, _type=r_bond.type)
                self.append_element(new_bond)
            _type, re = sort_types((self.Atoms[_c_head].type, _ru.Atoms[_ru.head].type))
            if re:
                head_bond = bsc.Bond(len(_ru.Bonds) + mb + 1, _ru.head + ma, _c_head, _type=_type)
            else:
                head_bond = bsc.Bond(len(_ru.Bonds) + mb + 1, _c_head, _ru.head + ma, _type=_type)
            self.append_element(head_bond)

        self.convert_to_default(rebuild=True)
        with open("/home/centos/CGMD/template/U1.pkl", "rb") as f:
            U = pkl.load(f)
        with open("/home/centos/CGMD/template/Ph.pkl", "rb") as f:
            Ph = pkl.load(f)
        U.set_type_define("fg")
        Ph.set_type_define("fg")
        # print(U)
        # U.rotate_n_random()
        # Ph.rotate_n_random()
        self.clear_all_joints()

        bb_list = list()
        atom_ids = list(self.Atoms.keys())
        atom_ids.sort()
        for atom_id in atom_ids:
            cg_atom = self.Atoms[atom_id]
            if cg_atom.type in ["Ph"]:
                head_end = self.get_connected_atom_id(cg_atom.id)
                if len(head_end) == 2:
                    chain_head = head_end.pop()
                    chain_end = head_end.pop()
                    # if chain_head < chain_end:
                    #     chain_head, chain_end = chain_end, chain_head
                    centroid = cg_atom.coordinate
                    self.delete_atom(atom_id)
                    # vec =
                    if cg_atom.type == "U":
                        repeat_unit = cp.deepcopy(U)
                    else:
                        repeat_unit = cp.deepcopy(Ph)
                    ma = self.get_max_atom_id()
                    mb = self.get_max_bond_id()
                    append_fg_atom(repeat_unit, centroid, self.Atoms[chain_end].coordinate,
                                self.Atoms[chain_head].coordinate, chain_head, chain_end, cg_atom.mole_id, cg_atom.type)
                elif len(head_end) == 1:
                    print("tail fg")
                    chain_head = head_end.pop()
                    centroid = cg_atom.coordinate
                    cent = self.delete_atom(atom_id)
                    # vec =
                    if cg_atom.type == "U":
                        repeat_unit = cp.deepcopy(U)
                    else:
                        repeat_unit = cp.deepcopy(Ph)
                    ma = self.get_max_atom_id()
                    mb = self.get_max_bond_id()
                    append_fg_atom_tail(repeat_unit, centroid, cent.coordinate,
                                self.Atoms[chain_head].coordinate, chain_head, cg_atom.mole_id, cg_atom.type)
                else:
                    raise KeyError("too many connected: %d" % (len(head_end)))
        # print("1")
        self.sort_atom_bond_id()
        # print("2")
        self.cal_all_joints(improper="ignore")

    def cg_to_cgu(self, u_orien=None):
        def append_fg_atom(_ru, _cent, end1, end2, _c_head, _c_end, _m_id):
            _ru.set_centroid(_cent)
            _ru.rotate_to_vector(end1 - end2)
            old_orient = _ru.Atoms[4].coordinate - _ru.Atoms[3].coordinate
            old_orient /= cso.calculate_distance(old_orient)

            if u_orien is not None:
                new_orient = np.asarray(u_orien)
            else:
                new_orient = (end1 + end2) / 2 - _cent
            _ru.rotate_orientation(old_orient, new_orient)

            for r_atom in _ru.Atoms.values():
                new_atom = bsc.Atom(r_atom.id + ma, r_atom.type, r_atom.coordinate.copy(), mole_id=_m_id)
                self.append_element(new_atom)
            for r_bond in _ru.Bonds.values():
                new_bond = bsc.Bond(r_bond.id + mb, r_bond.atoms[0] + ma, r_bond.atoms[1] + ma, _type=r_bond.type)
                self.append_element(new_bond)
            _type, re = sort_types((self.Atoms[_c_head].type, _ru.Atoms[_ru.head].type))
            if re:
                head_bond = bsc.Bond(len(_ru.Bonds) + mb + 1, _ru.head + ma, _c_head, _type=_type)
            else:
                head_bond = bsc.Bond(len(_ru.Bonds) + mb + 1, _c_head, _ru.head + ma, _type=_type)
            self.append_element(head_bond)
            _type, re = sort_types((self.Atoms[_c_end].type, _ru.Atoms[_ru.end].type))
            if re:
                end_bond = bsc.Bond(len(_ru.Bonds) + mb + 2, _ru.end + ma, _c_end, _type=_type)
            else:
                end_bond = bsc.Bond(len(_ru.Bonds) + mb + 2, _c_end, _ru.end + ma, _type=_type)
            self.append_element(end_bond)

        self.convert_to_default()
        with open("/home/centos/CGMD/template/U1.pkl", "rb") as f:
            U = pkl.load(f)

        # print(U)
        U.rotate_n_random()
        self.clear_all_joints()

        bb_list = list()
        atom_ids = list(self.Atoms.keys())
        atom_ids.sort()
        u_num = 0
        for atom_id in atom_ids:
            cg_atom = self.Atoms[atom_id]
            if cg_atom.type == "U":
                u_num += 1
                if u_num % 100 == 0: print(u_num)
                head_end = self.get_connected_atom_id(cg_atom.id)
                chain_head = head_end.pop()
                chain_end = head_end.pop()
                if chain_head < chain_end:
                    chain_head, chain_end = chain_end, chain_head
                centroid = cg_atom.coordinate
                self.delete_atom(atom_id)
                # vec =
                repeat_unit = cp.deepcopy(U)
                ma = self.get_max_atom_id()
                mb = self.get_max_bond_id()
                append_fg_atom(repeat_unit, centroid, self.Atoms[chain_end].coordinate,
                               self.Atoms[chain_head].coordinate, chain_head, chain_end, cg_atom.mole_id)
        print("1")
        self.sort_atom_bond_id()
        print("2")
        self.cal_all_joints(improper="ignore")

    def cg_to_aa(self):

        def append_cg_atoms(_atom, _aa3d, _ends, _tail_record, _ma, _mb):
            unit = cp.deepcopy(aa_unit[_atom.type])
            unit.set_centroid(_atom.coordinate)
            if len(_ends) == 1:
                end = _ends.pop()
                new_vec = self.Atoms[end].coordinate - _atom.coordinate
                unit.rotate_to_vector(new_vec)
                _tail_record[(_atom.id, end)] = unit.end + _ma
            else:
                head = _ends.pop()
                end = _ends.pop()
                if _atom.type == "Es":
                    if self.Atoms[head].type != "Ph":
                        head, end = end, head
                new_vec = self.Atoms[end].coordinate - self.Atoms[head].coordinate
                if _atom.type != "Me":
                    unit.rotate_to_vector(new_vec)
                if _atom.type == "U":
                    old_orient = unit.Atoms[4].coordinate - unit.Atoms[3].coordinate
                    old_orient /= cso.calculate_distance(old_orient)
                    new_orient = (self.Atoms[end].coordinate + self.Atoms[head].coordinate) / 2 - _atom.coordinate
                    unit.rotate_orientation(old_orient, new_orient)


                _tail_record[(_atom.id, head)] = unit.head + _ma
                _tail_record[(_atom.id, end)] = unit.end + _ma
            for atom in unit.Atoms.values():
                new_atom = bsc.Atom(atom.id + _ma, atom.type, atom.coordinate.copy(), mole_id=_atom.mole_id)
                _aa3d.append_element(new_atom)
            for bond in unit.Bonds.values():
                new_bond = bsc.Bond(bond.id + mb, bond.atoms[0] + _ma, bond.atoms[1] + _ma, _type=bond.type)
                _aa3d.append_element(new_bond)
            return _ma + len(unit.Atoms), _mb + len(unit.Bonds)

        def append_single_atom(_atom, _aa3d, _aa_map, _ma):
            _ma += 1
            new_atom = bsc.Atom(_ma, _atom.type, _atom.coordinate.copy(), mole_id=_atom.mole_id)
            _aa3d.append_element(new_atom)
            _aa_map[_atom.id] = _ma
            return _ma

        cg_types = ["U", "Ph", "TO(1)", "TO(2)", "Es", "Me"]
        aa_unit = dict()
        for ctp in cg_types:
            with open("/home/centos/CGMD/template/aa_%s.pkl" % ctp, "rb") as f:
                unit = pkl.load(f)
            unit.set_type_define("aa")
            aa_unit[ctp] = unit

        self.convert_to_default()
        aa_3d = Atom3D(box_l=self.box_l, box_h=self.box_h)
        aa_3d.set_type_define("aa")
        # (origin_id, origin_connected_id): new_tail_id
        tail_record = dict()
        aa_map = dict()

        # if cg bead, convert to aa; else, copy it.
        # record tail and aa_map
        ma = 0
        mb = 0
        for atom in self.Atoms.values():
            if atom.type in cg_types:
                ends = self.get_connected_atom_id(atom.id)
                ma, mb = append_cg_atoms(atom, aa_3d, ends, tail_record, ma, mb)
            else:
                ma = append_single_atom(atom, aa_3d, aa_map, ma)

        # connect bead ends and atoms
        for bond in self.Bonds.values():
            a0 = bond.atoms[0]
            a1 = bond.atoms[1]
            if self.Atoms[a0].type in cg_types:
                new_a0 = tail_record[(a0, a1)]
            else:
                new_a0 = aa_map[a0]
            if self.Atoms[a1].type in cg_types:
                new_a1 = tail_record[(a1, a0)]
            else:
                new_a1 = aa_map[a1]
            mb = aa_3d.create_bond(new_a0, new_a1, mb)

        print("1")
        aa_3d.sort_atom_bond_id()
        print("2")
        aa_3d.cal_all_joints(improper="COMPASS")
        return aa_3d


    def cgu_to_aa(self):

        def append_cg_atoms(_atom, _aa3d, _ends, _tail_record, _ma, _mb):
            unit = aa_unit[_atom.type]
            unit.set_centroid(_atom.coordinate)
            if len(_ends) == 1:
                end = _ends.pop()
                new_vec = self.Atoms[end].coordinate - _atom.coordinate
                unit.rotate_to_vector(new_vec)
                _tail_record[(_atom.id, end)] = unit.end + _ma
            else:
                head = _ends.pop()
                end = _ends.pop()
                if _atom.type == "Es":
                    if self.Atoms[head].type != "Ph":
                        head, end = end, head
                new_vec = self.Atoms[end].coordinate - self.Atoms[head].coordinate
                if _atom.type != "Me":
                    unit.rotate_to_vector(new_vec)
                if _atom.type == "U":
                    old_orient = unit.Atoms[4].coordinate - unit.Atoms[3].coordinate
                    old_orient /= cso.calculate_distance(old_orient)
                    new_orient = (self.Atoms[end].coordinate + self.Atoms[head].coordinate) / 2 - _atom.coordinate
                    unit.rotate_orientation(old_orient, new_orient)


                _tail_record[(_atom.id, head)] = unit.head + _ma
                _tail_record[(_atom.id, end)] = unit.end + _ma
            for atom in unit.Atoms.values():
                new_atom = bsc.Atom(atom.id + _ma, atom.type, atom.coordinate.copy(), mole_id=_atom.mole_id)
                _aa3d.append_element(new_atom)
            for bond in unit.Bonds.values():
                new_bond = bsc.Bond(bond.id + mb, bond.atoms[0] + _ma, bond.atoms[1] + _ma, _type=bond.type)
                _aa3d.append_element(new_bond)
            return _ma + len(unit.Atoms), _mb + len(unit.Bonds)

        def append_single_atom(_atom, _aa3d, _aa_map, _ma):
            _ma += 1
            new_atom = bsc.Atom(_ma, cgu_type_map[_atom.type], _atom.coordinate.copy(), mole_id=_atom.mole_id)
            _aa3d.append_element(new_atom)
            _aa_map[_atom.id] = _ma
            return _ma

        cg_types = ["U", "Ph", "TO(1)", "TO(2)", "Es", "Me"]
        cgu_type_map = {"c1": "c3\"", "n1" : "n3mh", "h1" : "h1n", "o1" : "o1="}
        aa_unit = dict()
        for ctp in cg_types:
            with open("/home/centos/CGMD/template/aa_%s.pkl" % ctp, "rb") as f:
                unit = pkl.load(f)
            unit.set_type_define("aa")
            aa_unit[ctp] = unit

        self.convert_to_default()
        aa_3d = Atom3D(box_l=self.box_l, box_h=self.box_h)
        aa_3d.set_type_define("aa")
        # (origin_id, origin_connected_id): new_tail_id
        tail_record = dict()
        aa_map = dict()

        # if cg bead, convert to aa; else, copy it.
        # record tail and aa_map
        ma = 0
        mb = 0
        for atom in self.Atoms.values():
            if atom.type in cg_types:
                ends = self.get_connected_atom_id(atom.id)
                ma, mb = append_cg_atoms(atom, aa_3d, ends, tail_record, ma, mb)
            else:
                ma = append_single_atom(atom, aa_3d, aa_map, ma)

        # connect bead ends and atoms
        for bond in self.Bonds.values():
            a0 = bond.atoms[0]
            a1 = bond.atoms[1]
            if self.Atoms[a0].type in cg_types:
                new_a0 = tail_record[(a0, a1)]
            else:
                new_a0 = aa_map[a0]
            if self.Atoms[a1].type in cg_types:
                new_a1 = tail_record[(a1, a0)]
            else:
                new_a1 = aa_map[a1]
            mb = aa_3d.create_bond(new_a0, new_a1, mb)

        print("1")
        aa_3d.sort_atom_bond_id()
        print("2")
        aa_3d.cal_all_joints(improper="COMPASS")
        return aa_3d


    def cgu_to_aa_draw(self):

        def append_cg_atoms(_atom, _aa3d, _ends, _tail_record, _ma, _mb):
            unit = aa_unit[_atom.type]
            unit.set_centroid(_atom.coordinate)
            if len(_ends) == 1:
                end = _ends.pop()
                new_vec = self.Atoms[end].coordinate - _atom.coordinate
                unit.rotate_to_vector(new_vec)
                _tail_record[(_atom.id, end)] = unit.end + _ma
            else:
                head = _ends.pop()
                end = _ends.pop()
                if _atom.type == "Es":
                    if self.Atoms[head].type != "Ph":
                        head, end = end, head
                new_vec = self.Atoms[end].coordinate - self.Atoms[head].coordinate
                if _atom.type != "Me":
                    unit.rotate_to_vector(new_vec)
                if _atom.type == "U":
                    old_orient = unit.Atoms[4].coordinate - unit.Atoms[3].coordinate
                    old_orient /= cso.calculate_distance(old_orient)
                    new_orient = (self.Atoms[end].coordinate + self.Atoms[head].coordinate) / 2 - _atom.coordinate
                    unit.rotate_orientation(old_orient, new_orient)


                _tail_record[(_atom.id, head)] = unit.head + _ma
                _tail_record[(_atom.id, end)] = unit.end + _ma
            for atom in unit.Atoms.values():
                new_atom = bsc.Atom(atom.id + _ma, atom.type, atom.coordinate.copy(), mole_id=_atom.mole_id)
                _aa3d.append_element(new_atom)
            for bond in unit.Bonds.values():
                new_bond = bsc.Bond(bond.id + mb, bond.atoms[0] + _ma, bond.atoms[1] + _ma, _type=bond.type)
                _aa3d.append_element(new_bond)
            return _ma + len(unit.Atoms), _mb + len(unit.Bonds)

        def append_single_atom(_atom, _aa3d, _aa_map, _ma):
            _ma += 1
            new_atom = bsc.Atom(_ma, cgu_type_map[_atom.type], _atom.coordinate.copy(), mole_id=_atom.mole_id)
            _aa3d.append_element(new_atom)
            _aa_map[_atom.id] = _ma
            return _ma

        cg_types = ["U", "Ph", "TO(1)", "TO(2)", "Es", "Me"]
        hard_types = ["U", "Ph", "Me"]
        cgu_type_map = {"c1": "h", "n1" : "h", "h1" : "h", "o1" : "h"}
        aa_unit = dict()
        for ctp in cg_types:
            with open("/home/centos/CGMD/template/aa_%s.pkl" % ctp, "rb") as f:
                unit = pkl.load(f)
                if ctp in hard_types:
                    for atom in unit.Atoms.values():
                        atom.type = "h"
                else:
                    for atom in unit.Atoms.values():
                        atom.type = "s"

            unit.set_type_define("aa")
            aa_unit[ctp] = unit

        self.convert_to_default()
        aa_3d = Atom3D(box_l=self.box_l, box_h=self.box_h)
        aa_3d.set_type_define("aa")
        # (origin_id, origin_connected_id): new_tail_id
        tail_record = dict()
        aa_map = dict()

        # if cg bead, convert to aa; else, copy it.
        # record tail and aa_map
        ma = 0
        mb = 0
        for atom in self.Atoms.values():
            if atom.type in cg_types:
                ends = self.get_connected_atom_id(atom.id)
                ma, mb = append_cg_atoms(atom, aa_3d, ends, tail_record, ma, mb)
            else:
                ma = append_single_atom(atom, aa_3d, aa_map, ma)

        # connect bead ends and atoms
        for bond in self.Bonds.values():
            a0 = bond.atoms[0]
            a1 = bond.atoms[1]
            if self.Atoms[a0].type in cg_types:
                new_a0 = tail_record[(a0, a1)]
            else:
                new_a0 = aa_map[a0]
            if self.Atoms[a1].type in cg_types:
                new_a1 = tail_record[(a1, a0)]
            else:
                new_a1 = aa_map[a1]
            mb = aa_3d.create_bond(new_a0, new_a1, mb)

        print("1")
        aa_3d.sort_atom_bond_id()
        return aa_3d

    # c1: 12, c3: 14, end connect wrong
    def convert_to_fg(self, h=True):
        def append_fg_atom(_ru, _cent, _vec, _c_head, _c_end, _m_id):
            _ru.set_centroid(_cent)
            _ru.rotate_to_vector(_vec)
            for r_atom in _ru.Atoms.values():
                new_atom = bsc.Atom(r_atom.id + ma, r_atom.type, r_atom.coordinate.copy(), mole_id=_m_id)
                self.append_element(new_atom)
            for r_bond in _ru.Bonds.values():
                new_bond = bsc.Bond(r_bond.id + mb, r_bond.atoms[0] + ma, r_bond.atoms[1] + ma, _type=r_bond.type)
                self.append_element(new_bond)
            _type, re = sort_types((self.Atoms[_c_head].type, _ru.Atoms[_ru.head].type))
            if re:
                head_bond = bsc.Bond(len(_ru.Bonds) + mb + 1, _ru.head + ma, _c_head, _type=_type)
            else:
                head_bond = bsc.Bond(len(_ru.Bonds) + mb + 1, _c_head, _ru.head + ma, _type=_type)
            self.append_element(head_bond)
            _type, re = sort_types((self.Atoms[_c_end].type, _ru.Atoms[_ru.end].type))
            if re:
                end_bond = bsc.Bond(len(_ru.Bonds) + mb + 2, _ru.end + ma, _c_end, _type=_type)
            else:
                end_bond = bsc.Bond(len(_ru.Bonds) + mb + 2, _c_end, _ru.end + ma, _type=_type)
            self.append_element(end_bond)

        self.convert_to_default(rebuild=True)
        if h:
            # with open("/home/centos/CGMD/template/Uh.pkl", "rb") as f:
            
            # c1
            with open("/home/centos/CGMD/template/U1.pkl", "rb") as f:
                U = pkl.load(f)
        else:
            with open("/home/centos/CGMD/template/U.pkl", "rb") as f:
                U = pkl.load(f)
        with open("/home/centos/CGMD/template/Ph.pkl", "rb") as f:
            Ph = pkl.load(f)
        self.clear_all_joints()
        # change_list = self.Atom_type_dict["Ph(0)"] + self.Atom_type_dict["U(0)"] + self.Atom_type_dict["B(0)"]
        last_type = ""
        bb_list = list()
        bb_head_end = list()
        atom_ids = list(self.Atoms.keys())
        atom_ids.sort()
        for atom_id in atom_ids:
            cg_atom = self.Atoms[atom_id]
            # if cg_atom.type in ["Ph(0)", "U(0)", "B(0)"]:
            if cg_atom.type in ["U", "Ph", "B"]:
                head_end = self.get_connected_atom_id(cg_atom.id)
                if len(head_end) < 2:
                    print(cg_atom.id,cg_atom.type)
                    chain_head = head_end.pop()
                    print(self.Atoms[chain_head].id, self.Atoms[chain_head].type)
                chain_head = head_end.pop()
                chain_end = head_end.pop()
                # if chain_head < chain_end:
                #     chain_head, chain_end = chain_end, chain_head
                centroid = cg_atom.coordinate
                self.delete_atom(atom_id)
                vec = self.Atoms[chain_end].coordinate - self.Atoms[chain_head].coordinate
                if cg_atom.type == "Ph" or cg_atom.type == "B(0)":
                    repeat_unit = Ph
                else:
                    repeat_unit = U
                ma = self.get_max_atom_id()
                mb = self.get_max_bond_id()
                append_fg_atom(repeat_unit, centroid, vec, chain_head, chain_end, cg_atom.mole_id)
                if cg_atom.type == "B(0)":
                    if last_type == "B(0)":
                        bb_head_end.append(repeat_unit.head + ma)
                        bb_list.append(cp.copy(bb_head_end))
                        bb_head_end = list()
                    else:
                        bb_head_end.append(repeat_unit.end + ma)
            last_type = cg_atom.type
        for bb_h_e in bb_list:
            c_h = bb_h_e[0]
            c_e = bb_h_e[1]
            b_id = 0
            for bond_id in self.Atoms[c_h].Bonds:
                if self.Bonds[bond_id].get_other_atom_id(c_h) == c_e:
                    b_id = bond_id
            self.delete_bond(b_id)
            co = (self.Atoms[c_h].coordinate + self.Atoms[c_e].coordinate) / 2
            ma = self.get_max_atom_id()
            mb = self.get_max_bond_id()
            new_atom = bsc.Atom(ma + 1, "c1", co, mole_id=self.Atoms[c_h].mole_id)
            self.append_element(new_atom)
            b1 = bsc.Bond(mb + 1, ma + 1, c_h, _type=("c1", "c2"))
            b2 = bsc.Bond(mb + 2, ma + 1, c_e, _type=("c1", "c2"))
            self.append_element(b1)
            self.append_element(b2)
        print("1")
        self.sort_atom_bond_id()
        print("2")
        self.cal_all_joints(improper="ignore")

    def cgu_to_cg(self):
        def get_u_ids(c1_id):
            u_set = {c1_id}
            ph_set = list()
            connects = self.get_connected_atom_id(c1_id)
            for a_id in connects:
                u_set.add(a_id)
                a_con = self.get_connected_atom_id(a_id)
                for ac in a_con:
                    if self.Atoms[ac].type == "Ph":
                        ph_set.append(ac)
                    else:
                        u_set.add(ac)
            return u_set, ph_set

        id_map = dict()
        self.convert_to_default()
        atom_id = self.get_max_atom_id() + 1
        bond_id = self.get_max_bond_id() + 1
        self.uss = dict()
        for am_id in list(self.Atoms.keys()):
            if am_id not in self.Atoms:
                continue
            atom = self.Atoms[am_id]
            if atom.type == "c1":
                u_s, ph_s = get_u_ids(am_id)
                # print(u_s, ph_s)
                ph_s.sort()
                u = bsc.Atom(atom_id, "U", atom.coordinate.copy(), atom.image.copy(), mole_id=atom.mole_id)
                self.append_element(u)
                b1 = bsc.Bond(bond_id, ph_s[0], atom_id, _type=("Ph", "U"))
                bond_id += 1
                b2 = bsc.Bond(bond_id, ph_s[1], atom_id, _type=("Ph", "U"))
                bond_id += 1
                self.append_element(b1)
                self.append_element(b2)
                self.uss[atom_id] = u_s
                id_map[atom_id] = u_s
                atom_id += 1
                for a_id in u_s:
                    self.delete_atom(a_id)
        self.clear_all_joints()
        atom_reflect = self.sort_atom_bond_id()
        new_id_map = dict()
        for old_id, new_id in atom_reflect.items():
            if old_id in id_map:
                new_id_map[new_id] = id_map[old_id]
            else:
                new_id_map[new_id] = old_id
        self.cal_all_joints()
        return new_id_map

    def renew_from_cgu(self, atom3d, only_u=False):
        def cal_centroid(ids):
            mass_sum = 0.
            cen = np.zeros(3)
            for idx in ids:
                atom = atom3d.Atoms[idx]
                mass = self.get_atom_type_mass(atom.type)
                mass_sum += mass
                cen += mass * atom.coordinate
            return cen / mass_sum

        if not hasattr(self, "uss"):
            raise KeyError("convert to cg first")
        atom3d.convert_to_default()
        self.lattice_parameter = atom3d.lattice_parameter
        self.box_l = atom3d.box_l
        self.box_h = atom3d.box_h
        # print(atom3d)
        for k, v in self.uss.items():
            atom = self.Atoms[k]
            atom.image = np.zeros(3, dtype=int)
            atom.coordinate = cal_centroid(v)
        if only_u:
            return 0
        else:
            for k, atom in self.Atoms.items():
                if k not in atom3d.Atoms:
                    if atom.type != "U":
                        print(self.Atoms[k])
                    continue
                # if atom.type != "U":
                atom.coordinate = atom3d.Atoms[k].coordinate.copy()
                atom.image = np.zeros(3, dtype=int)
        self.default = True
        self.in_cell = False
        pass

    def exchange_axis(self, axis1, axis2, v=True):
        def exchange_value(v):
            temp1, temp2 = v[axis1], v[axis2]
            v[axis1] = temp2
            v[axis2] = temp1
        
        exchange_value(self.lattice_parameter)
        exchange_value(self.box_h)
        exchange_value(self.box_l)

        for atom in self.Atoms.values():
            exchange_value(atom.coordinate)
            if v:
                exchange_value(atom.v)


def create_atom3d(root_path, mole_num):
    m1 = Molecule()
    m1.create_from_info(root_path + "struct_info.pkl")
    atom3d = Atom3D()
    atom3d.add_molecules(m1, mole_num)

    aa = AA3D()
    aa.get_mass(root_path + "1blk_40.data")
    aa.renew_coordinate(root_path + "1blk_40.lammpstrj", step=int(15000000))
    atom3d.renew_coordinate(aa)
    return atom3d


def create_atom3d_cg(root_path, mole_num):
    m1 = Molecule()
    m1.create_from_info(root_path + "cg_struct_info.pkl")
    atom3d = Atom3D()
    atom3d.add_molecules(m1, mole_num)

    aa = AA3D()
    aa.get_mass(root_path + "1blk_40.data")
    aa.renew_coordinate(root_path + "1blk_40.lammpstrj", step=int(15000000))
    atom3d.renew_coordinate(aa)
    return atom3d


if __name__ == '__main__':
    fp = "/home/centos/model/aa/SH6S/0/"

    atom3d = Atom3D()
    atom3d.create_from_info(fp + "fg_struct_info.pkl")
    aa = AA3D()
    aa.get_mass(fp + "SH6S.data")
    # aa.renew_coordinate(fp + "SH6S.lammpstrj", step=int(24005000))
    aa.renew_coordinate(fp + "SH6S.lammpstrj", step=int(16000000))
    atom3d.renew_coordinate(aa)
    atom3d.convert_to_in_cell()

    # # atom3d = create_atom3d(fp, 40)
    # atom3d = create_atom3d_cg(fp, 1)
    # print(atom3d)
    # atom3d.convert_to_fg()

    from pyMD.file_parser import LmpParser

    LmpParser.create_data_file(atom3d, "SH6S.data")
