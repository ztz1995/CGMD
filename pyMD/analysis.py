from pyMD.file_parser import LmpParser
from pyMD import functions as func
from pyMD import collective_structure_class as csc
from pyMD import basic_structure_class as bsc
from pyMD import cartesian_operation as cso
from pyMD import file_parser as fps
from scipy.spatial import cKDTree
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.signal import convolve2d
import copy as cp
from multiprocessing import Pool, Manager, Process
import os
from typing import Dict


# only contain atom type and coordinate information if fake
def atom3d_generator(data_file_path, trj_file_path, start, step, frames, one_file=False, _type="fg", fake=True):
    parser = LmpParser()
    atom3d = parser.load_data_file(data_file_path, _type=_type)
    frame = 0
    while frame < frames:
        atom3d = parser.renew_coordinate_file(atom3d, trj_file_path, start + frame * step, one_file=one_file)
        if fake:
            yield atom3d.get_fake_copy()
        else:
            yield atom3d
        frame += 1


# only contain atom type and coordinate information if fake
def atom3d_generator_from_aa(info_path, data_file_path, trj_file_path, start, step, frames, one_file=True, fake=True,
                             renew_v=False):
    atom3d = csc.Atom3D()
    atom3d.create_from_info(info_path)
    aa = csc.AA3D()
    aa.get_mass(data_file_path)
    frame = 0
    begin = 0
    lines = None
    while frame < frames:
        # print(frame)
        dump_idx = start + frame * step
        if one_file:
            begin, lines = aa.renew_coordinate(trj_file_path, step=int(dump_idx), begin=begin, lines=lines,
                                               renew_v=renew_v)
        else:
            aa.renew_coordinate(trj_file_path, step=int(dump_idx), one_file=one_file)
        atom3d.renew_coordinate(aa)
        if renew_v:
            atom3d.renew_mv2(aa)
        if fake:
            yield atom3d.get_fake_copy()
        else:
            yield atom3d
        frame += 1


class Atom3dGenerator:
    def __init__(self, data_file_path, trj_file_path, one_file, info_path=None, renew_v=False, pad=0):
        self.data_file_path = data_file_path
        self.trj_file_path = trj_file_path
        self.one_file = one_file
        if info_path is None:
            self.is_aa = False
        else:
            self.is_aa = True
        if self.is_aa:
            self.atom3d = csc.Atom3D()
            self.atom3d.create_from_info(info_path)
            self.aa = csc.AA3D()
            self.aa.get_mass(data_file_path)
        else:
            self.parser = LmpParser()
            self.atom3d = self.parser.load_data_file(data_file_path, _type="fg")
        self.renew_v = renew_v
        self.pad = pad
        self.file_lines = None
        self.file_end = 0

    def __call__(self, time_step=-1, fake=True, **kwargs):
        if time_step < 0:
            return self.atom3d
        else:
            if self.is_aa:
                if self.one_file:
                    if self.file_lines is None:
                        self.file_end, self.file_lines = self.aa.renew_coordinate(self.trj_file_path,
                                                                                  step=int(time_step),
                                                                                  one_file=self.one_file,
                                                                                  renew_v=self.renew_v, pad=self.pad,
                                                                                  **kwargs)
                    else:
                        self.file_end, _ = self.aa.renew_coordinate(self.trj_file_path, step=int(time_step),
                                                                    one_file=self.one_file,
                                                                    begin=self.file_end, lines=self.file_lines,
                                                                    renew_v=self.renew_v, pad=self.pad,
                                                                    **kwargs)
                else:
                    self.aa.renew_coordinate(self.trj_file_path, step=int(time_step), one_file=self.one_file,
                                             renew_v=self.renew_v, pad=self.pad, **kwargs)

                self.atom3d.renew_coordinate(self.aa, **kwargs)
            else:
                self.atom3d = self.parser.renew_coordinate_file(self.atom3d, self.trj_file_path, int(time_step),
                                                                one_file=self.one_file, pad=self.pad, **kwargs)
            if fake:
                return self.atom3d.get_fake_copy(**kwargs)
            else:
                return self.atom3d

    def get_name(self, time_step=-1):
        if time_step < 0:
            return self.data_file_path
        if self.one_file:
            return self.trj_file_path
        else:
            return self.trj_file_path + str(int(time_step)).zfill(self.pad)


class Atom3dGeneratorV2:
    def __init__(self, data_file_path, one_file, info_path=None, renew_v=False, pad=0):
        self.data_file_path = data_file_path
        self.one_file = one_file
        if info_path is None:
            self.is_aa = False
        else:
            self.is_aa = True
        if self.is_aa:
            self.atom3d = csc.Atom3D()
            self.atom3d.create_from_info(info_path)
            self.aa = csc.AA3D()
            self.aa.get_mass(data_file_path)
        else:
            self.parser = LmpParser()
            self.atom3d = self.parser.load_data_file(data_file_path, _type="fg")
        self.renew_v = renew_v
        self.pad = pad
        self.file_lines = None
        self.file_end = 0

    def __call__(self, time_step=-1, trj_path="", fake=True, **kwargs):
        if time_step < 0:
            return self.atom3d
        else:
            assert trj_path
            if self.is_aa:
                if self.one_file:
                    if self.file_lines is None:
                        self.file_end, self.file_lines = self.aa.renew_coordinate(trj_path,
                                                                                  step=int(time_step),
                                                                                  one_file=self.one_file,
                                                                                  renew_v=self.renew_v, pad=self.pad,
                                                                                  **kwargs)
                    else:
                        self.file_end, _ = self.aa.renew_coordinate(trj_path, step=int(time_step),
                                                                    one_file=self.one_file,
                                                                    begin=self.file_end, lines=self.file_lines,
                                                                    renew_v=self.renew_v, pad=self.pad,
                                                                    **kwargs)
                else:
                    self.aa.renew_coordinate(trj_path, step=int(time_step), one_file=self.one_file,
                                             renew_v=self.renew_v, pad=self.pad, **kwargs)

                self.atom3d.renew_coordinate(self.aa, **kwargs)
            else:
                self.atom3d = self.parser.renew_coordinate_file(self.atom3d, trj_path, int(time_step),
                                                                one_file=self.one_file, pad=self.pad, **kwargs)
            if fake:
                return self.atom3d.get_fake_copy(**kwargs)
            else:
                return self.atom3d


class Atom3dGenCGU2CG(Atom3dGenerator):
    def __init__(self, data_file_path, trj_file_path, one_file, renew_v=False, pad=0):
        super().__init__(data_file_path, trj_file_path, one_file, info_path=None, renew_v=renew_v, pad=pad)
        self.cg = cp.deepcopy(self.atom3d)
        self.cg.cgu_to_cg()

    def __call__(self, time_step=-1, fake=True, **kwargs):
        atom3d = super().__call__(time_step=time_step, fake=False, **kwargs)
        self.cg.renew_from_cgu(atom3d, **kwargs)
        if time_step < 0:
            return self.cg
        if fake:
            return self.cg.get_fake_copy(**kwargs)
        else:
            return self.cg


class Atom3dGenCM(Atom3dGenerator):
    def __init__(self, data_file_path, trj_file_path, one_file, info_path=None, pad=0):
        super().__init__(data_file_path, trj_file_path, one_file, info_path, pad=pad)
        self.cm_atom3d = csc.Atom3D()
        # print(self.atom3d.mole_dict)
        # print(self.atom3d)

        for k in self.atom3d.mole_dict:
            a = bsc.Atom(k, _type="cm", mole_id=k)
            self.cm_atom3d.append_element(a)

    def __call__(self, time_step=-1, **kwargs):
        if time_step < 0:
            return self.cm_atom3d
        else:
            new_a3d = super().__call__(time_step)
            new_a3d.convert_to_default()
            for a_id in self.cm_atom3d.Atoms:
                co = np.zeros(3)
                m = 0
                for aa_id in self.atom3d.mole_dict[a_id]:
                    atom = new_a3d.Atoms[aa_id]
                    mi = self.atom3d.type_mass_dict[atom.type]
                    co += atom.coordinate * mi
                    m += mi
                self.cm_atom3d.Atoms[a_id].coordinate = co / m
            self.cm_atom3d.default = True
            self.cm_atom3d.in_cell = False
            return self.cm_atom3d


class Atom3dGenSegCM(Atom3dGenerator):
    def __init__(self, data_file_path, trj_file_path, one_file, info_path=None, pad=0):
        super().__init__(data_file_path, trj_file_path, one_file, info_path, pad=pad)
        self.cm_atom3d = csc.Atom3D()
        td = {True: "soft", False: "hard"}
        self.seg_dict = dict()
        seg_id = 0
        # print(self.atom3d)
        # TODO realize this for 1 blk only.
        for k, mole in self.atom3d.mole_dict.items():
            seg_flag = False
            for a_id in mole:
                atom = self.atom3d.Atoms[a_id]
                soft_bool = atom.type[:2] == "TO" or atom.type[:2] == "Es"
                if atom.type[0].isupper():
                    if (soft_bool and (not seg_flag)) or ((not soft_bool) and seg_flag) or seg_id == 0:
                        seg_id += 1
                        seg_flag = not seg_flag
                        self.seg_dict[seg_id] = [a_id]
                        a = bsc.Atom(seg_id, _type=td[seg_flag], mole_id=k)
                        self.cm_atom3d.append_element(a)
                    else:
                        self.seg_dict[seg_id].append(a_id)
                else:
                    self.seg_dict[seg_id - 1].append(a_id)
        print(self.seg_dict)

    def __call__(self, time_step=-1, **kwargs):
        if time_step < 0:
            return self.cm_atom3d
        else:
            new_a3d = super().__call__(time_step)
            new_a3d.convert_to_default()
            for a_id in self.cm_atom3d.Atoms:
                co = np.zeros(3)
                m = 0
                for aa_id in self.seg_dict[a_id]:
                    atom = new_a3d.Atoms[aa_id]
                    mi = self.atom3d.type_mass_dict[atom.type]
                    co += atom.coordinate * mi
                    m += mi
                self.cm_atom3d.Atoms[a_id].coordinate = co / m
            self.cm_atom3d.default = True
            self.cm_atom3d.in_cell = False
            return self.cm_atom3d


class Atom3dGenSegCMV2(Atom3dGenerator):
    def __init__(self, data_file_path, trj_file_path, one_file, info_path=None, pad=0):
        super().__init__(data_file_path, trj_file_path, one_file, info_path, pad=pad)
        self.cm_atom3d = csc.Atom3D()

        seg_dict, head_tail = self.atom3d.set_seg_id()
        self.seg_dict = dict()

        for k, atoms in seg_dict.items():
            if k not in head_tail:
                seg_type = "hard"
                for a_id in atoms:
                    if self.atom3d.Atoms[a_id].type == "TO(2)":
                        seg_type = "soft"
                        break
                a = bsc.Atom(k, _type=seg_type, mole_id=k)
                self.cm_atom3d.append_element(a)
                self.seg_dict[k] = atoms

    def __call__(self, time_step=-1, **kwargs):
        if time_step < 0:
            return self.cm_atom3d
        else:
            new_a3d = super().__call__(time_step)
            new_a3d.convert_to_default()
            for a_id in self.cm_atom3d.Atoms:
                co = np.zeros(3)
                m = 0
                for aa_id in self.seg_dict[a_id]:
                    atom = new_a3d.Atoms[aa_id]
                    mi = self.atom3d.type_mass_dict[atom.type]
                    co += atom.coordinate * mi
                    m += mi
                self.cm_atom3d.Atoms[a_id].coordinate = co / m
            self.cm_atom3d.default = True
            self.cm_atom3d.in_cell = False
            return self.cm_atom3d


class DynamicAnalyzer:
    # @staticmethod
    # def mean_squared_displacement(atom3ds, n_repeat):
    #     def get_type_co_dict(a3d):
    #         a3d.convert_to_default()
    #         co_dict = dict()
    #         for _at in a3d.Atom_type_dict.keys():
    #             co_dict[_at] = list()
    #         # order not ensured, may cause bug
    #         for a in a3d.Atoms.values():
    #             co_dict[a.type].append(a.coordinate)
    #         for _at in a3d.Atom_type_dict.keys():
    #             co_dict[_at] = np.vstack(co_dict[_at])
    #         return co_dict
    #
    #     def cal_msd_1(last_co_d, new_co_d):
    #         msd_dict = dict()
    #         for k in last_co_d.keys():
    #             n = new_co_d[k]
    #             l = last_co_d[k]
    #             d = n - l
    #             d2 = np.square(d)
    #             msd_dict[k] = np.mean(np.sum(d2, axis=1))
    #         return msd_dict
    #
    #     atom3d = next(atom3ds)
    #     last_co_dict = get_type_co_dict(atom3d)
    #     msd = dict()
    #     for at in atom3d.Atom_type_dict:
    #         msd[at] = 0.
    #     for i in range(1, n_repeat + 1):
    #         new_atom3d = next(atom3ds)
    #         new_co_dict = get_type_co_dict(new_atom3d)
    #         n_msd = cal_msd_1(last_co_dict, new_co_dict)
    #         for at in atom3d.Atom_type_dict:
    #             msd[at] += n_msd[at]
    #         last_co_dict = new_co_dict
    #     for at in atom3d.Atom_type_dict:
    #         msd[at] /= n_repeat
    #     return msd

    @staticmethod
    def mean_squared_displacement(atom3d_gt, start, steps):
        def get_type_co_dict(a3d):
            a3d.convert_to_default()
            co_dict = dict()
            for _at in a3d.Atom_type_dict.keys():
                co_dict[_at] = list()

            atom_ids = list(a3d.Atoms.keys())
            atom_ids.sort()
            for a_id in atom_ids:
                a = a3d.Atoms[a_id]
                co_dict[a.type].append(a.coordinate)
            for _at in a3d.Atom_type_dict.keys():
                co_dict[_at] = np.vstack(co_dict[_at])
            return co_dict

        def cal_msd_1(last_co_d, new_co_d):
            msd_dict = dict()
            for k in last_co_d.keys():
                n = new_co_d[k]
                l = last_co_d[k]
                d = n - l
                d2 = np.square(d)
                msd_dict[k] = np.mean(np.sum(d2, axis=1))
            return msd_dict

        print(start)
        msd = dict()
        for at in atom3d_gt().Atom_type_dict:
            msd[at] = list()
        msd["time"] = list(steps)
        ori_atom3d = cp.deepcopy(atom3d_gt(start))
        ori_co_dict = get_type_co_dict(ori_atom3d)
        # print(ori_co_dict)
        for i in steps:
            # last_co_dict = get_type_co_dict(atom3d_gt(i))
            new_co_dict = get_type_co_dict(atom3d_gt(i + start))
            n_msd = cal_msd_1(ori_co_dict, new_co_dict)
            for at in atom3d_gt().Atom_type_dict:
                msd[at].append(n_msd[at])
        return msd

    @staticmethod
    def diffusion_coefficient(atom3ds, step, dt=5):
        last_atom3d = next(atom3ds)
        time_steps = [0]
        msd_dict = dict()
        temp_dict = dict()
        diffusion_coefficient_dict = dict()
        for _type in last_atom3d.Atom_type_dict:
            temp_dict[_type] = 0.
            msd_dict[_type] = [0.]
        for new_atom3d in atom3ds:
            time_steps.append(time_steps[-1] + dt * step)
            for atom_id, atom in last_atom3d.Atoms.items():
                co1 = atom.coordinate
                co2 = new_atom3d.Atoms[atom_id].get_image_coordinate(atom.image, new_atom3d.lattice_parameter)
                d = co1 - co2
                # if (d > 20.).any():
                #     print(d)
                temp_dict[atom.type] += np.inner(d, d)
            for _type in last_atom3d.Atom_type_dict:
                msd_dict[_type].append(temp_dict[_type] / len(last_atom3d.Atom_type_dict[_type]))
                temp_dict[_type] = 0.
            # last_atom3d = new_atom3d
        for _type in last_atom3d.Atom_type_dict:
            popt, _ = curve_fit(func.linear_func, time_steps[1:], msd_dict[_type][1:], p0=[1., 0.])
            # print(msd_dict[_type][:10])
            # if _type in ["U(0)", "TO(1)"]:
            plt.plot(time_steps[1:], msd_dict[_type][1:])
            plt.title(_type)
            plt.show()
            diffusion_coefficient_dict[_type] = popt[0] / 6
        return diffusion_coefficient_dict


class StaticAnalyzer:
    hard_type = {"Ph", "U(", "B(", "c1", "c2", "c3", "n1", "o1"}

    # @staticmethod
    # def grid_hard_domain(atom3d: csc.Atom3D, grid=0.5, r=3.):
    #     def ball_index():
    #         in_index_list = list()
    #         out_index_list = list()
    #         num = int(r / grid)
    #         n2 = num * num
    #         n2_out = (2 * num + 1) * (2 * num + 1)
    #         for i in range(-2 * num - 1, 2 * num + 2):
    #             for j in range(-2 * num - 1, 2 * num + 2):
    #                 for k in range(-2 * num - 1, 2 * num + 2):
    #                     d2 = np.inner((i, j, k), (i, j, k))
    #                     if d2 <= n2:
    #                         in_index_list.append((i, j, k))
    #                     elif d2 <= n2_out:
    #                         out_index_list.append((i, j, k))
    #         return in_index_list, out_index_list
    #
    #     def grid_atoms():
    #         hard_id_set = set()
    #         grid_atom_dict = dict()
    #         atom_grid_dict = dict()
    #         for atom_id, atom in atom3d.Atoms.items():
    #             if atom.type[:2] in ThermoAnalyzer.hard_type:
    #                 hard_id_set.add(atom_id)
    #                 atom_grid = tuple(int(co / grid) for co in (atom.coordinate - atom3d.box_l))
    #                 if atom_grid in grid_atom_dict:
    #                     print("more than one atom in one grid")
    #                 grid_atom_dict[atom_grid] = atom_id
    #
    #                 atom_grid_dict[atom_id] = atom_grid
    #         return hard_id_set, grid_atom_dict, atom_grid_dict
    #
    #     atom3d.convert_to_in_cell()
    #     in_ball, out_ball = ball_index()
    #     box_length = (atom3d.lattice_parameter / grid).astype(int)
    #     record_mat = np.zeros(box_length)
    #     volume_dict = dict()
    #     a_set, g_a_dict, a_g_dict = grid_atoms()
    #     print(a_set)
    #     # print(g_a_dict)(75, 76, 69)
    #     # print(a_g_dict)
    #     expand_set = set()
    #     current_id = 0
    #     while a_set:
    #         print("a_set,", a_set)
    #         print(len(a_set))
    #         a_id = a_set.pop()
    #         g_a_dict.pop(a_g_dict[a_id])
    #         expand_set.add(a_id)
    #         current_id += 1
    #         volume_dict[current_id] = 0
    #         while expand_set:
    #             # print("e_set,", expand_set)
    #             exp_id = expand_set.pop()
    #             cent = a_g_dict.pop(exp_id)
    #             for plus in in_ball:
    #                 co = (cent[0] + plus[0], cent[1] + plus[1], cent[2] + plus[2])
    #                 n_co = np.asarray(co)
    #                 if (0 <= n_co).all() and (n_co < box_length).all():
    #                     if record_mat[co] == 0:
    #                         record_mat[co] = current_id
    #                         volume_dict[current_id] += 1
    #                         if co in g_a_dict:
    #                             ne_id = g_a_dict.pop(co)
    #                             expand_set.add(ne_id)
    #                             a_set.remove(ne_id)
    #                     elif record_mat[co] != current_id:
    #                         print(co)
    #                         print(record_mat[co])
    #                         raise ValueError("bug1")
    #             for plus in out_ball:
    #                 co = (cent[0] + plus[0], cent[1] + plus[1], cent[2] + plus[2])
    #                 n_co = np.asarray(co)
    #                 if (0 <= n_co).all() and (n_co < box_length).all():
    #                     if co in g_a_dict:
    #                         ne_id = g_a_dict.pop(co)
    #                         expand_set.add(ne_id)
    #                         a_set.remove(ne_id)
    #                         if record_mat[co] != current_id and record_mat[co] != 0:
    #                             print(co)
    #                             print(record_mat[co])
    #                             raise ValueError("bug2")
    #     return volume_dict

    @staticmethod
    def grid_hard_domain(atom3d: csc.Atom3D, grid=0.5, r=3., out_path=None):
        class Node:
            def __init__(self, node_id, father=None):
                self.id = node_id
                self.father = father

            def get_grand(self):
                if self.father is None:
                    return self
                g = self
                f = self.father
                while f is not None:
                    g = f
                    f = f.father
                self.father = g
                return g

        # def ball_index():
        #     in_index_list = list()
        #     num = int(r / grid) + 1
        #     n2 = num * num
        #     for i in range(-num, num + 1):
        #         for j in range(-num, num + 1):
        #             for k in range(-num, num + 1):
        #                 d2 = np.inner((i, j, k), (i, j, k))
        #                 if d2 <= n2:
        #                     in_index_list.append((i, j, k))
        #     return in_index_list
        def ball_index():
            in_index_list = list()
            out_index_list = list()
            num = int(r / grid)
            n2 = (num - 1) * (num - 1)
            n2_out = num * num
            for i in range(-num, num + 1):
                for j in range(-num, num + 1):
                    for k in range(-num, num + 1):
                        d2 = np.inner((i, j, k), (i, j, k))
                        if d2 <= n2:
                            in_index_list.append((i, j, k))
                        elif d2 <= n2_out:
                            out_index_list.append((i, j, k))
            return in_index_list, out_index_list

        def grid_atoms():
            atom_grid_dict = dict()
            for atom_id, atom in atom3d.Atoms.items():
                if atom.type[:2] in StaticAnalyzer.hard_type:
                    atom_grid = tuple(int(co / grid) for co in (atom.coordinate - atom3d.box_l))
                    atom_grid_dict[atom_id] = atom_grid
            return atom_grid_dict

        def get_in_cell_co(_co):
            i, j, k = _co
            if i < 0:
                i += box_length[0]
            elif i >= box_length[0]:
                i -= box_length[0]
            if j < 0:
                j += box_length[1]
            elif j >= box_length[1]:
                j -= box_length[1]
            if k < 0:
                k += box_length[2]
            elif k >= box_length[2]:
                k -= box_length[2]
            return (i, j, k)

        atom3d.convert_to_in_cell()
        in_ball, out_ball = ball_index()
        box_length = (atom3d.lattice_parameter / grid).astype(int)
        record_mat = np.zeros(box_length, dtype=int)
        draw_mat = np.zeros(box_length, dtype=int)
        a_g_dict = grid_atoms()
        node_dict = dict()
        hard_id = 0
        for a_id, cent in a_g_dict.items():
            hard_id += 1
            if hard_id % 1000 == 0:
                print("hard,", hard_id)
            node_dict[hard_id] = Node(hard_id)
            for plus in in_ball:
                co = (cent[0] + plus[0], cent[1] + plus[1], cent[2] + plus[2])
                co = get_in_cell_co(co)
                if record_mat[co] == 0:
                    record_mat[co] = hard_id
                else:
                    n1 = node_dict[record_mat[co]].get_grand()
                    n2 = node_dict[hard_id].get_grand()
                    if n1.id != n2.id:
                        n2.father = n1
            for plus in out_ball:
                co = (cent[0] + plus[0], cent[1] + plus[1], cent[2] + plus[2])
                co = get_in_cell_co(co)
                if record_mat[co] == 0:
                    record_mat[co] = hard_id
                    draw_mat[co] = 1
                else:
                    n1 = node_dict[record_mat[co]].get_grand()
                    n2 = node_dict[hard_id].get_grand()
                    if n1.id != n2.id:
                        n2.father = n1

        hard_id_id = dict()
        for i in range(1, hard_id + 1):
            f_id = node_dict[i].get_grand().id
            hard_id_id[i] = f_id

        file_lines = list()
        f_id_type_id = dict()
        volume_dict = dict()
        type_id = 1
        atom_id = 1
        for i in range(box_length[0]):
            for j in range(box_length[1]):
                for k in range(box_length[2]):
                    h_id = record_mat[i, j, k]
                    if h_id != 0:
                        f_id = hard_id_id[h_id]
                        if f_id in f_id_type_id:
                            t_id = f_id_type_id[f_id]
                            volume_dict[t_id] += 1
                        else:
                            f_id_type_id[f_id] = type_id
                            volume_dict[type_id] = 1
                            t_id = type_id
                            type_id += 1
                        if draw_mat[i, j, k]:
                            file_lines.append("%d\t1\t%d\t%d\t%d\t%d\n" % (atom_id, t_id, i, j, k))
                        atom_id += 1
        pre_lines = ["LAMMPS data file\n\n", "\t%d atoms\n" % len(file_lines), "\t%d atom types\n" % len(volume_dict),
                     "\t%.2f  %.2f xlo xhi\n" % (-grid / 2, box_length[0] - 1 + grid / 2),
                     "\t%.2f  %.2f ylo yhi\n" % (-grid / 2, box_length[1] - 1 + grid / 2),
                     "\t%.2f  %.2f ylo yhi\n" % (-grid / 2, box_length[2] - 1 + grid / 2),
                     "\nAtoms\n\n"]
        if out_path:
            with open(out_path, "w") as f:
                f.writelines(pre_lines + file_lines)
        return volume_dict

    @staticmethod
    def grid_bead_2d(co_list, grid=0.5, bin_num=100, begin=np.zeros(3), drop_axis=0):
        begin = np.delete(begin, drop_axis, axis=0)
        xx = np.arange(0., bin_num + 1) * grid + begin[0]  # + box_l_da[0]
        yy = np.arange(0., bin_num + 1) * grid + begin[1]  # + box_l_da[1]
        co_list = np.delete(co_list, drop_axis, axis=1)
        H, bin_edges = np.histogramdd(co_list, bins=(xx, yy))
        return H

    @staticmethod
    def gaussian_kernel_convolve_2d(data, kernel_size, s2):
        def gauss():
            kernel = np.zeros((kernel_size, kernel_size))
            center = kernel_size // 2
            sum_val = 0
            for i in range(kernel_size):
                for j in range(kernel_size):
                    x, y = i - center, j - center
                    kernel[i, j] = np.exp(-(x ** 2 + y ** 2) / 2 * s2)
                    sum_val += kernel[i, j]
            kernel = kernel / sum_val
            return kernel

        kern = gauss()
        ret = convolve2d(data, kern, "same")
        return ret

    @staticmethod
    def beads_density(co_list, bead_size, grid=0.5, bin_num=100, begin=np.zeros(2), drop_axis=0):
        kernel_size = int(bead_size / 2 / grid)
        s2 = 0.01
        H = StaticAnalyzer.grid_bead_2d(co_list, grid, bin_num, begin, drop_axis)
        # return H
        density = StaticAnalyzer.gaussian_kernel_convolve_2d(H, kernel_size, s2)
        return density

    @staticmethod
    def saxs(atom3d: csc.Atom3D, grid=0.5, bin_num=120, iter_num=1, drop_axis=0):
        bead_size = {"TO(1)": 4.98, "TO(2)": 5.07, "Es": 4.65, "Me": 5.1, "Ph": 4.99, "U": 4.49,
                     "n1": 3.7, "h1": 1.45, "c1": 3.9, "o1": 3.43}
        electron = {"TO(1)": 41, "TO(2)": 40, "Es": 38, "Me": 8, "Ph": 40, "U": 30, "n1": 7, "h1": 1, "c1": 6, "o1": 6}
        assert (atom3d.lattice_parameter > grid * bin_num).all()
        atom3d.convert_to_in_cell()
        co_dict = {k: list() for k in atom3d.Atom_type_dict.keys()}
        for atom in atom3d.Atoms.values():
            co_dict[atom.type].append(atom.coordinate.copy())
        for k in atom3d.Atom_type_dict.keys():
            co_dict[k] = np.vstack(co_dict[k])
        fd_aver = np.zeros((bin_num, bin_num))
        np.random.seed(0)
        for i in range(iter_num):
            begin = (atom3d.lattice_parameter - grid * bin_num) * np.random.rand(1) + atom3d.box_l
            density_sum = np.zeros((bin_num, bin_num))
            for at, co_list in co_dict.items():
                if at == "TO(2)" or at == "TO(1)":
                    continue
                density_sum += StaticAnalyzer.beads_density(co_list, bead_size[at], grid, bin_num, begin, drop_axis) * \
                               electron[at]
                # atom3d.type_mass_dict[at]
            density_sum -= density_sum.mean()
            # fd_aver += np.abs(np.fft.fftshift(np.fft.fft2(density_sum))) ** 2
            fd_aver += np.abs(np.fft.fft2(density_sum)) ** 2
        return fd_aver / iter_num


class HydrogenBond:
    def __init__(self, atom3d, h_len=3.0, h_ang=120.):
        self.h_len = h_len
        self.h_ang = h_ang
        self.donors, self.acceptors, self.d1_c1, self.a_c1 = self.create_donors_acceptors(atom3d)

    @staticmethod
    def create_donors_acceptors(atom3d):
        _donors = [[], []]
        _acceptors = []
        d1_c1 = dict()
        a_c1 = dict()
        for atom_id, atom in atom3d.Atoms.items():
            if atom.type == "o1":
                _acceptors.append(atom_id)
                c1 = atom3d.get_connected_atom_id(atom_id).pop()
                a_c1[atom_id] = c1
            elif atom.type == "h1":
                h1 = atom_id
                n1 = atom3d.get_connected_atom_id(h1).pop()
                _donors[0].append(h1)
                _donors[1].append(n1)
                for con_id in atom3d.get_connected_atom_id(n1):
                    if atom3d.Atoms[con_id].type == "c1":
                        d1_c1[len(_donors[0]) - 1] = con_id
        return _donors, _acceptors, d1_c1, a_c1

    def hbond_judge(self, o1_co, h1_co, n1_co):
        hb_len = cso.calculate_distance(o1_co, h1_co)
        if hb_len > self.h_len:
            return False
        else:
            hb_ang = cso.calculate_angle(o1_co, h1_co, n1_co)
            if hb_ang < self.h_ang:
                return False
            else:
                return True

    @staticmethod
    def create_extended_acceptors(atom3d, _acceptors, extend_range=6.):
        def reflect_acceptor(ac, x, y, z):
            id_list.append(ac.id)
            d = np.asarray([x, y, z], dtype=int)
            co_list.append(ac.coordinate - d * atom3d.lattice_parameter)
            image_list.append(ac.image + d)

        id_list = list()
        co_list = list()
        image_list = list()
        # used to be a bug
        # upper = atom3d.lattice_parameter + atom3d.box_h - extend_range
        assert all((atom3d.box_h - atom3d.box_l - atom3d.lattice_parameter) == 0)
        upper = atom3d.box_h - extend_range
        lower = atom3d.box_l + extend_range
        for o1_id in _acceptors:
            o1 = atom3d.Atoms[o1_id]
            id_list.append(o1_id)
            co_list.append(o1.coordinate)
            image_list.append(o1.image)
            record = (o1.coordinate > upper).astype(int) + (o1.coordinate <= lower).astype(int) * -1
            if record[0] != 0:
                reflect_acceptor(o1, record[0], 0, 0)
                if record[1] != 0:
                    reflect_acceptor(o1, record[0], record[1], 0)
                    if record[2] != 0:
                        reflect_acceptor(o1, record[0], record[1], record[2])
                if record[2] != 0:
                    reflect_acceptor(o1, record[0], 0, record[2])
            if record[1] != 0:
                reflect_acceptor(o1, 0, record[1], 0)
                if record[2] != 0:
                    reflect_acceptor(o1, 0, record[1], record[2])
            if record[2] != 0:
                reflect_acceptor(o1, 0, 0, record[2])
        return id_list, co_list, image_list

    def cal_hbond(self, atom3d):
        # h bond saved in a dict of dict: {donor_id: {acceptor_atom_id: image_dif}}
        h_bond_dict = dict()
        atom3d.convert_to_in_cell(only_aa=True)
        id_list, co_list, image_list = self.create_extended_acceptors(atom3d, self.acceptors,
                                                                      extend_range=self.h_len + 2.)
        ac_tree = cKDTree(co_list)
        cut_off = self.h_len + 0.3
        for d_id, h1_id in enumerate(self.donors[0]):
            h1 = atom3d.Atoms[h1_id]
            h1_co = h1.coordinate
            n1_id = self.donors[1][d_id]
            n1_co = atom3d.Atoms[n1_id].get_image_coordinate(h1.image, atom3d.lattice_parameter)
            idx_list = ac_tree.query_ball_point(h1_co, cut_off)
            for idx in idx_list:
                if self.hbond_judge(co_list[idx], h1_co, n1_co):
                    if d_id not in h_bond_dict:
                        h_bond_dict[d_id] = dict()
                    # else:
                    #     raise KeyError("Multiple h bond formed!")
                    if id_list[idx] in h_bond_dict[d_id]:
                        raise KeyError("acceptor already exists!")
                    else:
                        h_bond_dict[d_id][id_list[idx]] = tuple(h1.image - image_list[idx])
        return h_bond_dict

    def distribution(self, atom3d, cut_off=5.0):
        def configure(o1_co, h1_co, n1_co):
            hb_len = cso.calculate_distance(o1_co, h1_co)
            hb_ang = cso.calculate_angle(o1_co, h1_co, n1_co)
            return hb_len, hb_ang

        hbonds = list()
        atom3d.convert_to_in_cell(only_aa=True)
        id_list, co_list, image_list = self.create_extended_acceptors(atom3d, self.acceptors,
                                                                      extend_range=cut_off + 1.)
        ac_tree = cKDTree(co_list)
        for d_id, h1_id in enumerate(self.donors[0]):
            h1 = atom3d.Atoms[h1_id]
            h1_co = h1.coordinate
            n1_id = self.donors[1][d_id]
            n1_co = atom3d.Atoms[n1_id].get_image_coordinate(h1.image, atom3d.lattice_parameter)
            idx_list = ac_tree.query_ball_point(h1_co, cut_off + 0.3)
            for idx in idx_list:
                hbonds.append(configure(co_list[idx], h1_co, n1_co))
        return hbonds

    def h_vec(self, atom3d, cut_off=5.0, judge=False):

        def configure(o1_co, h1_co):
            return o1_co - h1_co

        hbonds = list()
        atom3d.convert_to_in_cell(only_aa=True)
        id_list, co_list, image_list = self.create_extended_acceptors(atom3d, self.acceptors,
                                                                      extend_range=cut_off + 1.)
        ac_tree = cKDTree(co_list)
        for d_id, h1_id in enumerate(self.donors[0]):
            h1 = atom3d.Atoms[h1_id]
            h1_co = h1.coordinate
            n1_id = self.donors[1][d_id]
            n1_co = atom3d.Atoms[n1_id].get_image_coordinate(h1.image, atom3d.lattice_parameter)
            idx_list = ac_tree.query_ball_point(h1_co, cut_off + 0.3)
            for idx in idx_list:
                o1_co = co_list[idx]
                if judge:
                    flag = self.hbond_judge(o1_co, h1_co, n1_co)
                    if not flag:
                        continue
                hbonds.append(configure(o1_co, h1_co))
        return hbonds

    def hb1_hb2(self, h_dict):
        record_dict = dict()
        for d_id in h_dict:
            for a_id in h_dict[d_id]:
                c1, c2 = self.d1_c1[d_id], self.a_c1[a_id]
                h_tuple = tuple(sorted([c1, c2]))
                if h_tuple in record_dict:
                    record_dict[h_tuple] += 1
                else:
                    record_dict[h_tuple] = 1
        h1, h2, h3 = 0, 0, 0
        for k, v in record_dict.items():
            if v == 1:
                h1 += 1
            elif v == 2:
                h2 += 1
            else:
                h2 += 1
                h3 += 1
                print("Might be error,", k, v)
        return h1, h2, h3


class HBAutoCorrelation:

    def __init__(self, _type, step, frames, atom_gen, h_len=3., h_ang=130.):
        self.type = _type
        self.step = step
        self.frames = frames
        self.atom_gen = atom_gen
        # self.h_len = h_len
        # self.h_ang = h_ang
        self.h_record = dict()
        # self.donors = None
        # self.acceptors = None
        self.h_cal = HydrogenBond(atom_gen(), h_len, h_ang)

    def create_h_record_single(self, frames, screen):
        record = dict()
        tl = len(frames) // 20
        for i, frame in enumerate(frames):
            if screen:
                if not i % tl:
                    print("%d%%" % int(i // tl * 5))
            record[frame] = self.h_cal.cal_hbond(self.atom_gen(frame, only_aa=True))
        return record

    def create_h_record(self, frames, parallel=24):
        frames = list(frames)
        seg = len(frames) // parallel
        pool = Pool(parallel)
        ret = list()
        for i in range(parallel):
            if i < parallel - 1:
                ret.append(pool.apply_async(self.create_h_record_single, args=(frames[i * seg: (i + 1) * seg], False)))
            else:
                ret.append(pool.apply_async(self.create_h_record_single, args=(frames[i * seg:], True)))
        pool.close()
        pool.join()
        h_record = dict()
        for r in ret:
            record = r.get()
            h_record.update(record)
        return h_record

    def get_h(self, time_step):
        if time_step not in self.h_record:
            self.h_record[time_step] = self.h_cal.cal_hbond(self.atom_gen(time_step, only_aa=True))
        return self.h_record[time_step]

    def cal_single(self, start):
        def cross_hbond(l_hbond, n_hbond, drop=False):
            h_num = 0
            for donor in list(l_hbond.keys()):
                if donor in n_hbond:
                    flag = 0
                    for acceptor in list(l_hbond[donor].keys()):
                        if acceptor in n_hbond[donor]:
                            if l_hbond[donor][acceptor] == n_hbond[donor][acceptor]:
                                h_num += 1
                                flag = 1
                        if (not flag) and drop:
                            l_hbond[donor].pop(acceptor)
                            if not l_hbond[donor]:
                                l_hbond.pop(donor)
                elif drop:
                    l_hbond.pop(donor)
            return h_num

        init_h = self.get_h(start)
        # print(init_h)
        init_sum = 0
        for v in init_h.values():
            init_sum += len(v)
        print(init_sum)
        if self.type == "intermittent":
            result = np.zeros(self.frames)
            result[0] = 1.
            for i in range(1, self.frames):
                # print(i)
                new_h = self.get_h(self.step * i + start)
                # print(len(new_h))
                h_num = cross_hbond(init_h, new_h, drop=False)
                result[i] = h_num / init_sum
            return result

        elif self.type == "continuous":
            orig_h = cp.deepcopy(init_h)
            result = np.zeros(self.frames)
            result[0] = 1.
            for i in range(1, self.frames):
                # print(orig_h)
                if not orig_h:
                    break
                new_h = self.get_h(self.step * i + start)
                h_num = cross_hbond(orig_h, new_h, drop=True)
                result[i] = h_num / init_sum
            return result

    def cal_multiple(self, starts):
        result = np.zeros(self.frames)
        for start in starts:
            print(start)
            result += self.cal_single(start)
        result /= len(starts)
        return result

    def cal_multiple_single(self, starts):
        result = np.zeros(self.frames)
        for start in starts:
            print(start)
            result += self.cal_single(start)
        return result

    def cal_multiple_parallel(self, starts, parallel=24):
        pool = Pool(parallel)
        ret = list()
        seg = len(starts) // parallel
        for i in range(parallel):
            if i < parallel - 1:
                ret.append(pool.apply_async(self.cal_multiple_single, args=(starts[i * seg:(i + 1) * seg],)))
            else:
                ret.append(pool.apply_async(self.cal_multiple_single, args=(starts[i * seg:],)))
        pool.close()
        pool.join()
        result = np.zeros(self.frames)
        for r in ret:
            result += r.get()
        return result / len(starts)

    def life_time(self, start):
        def increase_time():
            for k in alive_dict:
                alive_dict[k] += 1

        def form_break(lh, nh):
            # break
            for donor in lh:
                if donor in nh:
                    for acceptor in lh[donor]:
                        if acceptor not in nh[donor]:
                            if (donor, acceptor) in alive_dict:
                                times.append(alive_dict.pop((donor, acceptor)) - 1)
                else:
                    for acceptor in lh[donor]:
                        if (donor, acceptor) in alive_dict:
                            times.append(alive_dict.pop((donor, acceptor)) - 1)
            # form
            for donor in nh:
                if donor in lh:
                    for acceptor in nh[donor]:
                        if acceptor not in lh[donor]:
                            alive_dict[(donor, acceptor)] = 1
                else:
                    for acceptor in nh[donor]:
                        alive_dict[(donor, acceptor)] = 1

        alive_dict = dict()
        times = list()
        last_h = self.get_h(start)
        for i in range(1, self.frames):
            if i % 100 == 0:
                print(i)
            new_h = self.get_h(self.step * i + start)
            increase_time()
            form_break(last_h, new_h)
            last_h = new_h
            # print(alive_dict)
        return times

    def break_form(self, start):
        def form_break(lh, nh):
            n_break = 0
            for donor in lh:
                if donor in nh:
                    for acceptor in lh[donor]:
                        if acceptor not in nh[donor]:
                            n_break += 1
                else:
                    n_break += len(lh[donor])
            # form
            n_form = 0
            for donor in nh:
                if donor in lh:
                    for acceptor in nh[donor]:
                        if acceptor not in lh[donor]:
                            n_form += 1
                else:
                    n_form += len(nh[donor])
            return n_break, n_form

        last_h = self.get_h(start)
        breaks = list()
        forms = list()
        for i in range(1, self.frames):
            new_h = self.get_h(self.step * i + start)
            n_break, n_form = form_break(last_h, new_h)
            breaks.append(n_break)
            forms.append(n_form)
            last_h = new_h
        return breaks, forms

    def life_time_da(self, start):
        def increase_time():
            for k in alive_dict:
                alive_dict[k] += 1

        def form_break(lh, nh):
            # break
            for donor in lh:
                if donor in nh:
                    for acceptor in lh[donor]:
                        if acceptor not in nh[donor]:
                            if (donor, acceptor) in alive_dict:
                                donor_id = self.h_cal.donors[0][donor]
                                if (donor_id, acceptor) not in da_time_dict:
                                    da_time_dict[(donor_id, acceptor)] = [alive_dict.pop((donor, acceptor)) - 1]
                                else:
                                    da_time_dict[(donor_id, acceptor)].append(alive_dict.pop((donor, acceptor)) - 1)
                else:
                    for acceptor in lh[donor]:
                        if (donor, acceptor) in alive_dict:
                            donor_id = self.h_cal.donors[0][donor]
                            if (donor_id, acceptor) not in da_time_dict:
                                da_time_dict[(donor_id, acceptor)] = [alive_dict.pop((donor, acceptor)) - 1]
                            else:
                                da_time_dict[(donor_id, acceptor)].append(alive_dict.pop((donor, acceptor)) - 1)
            # form
            for donor in nh:
                if donor in lh:
                    for acceptor in nh[donor]:
                        if acceptor not in lh[donor]:
                            alive_dict[(donor, acceptor)] = 1
                else:
                    for acceptor in nh[donor]:
                        alive_dict[(donor, acceptor)] = 1

        alive_dict = dict()
        da_time_dict = dict()
        last_h = self.get_h(start)
        for i in range(1, self.frames):
            if i % 100 == 0:
                print(i)
            new_h = self.get_h(self.step * i + start)
            increase_time()
            form_break(last_h, new_h)
            last_h = new_h
            # print(alive_dict)
        return da_time_dict


class HBAutoCorrelationV2:
    def __init__(self, atom3d_gen:Atom3dGeneratorV2, h_len=3., h_ang=130.):
        # only used for intermittent type
        self.h_record = dict()
        self.h_cal = HydrogenBond(atom3d_gen(), h_len, h_ang)
        self.atom3d_gen = atom3d_gen

    def create_h_record_single(self, frames):
        record = dict()
        # tl = len(frames) // 20
        for frame in frames:
            # if screen:
            #     if not i % tl:
            #         print("%d%%" % int(i // tl * 5))
            timestep, trj_path = frame
            print(timestep)
            if timestep in self.h_record:
                continue
            record[timestep] = self.h_cal.cal_hbond(self.atom3d_gen(timestep, trj_path, only_aa=True))
        return record

    def create_h_record(self, frames, parallel=24):
        frames = list(frames)
        if parallel <= 0:
            h_record = self.create_h_record_single(frames)
            return h_record
        else:
            seg = len(frames) // parallel
            pool = Pool(parallel)
            ret = list()
            for i in range(parallel):
                if i < parallel - 1:
                    ret.append(pool.apply_async(self.create_h_record_single, args=(frames[i * seg: (i + 1) * seg],)))
                else:
                    ret.append(pool.apply_async(self.create_h_record_single, args=(frames[i * seg:], )))
            pool.close()
            pool.join()
            h_record = dict()
            for r in ret:
                record = r.get()
                h_record.update(record)
            return h_record

    def get_h(self, time_step):
        return self.h_record[time_step]

    def cal_single(self, t0_t1_pairs):
        def cross_hbond(l_hbond:Dict, n_hbond:Dict, drop=False):
            h_num = 0
            for donor in list(l_hbond.keys()):
                if donor in n_hbond:
                    flag = 0
                    for acceptor in list(l_hbond[donor].keys()):
                        if acceptor in n_hbond[donor]:
                            if l_hbond[donor][acceptor] == n_hbond[donor][acceptor]:
                                h_num += 1
                                flag = 1
                        if (not flag) and drop:
                            l_hbond[donor].pop(acceptor)
                            if not l_hbond[donor]:
                                l_hbond.pop(donor)
                elif drop:
                    l_hbond.pop(donor)
            sum_num = [len(_) for _ in n_hbond.values()]
            n_h_num = np.sum(sum_num) - h_num
            return h_num, n_h_num

        h_nums = list()
        for t0, t1 in t0_t1_pairs:
            print(t0, t1)
            init_h = self.get_h(t0)
            new_h = self.get_h(t1)
            h_num = cross_hbond(init_h, new_h, drop=False)
            h_nums.append(h_num)
        return h_nums

    def cal_multiple_parallel(self, t0_t1_pairs, parallel=24):
        pool = Pool(parallel)
        ret = list()
        seg = len(t0_t1_pairs) // parallel
        for i in range(parallel):
            if i < parallel - 1:
                ret.append(pool.apply_async(self.cal_single, args=(t0_t1_pairs[i * seg:(i + 1) * seg],)))
            else:
                ret.append(pool.apply_async(self.cal_single, args=(t0_t1_pairs[i * seg:],)))
        pool.close()
        pool.join()
        result = []
        for r in ret:
            result += r.get()
        return result


class PropertyCalculator:
    def pre_calculation(self, frame_set):
        pass

    def pre_calculation_parallel(self, frame_set, parallel=24):
        pass

    def get_correlation(self, frame1, frame2):
        return 0

    def pre_calculation_shared(self, val):
        return 0

    def get_correlation_shared(self, frame1, frame2, val):
        return 0


class URecord:
    def __init__(self, _type):
        self.u_record = dict()
        self.type = _type

    def create(self, atom3d, frame, cut_off):
        if frame not in self.u_record:
            if self.type == "cg":
                atom_type = "U"
            else:
                atom_type = "c1"
            central, extended = atom3d.get_central_and_extended_in_range(atom_type, cut_off)
            all_list = central + extended
            all_co = [atom.coordinate for atom in all_list]
            # print(all_co)
            tree = cKDTree(all_co)
            uu = dict()
            for u1 in central:
                idx_list = tree.query_ball_point(u1.coordinate, cut_off)
                if len(idx_list) > 1:
                    for idx in idx_list:
                        u2 = all_list[idx]
                        if u1.id == u2.id:
                            continue
                        if u2.id > u1.id:
                            uu[(u1.id, u2.id)] = u2.image - u1.image
                        else:
                            uu[(u2.id, u1.id)] = u1.image - u2.image
            self.u_record[frame] = uu

    def get(self, frame):
        if frame not in self.u_record:
            raise KeyError("not such frame")
        else:
            return self.u_record[frame]


class UreaCalculator(PropertyCalculator):
    def __init__(self, atom3d_gen: Atom3dGenerator, _type, cutoff, urecord=None):
        self.atom3d_gen = atom3d_gen
        self.cutoff = cutoff
        self.type = _type
        if urecord is None:
            self.urecord = URecord(_type)
        else:
            self.urecord = urecord
        self.core = dict()

    def pre_calculation(self, frame_set):
        for frame in frame_set:
            atom3d = self.atom3d_gen(frame, only_u=True)
            atom3d.convert_to_in_cell()
            self.urecord.create(atom3d, frame, self.cutoff)
        return self.urecord.u_record

    def pre_calculation_parallel(self, frame_set, parallel=24):
        starts = list(frame_set)
        if parallel < 1:
            result = self.pre_calculation(starts)
        else:
            pool = Pool(parallel)
            ret = list()
            seg = len(starts) // parallel
            for i in range(parallel):
                if i < parallel - 1:
                    ret.append(pool.apply_async(self.pre_calculation, args=(starts[i * seg:(i + 1) * seg],)))
                else:
                    ret.append(pool.apply_async(self.pre_calculation, args=(starts[i * seg:],)))
            pool.close()
            pool.join()
            result = dict()
            for r in ret:
                result.update(r.get())
        self.urecord.u_record = result
        # print(result)

    def get_correlation(self, frame1, frame2):
        if (frame1, frame2) not in self.core:
            uu1 = self.urecord.get(frame1)
            uu2 = self.urecord.get(frame2)
            up_num = 0
            for up in uu1:
                if up in uu2:
                    if (uu1[up] == uu2[up]).all():
                        up_num += 1
                        # print(up)
                        # print(uu1[up], uu2[up])
            self.core[(frame1, frame2)] = up_num / len(uu1)
        return self.core[(frame1, frame2)]


class StressCalculator(PropertyCalculator):
    def __init__(self, val, step):
        self.val = val
        self.step = step

    def get_correlation(self, frame1, frame2):
        # print(self.val[frame1//self.step])
        # print(self.val[frame2//self.step])
        return self.val[frame1 // self.step] * self.val[frame2 // self.step]

    def pre_calculation_shared(self, val):
        manager = Manager()
        shared_val = manager.Array('f', val)

        # shared_val = Array('f', val)
        return shared_val

    def get_correlation_shared(self, frame1, frame2, val):
        # print(self.val[frame1//self.step])
        # print(self.val[frame2//self.step])
        return val[frame1 // self.step] * val[frame2 // self.step]


class AutoCorrelation:
    def __init__(self, calculator: PropertyCalculator, frames, steps, starts):
        self.calculator = calculator
        self.frames = frames
        self.steps = steps
        self.starts = starts
        self.frame_set = set()

    def cal(self):
        if not self.frame_set:
            for start in self.starts:
                self.frame_set.update(set(np.arange(self.frames) * self.steps + start))
            pass
        # self.calculator.pre_calculation(self.frame_set)
        self.calculator.pre_calculation_parallel(self.frame_set, 20)
        print("finish pre calculation")
        result = np.zeros(self.frames)
        for start in self.starts:
            print(start)
            for frame in range(self.frames):
                if frame * self.steps < 10000000:
                    continue
                print(start, start + frame * self.steps)
                result[frame] += self.calculator.get_correlation(start, start + frame * self.steps)
        return result / len(self.starts)

    def cal_single(self, starts):
        result = np.zeros(self.frames)
        for start in starts:
            print(start)
            for frame in range(self.frames):
                result[frame] += self.calculator.get_correlation(start, start + frame * self.steps)
        return result

    def cal_single_shared(self, starts, shared_val, i, return_dict):
        result = np.zeros(self.frames)
        for j, start in enumerate(starts):
            print(j)
            for frame in range(self.frames):
                result[frame] += self.calculator.get_correlation_shared(start, start + frame * self.steps, shared_val)
        return_dict[i] = result

    def cal_parallel(self, parallel=24, p2=5):
        if not self.frame_set:
            for start in self.starts:
                self.frame_set.update(set(np.arange(self.frames) * self.steps + start))
            pass
        self.calculator.pre_calculation_parallel(self.frame_set, parallel=parallel)
        print("finish pre calculation")
        starts = list(self.starts)
        if parallel < 1:
            result = self.cal_single(starts)
        else:
            pool = Pool(p2)
            ret = list()
            seg = len(starts) // p2
            for i in range(p2):
                if i < p2 - 1:
                    ret.append(pool.apply_async(self.cal_single, args=(starts[i * seg:(i + 1) * seg],)))
                else:
                    ret.append(pool.apply_async(self.cal_single, args=(starts[i * seg:],)))
            pool.close()
            pool.join()
            result = np.zeros(self.frames)
            for r in ret:
                result += r.get()
        return result / len(self.starts)

    def cal_parallel_shared(self, val, parallel=24):
        # shared_val = self.calculator.pre_calculation_shared(val)
        manager = Manager()
        shared_val = manager.Array("f", val)
        return_dict = manager.dict()
        print("finish 2")

        starts = list(self.starts)
        seg = len(starts) // parallel
        jobs = list()
        for i in range(parallel):
            if i < parallel - 1:
                p = Process(target=self.cal_single_shared,
                            args=(starts[i * seg:(i + 1) * seg], shared_val, i, return_dict))
            else:
                p = Process(target=self.cal_single_shared, args=(starts[i * seg:], shared_val, i, return_dict))
            jobs.append(p)
        for p in jobs:
            p.start()
            p.join()

        result = np.zeros(self.frames)
        # print(return_dict)
        for r in return_dict.values():
            # print(r)
            result += r
        return result / len(self.starts)


class Urea:
    def __init__(self, c1_id):
        self.c1 = c1_id
        self.o1 = 0
        self.n1 = list()
        self.h1 = list()
        # c1, o1, n11, n12, h11, h12
        self.co_matrix = np.zeros((6, 3))
        self.idx = -1

    def get_all_id_except_c(self):
        return [self.o1] + self.n1 + self.h1

    def get_copy(self):
        nu = Urea(self.c1)
        nu.o1 = self.o1
        nu.h1 = self.h1[:]
        nu.n1 = self.n1[:]
        nu.co_matrix = self.co_matrix.copy()
        nu.idx = self.idx
        return nu

    def set_co(self, atom3d: csc.Atom3D):
        self.co_matrix[0] = atom3d.Atoms[self.c1].coordinate
        self.co_matrix[1] = atom3d.Atoms[self.o1].coordinate
        self.co_matrix[2] = atom3d.Atoms[self.n1[0]].coordinate
        self.co_matrix[3] = atom3d.Atoms[self.n1[1]].coordinate
        self.co_matrix[4] = atom3d.Atoms[self.h1[0]].coordinate
        self.co_matrix[5] = atom3d.Atoms[self.h1[1]].coordinate


class UreaAnalyzer:
    def __init__(self, atom3d: csc.Atom3D):
        self.ureas = list()
        self.create_ureas(atom3d)
        # self.only

    def create_ureas(self, a3d):
        i = 0
        for atom_id, atom in a3d.Atoms.items():
            if atom.type == "c1":
                u = Urea(atom_id)
                u.idx = i
                i += 1
                connects = a3d.get_connected_atom_id(atom_id)
                for a_id in connects:
                    a = a3d.Atoms[a_id]
                    if a.type == "o1":
                        u.o1 = a_id
                    elif a.type == "n1":
                        u.n1.append(a_id)
                        con = a3d.get_connected_atom_id(a_id)
                        for n_id in con:
                            nc = a3d.Atoms[n_id]
                            if nc.type == "h1":
                                u.h1.append(n_id)
                self.ureas.append(u)

    def convert_u_in_cell(self, atom3d: csc.Atom3D):
        for u in self.ureas:
            c1 = atom3d.Atoms[u.c1]
            c1.set_in_cell_coordinate(atom3d.lattice_parameter, atom3d.box_l)
            for a_id in u.get_all_id_except_c():
                a = atom3d.Atoms[a_id]
                a.set_image(c1.image, atom3d.lattice_parameter)
        return atom3d

    def set_u_co(self, atom3d: csc.Atom3D):
        for u in self.ureas:
            u.set_co(atom3d)

    # must convert_u_in_cell() and set_u_co() first
    def get_extended_u(self, atom3d: csc.Atom3D, extend_range=6.):
        def get_moved_u(_u, x, y, z):
            new_u = _u.get_copy()
            d = np.asarray([x, y, z], dtype=int)
            new_u.co_matrix -= d * atom3d.lattice_parameter
            return new_u

        extended_list = list()
        upper = atom3d.box_h - extend_range
        lower = atom3d.box_l + extend_range
        for u in self.ureas:
            c1 = atom3d.Atoms[u.c1]
            record = (c1.coordinate > upper).astype(int) + (c1.coordinate <= lower).astype(int) * -1
            if record[0] != 0:
                extended_list.append(get_moved_u(u, record[0], 0, 0))
                if record[1] != 0:
                    extended_list.append(get_moved_u(u, record[0], record[1], 0))
                    if record[2] != 0:
                        extended_list.append(get_moved_u(u, record[0], record[1], record[2]))
                if record[2] != 0:
                    extended_list.append(get_moved_u(u, record[0], 0, record[2]))
            if record[1] != 0:
                extended_list.append(get_moved_u(u, 0, record[1], 0))
                if record[2] != 0:
                    extended_list.append(get_moved_u(u, 0, record[1], record[2]))
            if record[2] != 0:
                extended_list.append(get_moved_u(u, 0, 0, record[2]))
        return extended_list

    @staticmethod
    def cal_uu_carbonyl_hbond(u1: Urea, u2: Urea, h_len=3.0, h_ang=120):
        def hbond_judge(o1_co, h1_co, n1_co):
            hb_len = cso.calculate_distance(o1_co, h1_co)
            if hb_len > h_len:
                return False
            else:
                hb_ang = cso.calculate_angle(o1_co, h1_co, n1_co)
                if hb_ang < h_ang:
                    return False
                else:
                    return True

        h_sum = 0
        if hbond_judge(u1.co_matrix[1], u2.co_matrix[4], u2.co_matrix[2]):
            h_sum += 1
        if hbond_judge(u1.co_matrix[1], u2.co_matrix[5], u2.co_matrix[3]):
            h_sum += 1
        return h_sum

    @staticmethod
    def hbond_judge(o1_co, h1_co, n1_co, h_len=3.0, h_ang=120):
        hb_len = cso.calculate_distance(o1_co, h1_co)
        if hb_len > h_len:
            return False
        else:
            hb_ang = cso.calculate_angle(o1_co, h1_co, n1_co)
            if hb_ang < h_ang:
                return False
            else:
                return True

    @staticmethod
    def cal_uu_hbond(u1: Urea, u2: Urea, h_len=3.0, h_ang=120):
        h_sum = 0
        if UreaAnalyzer.hbond_judge(u1.co_matrix[1], u2.co_matrix[4], u2.co_matrix[2], h_len, h_ang):
            h_sum += 1
        if UreaAnalyzer.hbond_judge(u1.co_matrix[1], u2.co_matrix[5], u2.co_matrix[3], h_len, h_ang):
            h_sum += 1
        if UreaAnalyzer.hbond_judge(u2.co_matrix[1], u1.co_matrix[4], u1.co_matrix[2], h_len, h_ang):
            h_sum += 1
        if UreaAnalyzer.hbond_judge(u2.co_matrix[1], u1.co_matrix[5], u1.co_matrix[3], h_len, h_ang):
            h_sum += 1
        if h_sum < 2:
            return h_sum
        else:
            return 2

    def cal_all_uu_vec(self, atom3d: csc.Atom3D, search_range=6.0):
        self.convert_u_in_cell(atom3d)
        self.set_u_co(atom3d)
        central_u = self.ureas
        extended_u = self.get_extended_u(atom3d)
        central_u_co = list(map(lambda u: u.co_matrix[0], central_u))
        # extended_u_co = list(map(lambda u: u.co_matrix[0], extended_u))
        all_u = central_u + extended_u
        all_u_co = list(map(lambda u: u.co_matrix[0], all_u))
        tree = cKDTree(all_u_co)
        uu_h1 = list()
        uu_h2 = list()
        for i, u_co in enumerate(central_u_co):
            u1 = central_u[i]
            for idx in tree.query_ball_point(u_co, search_range):
                if idx == i:
                    continue
                u2 = all_u[idx]
                h_num = self.cal_uu_hbond(u1, u2)
                if h_num == 1:
                    uu_h1.append(u1.co_matrix[0] - u2.co_matrix[0])
                elif h_num == 2:
                    uu_h2.append(u1.co_matrix[0] - u2.co_matrix[0])
        return uu_h1, uu_h2

    def cal_all_uu_carbonyl(self, atom3d: csc.Atom3D, search_range=6.0):
        self.convert_u_in_cell(atom3d)
        self.set_u_co(atom3d)
        central_u = self.ureas
        extended_u = self.get_extended_u(atom3d)
        central_u_co = list(map(lambda u: u.co_matrix[0], central_u))
        # extended_u_co = list(map(lambda u: u.co_matrix[0], extended_u))
        all_u = central_u + extended_u
        all_u_co = list(map(lambda u: u.co_matrix[0], all_u))
        tree = cKDTree(all_u_co)
        co_h = [[], [], []]
        for i, u_co in enumerate(central_u_co):
            u1 = central_u[i]
            h_num = 0
            for idx in tree.query_ball_point(u_co, search_range):
                if idx == i:
                    continue
                u2 = all_u[idx]
                h_num += self.cal_uu_carbonyl_hbond(u1, u2)
            if h_num > 2:
                print(h_num)
                h_num = 2
            co_h[h_num].append(u1.co_matrix[1] - u1.co_matrix[0])
        return co_h

    @staticmethod
    def uu_statistics(uu_h1, uu_h2, direction=np.array([1, 0, 0])):
        uu_num = [len(uu_h1) / 2, len(uu_h2) / 2]
        uu_d = [0, 0]
        d_h1 = 0.
        for uu_vec in uu_h1:
            d_h1 += np.square(np.dot(uu_vec, direction)) / np.dot(uu_vec, uu_vec)
        uu_d[0] += d_h1 / uu_num[0] / 2
        d_h2 = 0.
        for uu_vec in uu_h2:
            d_h2 += np.square(np.dot(uu_vec, direction)) / np.dot(uu_vec, uu_vec)
        uu_d[1] += d_h2 / uu_num[1] / 2
        return uu_num, uu_d

    @staticmethod
    def carbonyl_statistics(co_h, direction=np.array([1, 0, 0])):
        # print(co_h)
        co_num = list()
        co_direction = list()
        for i in range(3):
            co_num.append(len(co_h[i]))
            d = 0.
            for co_vec in co_h[i]:
                d += np.square(np.dot(co_vec, direction)) / np.dot(co_vec, co_vec)
            co_direction.append(d / co_num[i])
        return co_num, co_direction

    def cal_all_bonds(self, frames, atom3ds):
        co_dict = dict()
        nh_dict = dict()
        for urea in self.ureas:
            co_dict[(urea.c1, urea.o1)] = []
            nh_dict[(urea.n1[0], urea.h1[0])] = []
            nh_dict[(urea.n1[1], urea.h1[1])] = []
        for frame in frames:
            if frame % 100 == 0:
                print(frame)
            atom3d = atom3ds(frame)
            atom3d.convert_to_default()
            for co in co_dict:
                co_dict[co].append(
                    cso.calculate_distance(atom3d.Atoms[co[0]].coordinate, atom3d.Atoms[co[1]].coordinate))
            for nh in nh_dict:
                nh_dict[nh].append(
                    cso.calculate_distance(atom3d.Atoms[nh[0]].coordinate, atom3d.Atoms[nh[1]].coordinate))
        return co_dict, nh_dict

    def visit_all_uu(self, atom3d: csc.Atom3D, visit, search_range=6.0, extended=True):
        self.convert_u_in_cell(atom3d)
        self.set_u_co(atom3d)
        central_u = self.ureas
        if extended:
            extended_u = self.get_extended_u(atom3d)
            # extended_u_co = list(map(lambda u: u.co_matrix[0], extended_u))
            all_u = central_u + extended_u
        else:
            all_u = central_u
        central_u_co = list(map(lambda u: u.co_matrix[0], central_u))
        all_u_co = list(map(lambda u: u.co_matrix[0], all_u))
        tree = cKDTree(all_u_co)
        for i, u_co in enumerate(central_u_co):
            u1 = central_u[i]
            for idx in tree.query_ball_point(u_co, search_range):
                if idx == i:
                    continue
                u2 = all_u[idx]
                visit(u1, u2)


class UUVisitor:

    def visit(self, u1: Urea, u2: Urea):
        pass

    def get_result(self):
        return 0.

    def clear(self):
        pass


class UUStress(UUVisitor):

    def __init__(self, force, atom3d: csc.Atom3D):
        self.force = force
        self.stress_sum = np.zeros(3)
        self.atom3d = atom3d
        self.convert = 4.18585 / 6.02 / 1.01325 * np.power(10, 2 - 12 + 20 - 5)

    def visit(self, u1, u2):
        c1 = self.atom3d.Atoms[u1.c1]
        c2 = self.atom3d.Atoms[u2.c1]
        if c1.block == c2.block:
            return 0
        types = ["c1", "o1", "n1", "n1", "h1", "h1"]
        for i in range(6):
            for j in range(6):
                type_tuple = tuple(sorted([types[i], types[j]]))
                dr = u2.co_matrix[j] - u1.co_matrix[i]
                dist = cso.calculate_distance(dr)
                f = self.force.cal_force_tt(dist, type_tuple)
                self.stress_sum -= np.square(dr) * f / dist

    def get_result(self):
        return self.stress_sum.copy() * self.convert / 2

    def clear(self):
        self.stress_sum = np.zeros(3)


class UUForce(UUVisitor):

    def __init__(self, force, atom3d: csc.Atom3D):
        self.force = force
        self.tangential = list()
        self.axial = list()
        self.atom3d = atom3d
        self.convert = 4.18585 / 6.02 * 100

    def visit(self, u1, u2):
        c1 = self.atom3d.Atoms[u1.c1]
        c2 = self.atom3d.Atoms[u2.c1]
        if c1.block == c2.block:
            return 0
        dc = u2.co_matrix[0] - u1.co_matrix[0]
        v_axial = dc / cso.calculate_distance(dc)
        dn = u1.co_matrix[2] - u1.co_matrix[3]
        v_tange = dn / cso.calculate_distance(dn)
        types = ["c1", "o1", "n1", "n1", "h1", "h1"]
        f_axial = 0.
        f_tange = np.zeros(3)
        for i in range(6):
            for j in range(6):
                type_tuple = tuple(sorted([types[i], types[j]]))
                dr = -(u2.co_matrix[j] - u1.co_matrix[i])
                dist = cso.calculate_distance(dr)
                f_ = self.force.cal_force_tt(dist, type_tuple)
                f_ = dr * f_ / dist
                # print(f)
                f_a = np.inner(f_, v_axial)
                f_t = f_ - f_a * v_axial
                f_axial += f_a
                f_tange += f_t
        self.axial.append(f_axial)
        self.tangential.append(cso.calculate_distance(f_tange))

    def get_result(self):
        return np.asarray(self.axial) * self.convert, np.asarray(self.tangential) * self.convert

    def clear(self):
        self.axial = list()
        self.tangential = list()


# class UUStressDecompose(UUVisitor):
#
#     def __init__(self, ff, atom3d: csc.Atom3D):
#         self.ff = ff
#         self.tangential = list()
#         self.axial = list()
#         self.atom3d = atom3d
#         self.convert = 4.18585 / 6.02 * 100
#
#     def visit(self, u1, u2):
#         c1 = self.atom3d.Atoms[u1.c1]
#         c2 = self.atom3d.Atoms[u2.c1]
#         if c1.block == c2.block:
#             return 0
#         dc = u2.co_matrix[0] - u1.co_matrix[0]
#         v_axial = dc / cso.calculate_distance(dc)
#         dn = u1.co_matrix[2] - u1.co_matrix[3]
#         v_tange = dn / cso.calculate_distance(dn)
#         types = ["c1", "o1", "n1", "n1", "h1", "h1"]
#         f_axial = 0.
#         f_tange = 0.
#         for i in range(6):
#             for j in range(6):
#                 type_tuple = tuple(sorted([types[i], types[j]]))
#                 dr = -(u2.co_matrix[j] - u1.co_matrix[i])
#                 dist = cso.calculate_distance(dr)
#                 f_ = self.ff.cal_force_tt(dist, type_tuple)
#                 f_ = dr * f_ / dist
#                 # print(f)
#                 f_a = np.inner(f_, v_axial)
#                 f_t = f_ - f_a * v_axial
#                 f_axial += f_a
#                 f_tange += cso.calculate_distance(f_t)
#         self.axial.append(f_axial)
#         self.tangential.append(f_tange)
#
#     def get_result(self):
#         return np.asarray(self.axial) * self.convert, np.asarray(self.tangential) * self.convert
#
#     def clear(self):
#         self.axial = list()
#         self.tangential = list()


class UUPotential(UUVisitor):

    def __init__(self, force, atom3d: csc.Atom3D):
        self.force = force
        self.potential = list()
        self.atom3d = atom3d
        self.convert = 1

    def visit(self, u1, u2):
        c1 = self.atom3d.Atoms[u1.c1]
        c2 = self.atom3d.Atoms[u2.c1]
        if c1.block == c2.block:
            return 0
        types = ["c1", "o1", "n1", "n1", "h1", "h1"]
        potential = 0.
        for i in range(6):
            for j in range(6):
                type_tuple = tuple(sorted([types[i], types[j]]))
                dr = u2.co_matrix[j] - u1.co_matrix[i]
                dist = cso.calculate_distance(dr)
                pe = self.force.cal_energy_tt(dist, type_tuple)
                potential += pe
        self.potential.append(potential)

    def get_result(self):
        return np.asarray(self.potential) * self.convert

    def clear(self):
        self.potential = list()


class EndToEndDistance:
    def __init__(self, atom3d_gen: Atom3dGenerator, chain_len, chain_num, ns):
        self.atom3d_gen = atom3d_gen
        self.ee_dict = self.get_ee_dict(chain_len, chain_num, ns)
        self.ns = ns

    def get_ee_dict(self, chain_len, chain_num, ns):
        atom3d = self.atom3d_gen()
        ee_dict = dict()
        for n in ns:
            ee_dict[n] = list()
            if n == 0 or n >= chain_len - 1:
                for i in range(chain_num):
                    ee_dict[n].append((i * chain_len + 1, (i + 1) * chain_len))
            else:
                for i in range(chain_num):
                    ee_dict[n].append((i * chain_len + 1, i * chain_len + 1 + n))
                    ee_dict[n].append(((i + 1) * chain_len - n, (i + 1) * chain_len))
        return ee_dict

    def cal_dists(self, frames):
        dists = dict()
        for n in self.ns:
            dists[n] = list()
        # print(self.ee_dict)
        for frame in frames:
            print(frame)
            atom3d = self.atom3d_gen(frame)
            # print(atom3d)
            atom3d.convert_to_default()
            for n in self.ns:
                for pair in self.ee_dict[n]:
                    # print(atom3d.Atoms[pair[0]])
                    # print(atom3d.Atoms[pair[1]])
                    dists[n].append(
                        cso.calculate_distance(atom3d.Atoms[pair[0]].coordinate, atom3d.Atoms[pair[1]].coordinate))
        return dists


class EndToEndDistanceType:
    def __init__(self, atom3d_gen: Atom3dGenerator):
        self.atom3d_gen = atom3d_gen
        self.ee_dict = self.get_ee_dict()

    def get_ee_dict(self):
        atom3d = self.atom3d_gen()
        ee_dict = dict()
        for atom in atom3d.Atoms.values():
            if atom.type == "TO(1)":
                if atom.mole_id in ee_dict:
                    ee_dict[atom.mole_id].append(atom.id)
                else:
                    ee_dict[atom.mole_id] = [atom.id]
        return ee_dict

    def cal_dists(self, frames):
        dists = list()
        for frame in frames:
            # print(frame)
            atom3d = self.atom3d_gen(frame)
            # print(atom3d)
            atom3d.convert_to_default()
            for pair in self.ee_dict.values():
                # print(atom3d.Atoms[pair[0]])
                # print(atom3d.Atoms[pair[1]])
                dists.append(cso.calculate_distance(atom3d.Atoms[pair[0]].coordinate, atom3d.Atoms[pair[1]].coordinate))
        return dists

    def cal_dists_parallel(self, frames, parallel=24):
        pool = Pool(parallel)
        dists_all = list()
        seg = len(frames) // parallel
        ret = list()
        for i in range(parallel):
            if i < parallel - 1:
                ret.append(pool.apply_async(self.cal_dists, args=(frames[i * seg:(i + 1) * seg],)))
            else:
                ret.append(pool.apply_async(self.cal_dists, args=(frames[i * seg:],)))
        pool.close()
        pool.join()
        for r in ret:
            dists_all += r.get()
        return dists_all


class SAXS:
    def __init__(self, atom3d_gen: Atom3dGenerator):
        self.atom3ds = atom3d_gen

    def density_matrix(self, time_step, n_xyz, grid=1., ex_soft=False):
        atom3d = self.atom3ds(time_step)
        atom3d.convert_to_in_cell()
        xa, ya, za = atom3d.lattice_parameter
        xl, yl, zl = atom3d.box_l
        n_xa, n_ya, n_za = n_xyz
        xxa = np.linspace((xa - n_xa * grid) / 2 + xl, xa - (xa - n_xa * grid) / 2 + xl, n_xa + 1)
        yya = np.linspace((ya - n_ya * grid) / 2 + yl, ya - (ya - n_ya * grid) / 2 + yl, n_ya + 1)
        zza = np.linspace((za - n_za * grid) / 2 + zl, za - (za - n_za * grid) / 2 + zl, n_za + 1)
        H_sum = np.zeros((n_xa, n_ya, n_za))
        # TODO consider period boundary condition
        for atom_type in self.atom3ds().Atom_type_dict:
            # exclude soft for better contrast
            if ex_soft and (atom_type[:2] == "TO" or atom_type[:2] == "Es"):
                continue
            coordinate_list = list()
            for a_id in self.atom3ds().Atom_type_dict[atom_type]:
                # print(atom3d.Atoms[a_id].coordinate)
                coordinate_list.append(atom3d.Atoms[a_id].coordinate)
            coordinate_list = np.vstack(coordinate_list)

            H, bin_edges = np.histogramdd(coordinate_list, bins=(xxa, yya, zza))
            # print(bin_edges)
            H_sum += H * atom3d.type_mass_dict[atom_type]
        return H_sum

    def density_multi_step(self, time_steps, grid=1.):
        atom3d = self.atom3ds(time_steps[0])
        xa, ya, za = atom3d.lattice_parameter
        n_xyz = [int(xa / grid), int(ya / grid), int(za / grid)]
        H_sum = self.density_matrix(time_steps[0], n_xyz, grid)
        for time_step in time_steps[1:]:
            print(time_step)
            H_sum += self.density_matrix(time_step, n_xyz, grid)
        H_sum /= len(time_steps)
        return H_sum

    @staticmethod
    def sum_box_density(H, sum_size):
        assert H.shape[0] % sum_size == 0
        number = int(H.shape[0] / sum_size)
        new_H = np.split(H[:int(number * sum_size)], number, axis=0)
        new_H = np.sum(new_H, axis=1)
        new_H = np.split(new_H[:int(number * sum_size)], number, axis=1)
        new_H = np.sum(new_H, axis=2)
        new_H = np.split(new_H[:int(number * sum_size)], number, axis=2)
        new_H = np.sum(new_H, axis=3)
        return new_H

    @staticmethod
    def random_sub_box(H, num_grid):
        assert num_grid < min(H.shape), "num_grid must less than %d" % min(H.shape)
        xb = np.random.randint(0, H.shape[0] - num_grid, dtype=int)
        yb = np.random.randint(0, H.shape[1] - num_grid, dtype=int)
        zb = np.random.randint(0, H.shape[2] - num_grid, dtype=int)
        return H[xb:xb + num_grid, yb:yb + num_grid, zb:zb + num_grid].copy()

    @staticmethod
    def fft3d(H, sum_size):
        fft_H = np.abs(np.fft.fftn(H)) / np.power(sum_size * H.shape[0], 3)
        fft_freq = np.fft.fftfreq(H.shape[0], d=sum_size) * 2 * np.pi
        return fft_freq, fft_H

    @staticmethod
    def fft2d(H, sum_axis):
        H = np.sum(H, axis=sum_axis)
        fft_H = np.abs(np.fft.fftn(H)) / np.power(H.shape[0], 3)
        fft_freq = np.fft.fftfreq(H.shape[0]) * 2 * np.pi
        return np.fft.fftshift(fft_freq), np.fft.fftshift(fft_H)

    @staticmethod
    def fft1d(H, sum_axis):
        H = np.sum(H, axis=sum_axis[0])
        H = np.sum(H, axis=sum_axis[1] - 1)
        fft_H = np.abs(np.fft.fftn(H)) / np.power(H.shape[0], 3)
        fft_freq = np.fft.fftfreq(H.shape[0]) * 2 * np.pi
        return np.fft.fftshift(fft_freq), np.fft.fftshift(fft_H)

    def fft3d_aver(self, H, num_grid, sum_size, num=10):
        H_sub = self.random_sub_box(H, num_grid)
        H_sub = self.sum_box_density(H_sub, sum_size)
        fft_freq, fft_h = self.fft3d(H_sub, sum_size)
        for i in range(num - 1):
            H_sub = self.random_sub_box(H, num_grid)
            H_sub = self.sum_box_density(H_sub, sum_size)
            fft_freq_1, fft_h_1 = self.fft3d(H_sub, sum_size)
            fft_h += fft_h_1
        return fft_freq, fft_h / num

    def fft2d_aver(self, H, num_grid, num=10, sum_axis=1):
        H_sub = self.random_sub_box(H, num_grid)
        fft_freq, fft_h = self.fft2d(H_sub, sum_axis=sum_axis)
        for i in range(num - 1):
            H_sub = self.random_sub_box(H, num_grid)
            fft_freq_1, fft_h_1 = self.fft2d(H_sub, sum_axis=sum_axis)
            fft_h += fft_h_1
        return fft_freq, fft_h / num

    def fft1d_aver(self, H, num_grid, num=10, sum_axis=(1, 2)):
        H_sub = self.random_sub_box(H, num_grid)
        fft_freq, fft_h = self.fft1d(H_sub, sum_axis=sum_axis)
        for i in range(num - 1):
            H_sub = self.random_sub_box(H, num_grid)
            fft_freq_1, fft_h_1 = self.fft1d(H_sub, sum_axis=sum_axis)
            fft_h += fft_h_1
        return fft_freq, fft_h / num

    @staticmethod
    def cylinder_average(freqs, h_3ds, out_grid=0.02, out_max=0.4):
        grid_num = int(out_max / out_grid)
        record_sum = np.zeros((grid_num, grid_num))
        record_num = np.zeros((grid_num, grid_num))
        for _ in range(len(h_3ds)):
            h_3d = h_3ds[_]
            freq = freqs[_]
            for i in range(h_3d.shape[0]):
                for j in range(h_3d.shape[1]):
                    for k in range(h_3d.shape[2]):
                        x = freq[i]
                        y = freq[j]
                        z = freq[k]
                        yz = np.sqrt(y * y + z * z)
                        ix = int(abs(x) / out_grid)
                        iyz = int(yz / out_grid)
                        if ix >= grid_num or iyz >= grid_num:
                            continue
                        record_sum[ix, iyz] += h_3d[i, j, k]
                        record_num[ix, iyz] += 1
        record_num += np.where(record_num == 0, 1, 0)
        return record_sum / record_num

    def fft_multi_aver_3d(self, time_steps, random_num=10, out_grid=0.02, out_max=0.4):
        H = self.density_multi_step(time_steps)
        sum_size = 1
        all_freq = list()
        all_fft = list()
        for num_grid in range(80, 165, 5):
            fft_freq, fft_h = self.fft3d_aver(H, num_grid, sum_size, num=random_num)
            all_freq.append(fft_freq)
            all_fft.append(fft_h)
        saxs = self.cylinder_average(all_freq, all_fft, out_grid, out_max)
        return saxs

    # sum_size = 1
    def fft_multi_aver_2d(self, time_steps, random_num=10, out_grid=0.02, out_max=0.4, num_grid=180):
        print("calculate density matrix")
        H = self.density_multi_step(time_steps)
        print(H.shape)
        # num_grid = 150
        
        # num_grid = 300
        print("y direction")
        fft_freq_1, fft_h_1 = self.fft2d_aver(H, num_grid, num=random_num, sum_axis=1)
        print("z direction")
        fft_freq_2, fft_h_2 = self.fft2d_aver(H, num_grid, num=random_num, sum_axis=2)
        fft_h = (fft_h_1 + fft_h_2) / 2
        l = np.where(fft_freq_1 > -1 * out_max)[0][0]
        r = np.where(fft_freq_1 > out_max)[0][0]

        # all_freq = list()
        # all_fft = list()
        # for num_grid in range(80, 165, 5):
        #     fft_freq, fft_h = self.fft3d_aver(H, num_grid, sum_size, num=random_num)
        #     all_freq.append(fft_freq)
        #     all_fft.append(fft_h)
        # saxs = self.cylinder_average(all_freq, all_fft, out_grid, out_max)
        return fft_freq_1[l:r], fft_h[l:r, l:r]

    def fft_multi_aver_1d(self, time_steps, random_num=10, out_grid=0.02, out_max=0.4):
        print("calculate density matrix")
        H = self.density_multi_step(time_steps)
        print(H.shape)
        # num_grid = 150
        fft_freq = list()
        fft_h = list()
        num_grids = np.arange(0.045, 0.2, 0.005)
        num_grids = 1 / (num_grids / 2 / np.pi)
        num_grids = num_grids.astype(int)
        for num_grid in num_grids:
            print("yz direction")
            fft_freq_1, fft_h_1 = self.fft1d_aver(H, num_grid, num=random_num, sum_axis=(1, 2))
            # fft_freq += list(fft_freq_1)
            # fft_h += list(fft_h_1)
            fft_freq.append(fft_freq_1)
            fft_h.append(fft_h_1)
        #
        # fft_freq = np.asarray(fft_freq)
        # fft_h = np.asarray(fft_h)
        # ind = np.argsort(fft_freq)
        # fft_freq = np.sort(fft_freq)
        # fft_h = np.take_along_axis(fft_h, ind, axis=0)
        #
        # l = np.where(fft_freq > -1 * out_max)[0][0]
        # r = np.where(fft_freq > out_max)[0][0]
        #
        # return fft_freq[l:r], fft_h[l:r]
        return fft_freq, fft_h


class GenerateFrames:
    def __init__(self, working_path):
        self.wp = working_path
        self.clear_all()
        os.makedirs(self.wp, exist_ok=True)

    def set_up(self, start, steps, devery):
        pass

    def generate_trj(self, num_parallel=24, screen=False):
        p_lmp = Process(target=fps.invoke_lmp, args=(self.wp, "in.gen", num_parallel, screen))
        p_lmp.start()
        p_lmp.join(timeout=100)
        return self.wp

    def clear_all(self):
        os.system("rm -r %s" % self.wp)


class GenerateFramesStrain(GenerateFrames):
    def __init__(self, working_path, atom3d_gen: Atom3dGenerator, strain_rate):
        super().__init__(working_path)
        self.atom3ds = atom3d_gen
        self.srate = strain_rate

    def set_up(self, start, steps, devery):
        super().set_up(start, steps, devery)
        # prepare files
        data_path = self.atom3ds.get_name()
        origin_path = data_path[:data_path.rfind("/") + 1]
        trj_path = self.atom3ds.get_name(start)
        os.system("cp %s %s" % (data_path, self.wp))
        os.system("cp %s %s" % (trj_path, self.wp))
        os.system("cp %s*.param %s" % (origin_path, self.wp))

        # prepare in file
        template_path = "/home/centos/Projects/CGMD/template/in.strain"
        data_name = data_path[data_path.rfind("/") + 1:]
        trj_name = trj_path[trj_path.rfind("/") + 1:]
        with open(template_path, "r") as file:
            lines = file.readlines()
        v_lines = ["variable    data_name string %s\n" % data_name,
                   "variable    trj_name string %s\n" % trj_name,
                   "variable    trj_step equal %d\n" % start,
                   "variable	srate equal %.8f\n" % (self.srate / (1 + self.srate * start)),
                   "variable 	steps equal %d\n" % steps,
                   "variable    devery equal %d\n" % devery]
        with open(self.wp + "in.gen", "w") as file:
            file.writelines(v_lines + lines)


class ImageHelper:
    def __init__(self, image, d=0, pbc=True):
        self.image = image
        self.pbc = pbc
        self.eimage = self.expand()
        self.d = d

    def expand(self):
        newi = np.ones((self.image.shape[0] + 2, self.image.shape[1] + 2, self.image.shape[2] + 2))
        newi[1:-1, 1:-1, 1:-1] = self.image
        if self.pbc:
            newi[1:-1, 1:-1, 0] = self.image[:, :, -1]
            newi[1:-1, 1:-1, -1] = self.image[:, :, 0]
            newi[1:-1, 0, 1:-1] = self.image[:, -1, :]
            newi[1:-1, -1, 1:-1] = self.image[:, 0, :]
            newi[0, 1:-1, 1:-1] = self.image[-1, :, :]
            newi[-1, 1:-1, 1:-1] = self.image[0, :, :]
            newi[0, 0, 1:-1] = self.image[-1, -1, :]
            newi[0, -1, 1:-1] = self.image[-1, 0, :]
            newi[-1, 0, 1:-1] = self.image[0, -1, :]
            newi[-1, -1, 1:-1] = self.image[0, 0, :]
            newi[0, 1:-1, 0] = self.image[-1, :, -1]
            newi[0, 1:-1, -1] = self.image[-1, :, 0]
            newi[-1, 1:-1, -1] = self.image[0, :, 0]
            newi[-1, 1:-1, 0] = self.image[0, :, -1]
            newi[1:-1, 0, 0] = self.image[:, -1, -1]
            newi[1:-1, -1, 0] = self.image[:, 0, -1]
            newi[1:-1, -1, -1] = self.image[:, 0, 0]
            newi[1:-1, 0, -1] = self.image[:, -1, 0]
            newi[0, 0, 0] = self.image[-1, -1, -1]
            newi[0, 0, -1] = self.image[-1, -1, 0]
            newi[0, -1, -1] = self.image[-1, 0, 0]
            newi[0, -1, 0] = self.image[-1, 0, -1]
            newi[-1, 0, 0] = self.image[0, -1, -1]
            newi[-1, 0, -1] = self.image[0, -1, 0]
            newi[-1, -1, -1] = self.image[0, 0, 0]
            newi[-1, -1, 0] = self.image[0, 0, -1]
        return newi

    def ix(self, a, b, c):
        if self.d == 0:
            return self.image[a, b, c]
        elif self.d == 1:
            return self.image[c, a, b]
        else:
            return self.image[b, c, a]

    def mat(self, a, b, c):
        if self.d == 0:
            return self.eimage[a:a + 3, b:b + 3, c:c + 3]
        elif self.d == 1:
            return self.eimage[c:c + 3, a:a + 3, b:b + 3]
        else:
            return self.eimage[b:b + 3, c:c + 3, a:a + 3]

    def set_i(self, a, b, c, v):
        if self.d == 0:
            self.image[a, b, c] = v
        elif self.d == 1:
            self.image[c, a, b] = v
        else:
            self.image[b, c, a] = v

    def set_e(self, a, b, c, v):
        if self.d == 0:
            self.eimage[a + 1, b + 1, c + 1] = v
        elif self.d == 1:
            self.eimage[c + 1, a + 1, b + 1] = v
        else:
            self.eimage[b + 1, c + 1, a + 1] = v


class SegmentJudger:
    @staticmethod
    def load_seg_id(_fp):
        mole_seg = dict()
        with open(_fp, "r") as file:
            seg_lines = file.readlines()
        for line in seg_lines:
            if line:
                if line[0] != "#":
                    args = line.split(" ")
                    mole_seg[int(args[0])] = int(args[1])
        return mole_seg

    # if in the same segment, return True, else return False
    def judge(self, atom1: bsc.Atom, atom2: bsc.Atom, image):
        return True


class SJPara(SegmentJudger):
    def __init__(self, mole_seg, _type="x"):
        if isinstance(mole_seg, dict):
            self.mole_seg = mole_seg
        else:
            self.mole_seg = self.load_seg_id(mole_seg)
        self.type = _type

    def judge(self, atom1: bsc.Atom, atom2: bsc.Atom, image):
        mol1, mol2 = atom1.mole_id, atom2.mole_id
        seg1, seg2 = self.mole_seg[mol1], self.mole_seg[mol2]
        if self.type == "x":
            if seg1 == seg2 and image[1] == 0 and image[2] == 0:
                return True
        else:
            raise KeyError("type not defined")
        return False


class Graph:
    class Node:
        def __init__(self, nid, value, *eid):
            self.n_id = nid
            self.value = value
            self.eid = list(eid)
            # self.visited = False

    class Edge:
        def __init__(self, eid, value, nid1, nid2):
            self.eid = eid
            self.value = value
            self.nids = {nid1:nid2, nid2:nid1}
            # self.visited = False

    def __init__(self):
        self.nodes = dict()
        self.edges = dict()

    def get_next(self, n_id, e_id):
        return self.edges[e_id].nids[n_id]

    def get_connected(self, n_id, level=1):
        connected = set()
        for e_id in self.nodes[n_id].eid:
            if self.edges[e_id].value >= level:
                connected.add(self.get_next(n_id, e_id))
        return connected

    # def clear_visit(self):
    #     for node in self.nodes:
    #         node.visited = False

    def traverse_bfs_single(self, b_id, visit, level=1):
        waiting = {b_id}
        visited = set()
        while len(waiting) > 0:
            v_id = waiting.pop()
            if v_id in visited:
                continue
            visit(v_id)
            visited.add(v_id)
            nw = self.get_connected(v_id, level=level)
            nw -= visited
            waiting = waiting | nw
        return visited


class UreaGraph(Graph):
    def create_from_h_dict(self, h_dict, h_cal):
        e_id = 0
        record_dict = dict()
        for d_id in h_dict:
            for a_id in h_dict[d_id]:
                c1, c2 = h_cal.d1_c1[d_id], h_cal.a_c1[a_id]
                h_tuple = tuple(sorted([c1, c2]))
                if h_tuple not in record_dict:
                    record_dict[h_tuple] = e_id
                    edge = Graph.Edge(e_id, 1, c1, c2)
                    self.edges[e_id] = edge
                    if c1 not in self.nodes:
                        self.nodes[c1] = Graph.Node(c1, 0, e_id)
                    else:
                        self.nodes[c1].eid.append(e_id)
                    if c2 not in self.nodes:
                        self.nodes[c2] = Graph.Node(c2, 0, e_id)
                    else:
                        self.nodes[c2].eid.append(e_id)
                    e_id += 1
                else:
                    self.edges[record_dict[h_tuple]].value += 1

    def connected_by_block(self, atom3d: csc.Atom3D):
        atom3d.set_block_id()
        block_record = dict()
        for atom in atom3d.Atoms.values():
            if atom.type == "c1":
                if atom.block in block_record:
                    block_record[atom.block].append(atom.id)
                else:
                    block_record[atom.block] = [atom.id]
        e_id = max(self.edges.keys()) + 1
        for pair in block_record.values():
            edge = Graph.Edge(e_id, 3, pair[0], pair[1])
            self.edges[e_id] = edge
            if pair[0] not in self.nodes:
                self.nodes[pair[0]] = Graph.Node(pair[0], 0, e_id)
            else:
                self.nodes[pair[0]].eid.append(e_id)
            if pair[1] not in self.nodes:
                self.nodes[pair[1]] = Graph.Node(pair[1], 0, e_id)
            else:
                self.nodes[pair[1]].eid.append(e_id)
            e_id += 1





if __name__ == '__main__':
    fp = "/home/centos/model/aa/1blk_40/0/"
    # a3d = csc.create_atom3d(fp, 40)
    # a3d.convert_to_in_cell()
    # LmpParser.create_data_file(a3d, "test0.data")
    # v = StaticAnalyzer.grid_hard_domain(a3d, r=3., out_path="test.data")
    # print(v)
