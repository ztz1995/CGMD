import numpy as np
from pyMD import collective_structure_class as csc
from pyMD import basic_structure_class as bsc
from pyMD import cartesian_operation as cso
from pyMD import functions as f
from pyMD.file_parser import LmpParser
from pyMD import ibi
import random
import copy as cp


# class HardBox(csc.Atom3D):
#     def __init__(self):
#         super().__init__()
#         self.es_pairs = list()
#         self.es_pairs_last = list()
#         self.chain_creator = None
#         self.super_to_mole = dict()
#         self.mole_to_super = dict()
#         self.super_to_pair = dict()
#         self.repeat_num = 0
#
#     def initialize(self, chain_num=204, repeat_num=4):
#         self.repeat_num = repeat_num
#         fp1 = "../data/20200420_lr=0.05_1blk_40_cg_ex_hard/record/"
#         r1 = ibi.Record(fp1, restart=True)
#         param1 = r1.get_param(22, 1)
#
#         pattern = ["TO(2)"] * 6 + ["Es", "Ph", "U", "Ph", "Me", "Ph", "U", "Ph", "Es"] + ["TO(2)"] * 7
#         pattern[0] = "TO(1)"
#         pattern[-1] = "TO(1)"
#
#         self.create_from_pattern_fix2(pattern, chain_num, param=param1, density=1., center=10, soft=False)
#
#         Es_ids = list()
#         mole_es = dict()
#         for a_id, atom in self.Atoms.items():
#             if atom.type == "Es":
#                 Es_ids.append(a_id)

#
#         random.shuffle(Es_ids)
#
#         pair_mole_id = set()
#         for _ in range(round(chain_num / 2)):
#             e1_id = Es_ids.pop()
#             while self.Atoms[e1_id].mole_id in pair_mole_id:
#                 e1_id = Es_ids.pop()
#             e2_id = Es_ids.pop()
#             while self.Atoms[e2_id].mole_id in pair_mole_id or self.Atoms[e2_id].mole_id == self.Atoms[e1_id].mole_id:
#                 e2_id = Es_ids.pop()
#             self.es_pairs.append((e1_id, e2_id))
#             pair_mole_id.add(self.Atoms[e1_id].mole_id)
#             pair_mole_id.add(self.Atoms[e2_id].mole_id)
#
#     def rule(self, p1, p2):
#         if self.Atoms[p1[0]].mole_id == self.Atoms[p1[1]].mole_id:
#             return False
#         if self.Atoms[p2[0]].mole_id == self.Atoms[p2[1]].mole_id:
#             return False
#         return True
#
#     def min_r2(self, r1, r2):
#         dr = r2 - r1
#         judge_p = np.where(dr > self.lattice_parameter / 2, 1, 0)
#         judge_n = np.where(dr < -self.lattice_parameter / 2, -1, 0)
#         judge = judge_n + judge_p
#         dr -= judge * self.lattice_parameter
#         return np.inner(dr, dr), r2 - judge * self.lattice_parameter
#
#     def pair_energy(self, pair, T, N=14, b=5.13):
#         r_2 = self.min_r2(self.Atoms[pair[0]].coordinate, self.Atoms[pair[1]].coordinate)[0]
#         return 1.5 * T * r_2 / N / b / b
#
#     def energy_dif(self, op1, op2, np1, np2, T):
#         return self.pair_energy(np1, T) + self.pair_energy(np2, T) - self.pair_energy(op1, T) - self.pair_energy(op2, T)
#
#     def move(self, T):
#         def new_pair(_p1, _p2):
#             if np.random.rand() < 0.5:
#                 _new_p1 = (_p1[0], _p2[0])
#                 _new_p2 = (_p1[1], _p2[1])
#             else:
#                 _new_p1 = (_p1[0], _p2[1])
#                 _new_p2 = (_p1[1], _p2[0])
#             return _new_p1, _new_p2
#
#         self.es_pairs_last = cp.copy(self.es_pairs)
#         random.shuffle(self.es_pairs)
#         p1 = self.es_pairs.pop()
#         p2 = self.es_pairs.pop()
#         new_p1, new_p2 = new_pair(p1, p2)
#         while not self.rule(new_p1, new_p2):
#             self.es_pairs.insert(p2, 0)
#             p2 = self.es_pairs.pop()
#             new_p1, new_p2 = new_pair(p1, p2)
#         self.es_pairs.append(new_p1)
#         self.es_pairs.append(new_p2)
#         return self.energy_dif(p1, p2, new_p1, new_p2, T)
#
#     def back(self):
#         self.es_pairs = cp.copy(self.es_pairs_last)
#
#     def total_energy(self):
#         te = 0
#         for pair in self.es_pairs:
#             te += self.pair_energy(pair, 298)
#         return te
#
#     def create_soft(self, pair):
#         if self.chain_creator is None:
#             self.chain_creator = f.ChainCreator(15, 5.13)
#         es1 = self.Atoms[pair[0]]
#         es2 = self.Atoms[pair[1]]
#         r1 = es1.coordinate
#         r2 = es2.coordinate
#         r_2, r2 = self.min_r2(r1, r2)
#         co_list = self.chain_creator.create_chain(r1, r2, 13)
#         atom_id = len(self.Atoms) + 1
#         bond_id = len(self.Bonds) + 1
#         for i, co in enumerate(co_list):
#             if i == 0:
#                 self.append_element(bsc.Atom(atom_id, "TO(2)", co % self.lattice_parameter))
#                 self.append_element(bsc.Bond(bond_id, pair[0], atom_id, _type=("Es", "TO(2)")))
#             else:
#                 self.append_element(bsc.Atom(atom_id, "TO(2)", co % self.lattice_parameter))
#                 self.append_element(bsc.Bond(bond_id, atom_id - 1, atom_id, _type=("TO(2)", "TO(2)")))
#             atom_id += 1
#             bond_id += 1
#         self.append_element(bsc.Bond(bond_id, pair[1], atom_id - 1, _type=("TO(2)", "TO(2)")))
#
#     def create_data_file(self, path):
#         # TODO mole id need fix
#         atom3d = cp.deepcopy(self)
#         Es_ids = set()
#         for a_id, atom in self.Atoms.items():
#             if atom.type == "Es":
#                 Es_ids.add(a_id)
#         # bond_id = max(atom3d.Bonds) + 1
#         for i, pair in enumerate(self.es_pairs):
#             # atom3d.append_element(bsc.Bond(bond_id + i, pair[0], pair[1], _type=("Es", "Es")))
#             atom3d.create_soft(pair)
#             Es_ids.remove(pair[0])
#             Es_ids.remove(pair[1])
#         Es_ids = list(Es_ids)
#         random.shuffle(Es_ids)
#
#         LmpParser.create_data_file(atom3d, path, q=False, improper=False)


class ChainModel:
    def energy(self, n, r2, T):
        return 0

    def probability(self, n, r2):
        return 0


class FreelyJointChain(ChainModel):
    def __init__(self, b):
        self.b = b

    def energy(self, n, r2, T):
        return T * 1.5 * r2 / n / self.b / self.b

    def probability(self, n, r2):
        p = np.power(3 / 2 / np.pi / n / self.b / self.b, 1.5) * np.exp(-3 / 2 / n / self.b / self.b * r2)
        return p


class TOChain(ChainModel):
    def __init__(self):
        import pickle as pkl
        # with open("/home/centos/Projects/CGMD/analysis/overall/data/EndToEnd_param", "rb") as file:
        #     self.param = pkl.load(file)
        with open("/home/centos/Projects/CGMD/analysis/overall/data/EndToEnd_kernel.pkl", "rb") as file:
            self.param = pkl.load(file)
        self.log_dens = dict()
        self.x = np.arange(0, 100, 0.01)
        fp = "/home/centos/Projects/CGMD/analysis/overall/data/kernels.pkl"
        import os
        if os.path.exists(fp):
            with open("/home/centos/Projects/CGMD/analysis/overall/data/kernels.pkl", "rb") as file:
                self.log_dens = pkl.load(file)
        else:
            for n in self.param:
                print(n)
                self.log_dens[n] = self.param[n].score_samples(self.x[:, np.newaxis])
            with open("/home/centos/Projects/CGMD/analysis/overall/data/kernels.pkl", "wb") as file:
                pkl.dump(self.log_dens, file)
        print("WARNING: box length should larger than %.2f" % (self.x[self.log_dens[14].argmax()] * 2))

    def energy(self, n, r2, T):
        r = np.sqrt(r2)
        r_id = int(r / 0.01)
        if r_id >= 10000 - 1:
            log_p = self.log_dens[n][-1]
        else:
            dif = (r - 0.01 * r_id) / 0.01
            log_p = self.log_dens[n][r_id] * (1 - dif) + self.log_dens[n][r_id + 1] * dif
        # return T * (-np.log(k / np.sqrt(2 * np.pi * s2)) + np.square(np.sqrt(r2) - mu) / s2)
        return - T * log_p

    def probability(self, n, r2):
        # k, mu, s2 = self.param[n]
        # return f.single_gaussian(np.sqrt(r2), k, mu, s2)
        r = np.sqrt(r2)
        r_id = int(r / 0.01)
        if r_id >= 9999:
            log_p = self.log_dens[n][-1]
        else:
            dif = (r - 0.01 * r_id) / 0.01
            log_p = self.log_dens[n][r_id] * (1 - dif) + self.log_dens[n][r_id + 1] * dif
        return np.exp(log_p)


class ChainCreator:
    def __init__(self, chain_model: ChainModel, b, exclude_func=None):
        self.chain_model = chain_model
        self.b = b
        self.exclude_func = exclude_func

    def new_co(self, last_co):
        # n_co = last_co + self.b * cso.random_unit() * np.random.uniform(0.8, 1.2)
        n_co = last_co + self.b * cso.random_unit() * np.random.uniform(0.9, 1.1)
        if self.exclude_func is not None:
            while self.exclude_func(n_co):
                # n_co = last_co + self.b * cso.random_unit() * np.random.uniform(0.8, 1.2)
                n_co = last_co + self.b * cso.random_unit() * np.random.uniform(0.9, 1.1)
        return n_co

    def cal_p(self, co, s, N, r1, r2):
        dr1_2 = np.inner(r1 - co, r1 - co)
        p1 = self.chain_model.probability(s, dr1_2)
        dr2_2 = np.inner(r2 - co, r2 - co)
        p2 = self.chain_model.probability(N - s, dr2_2)
        p = p1 * p2
        return p

    def cal_p_single(self, co, s, N, r1):
        dr1_2 = np.inner(r1 - co, r1 - co)
        p1 = self.chain_model.probability(s, dr1_2)
        return p1

    def evaluate(self, l_co, s, N, r1, r2):
        p_max = -1
        co = None
        for i in range(100):
            new_co = self.new_co(l_co)
            p = self.cal_p(new_co, s, N, r1, r2)
            # print(p)
            if p > p_max:
                p_max = p
                co = new_co
        return p_max, co

    def evaluate_single(self, l_co, s, N, r1):
        p_max = -1
        co = None
        for i in range(100):
            new_co = self.new_co(l_co)
            p = self.cal_p_single(new_co, s, N, r1)
            # print(p)
            if p > p_max:
                p_max = p
                co = new_co
        return p_max, co

    def gen_new_co(self, l_co, s, N, p_max, r1, r2):
        new_co = self.new_co(l_co)
        p = self.cal_p(new_co, s, N, r1, r2)
        i = 0
        while p < p_max * np.random.rand():
            new_co = self.new_co(l_co)
            p = self.cal_p(new_co, s, N, r1, r2)
            i += 1
            if i > 1000:
                return [-1]
        return new_co

    def gen_new_co_single(self, l_co, s, N, p_max, r1):
        new_co = self.new_co(l_co)
        p = self.cal_p_single(new_co, s, N, r1)
        i = 0
        while p < p_max * np.random.rand():
            new_co = self.new_co(l_co)
            p = self.cal_p_single(new_co, s, N, r1)
            i += 1
            if i > 10000:
                return [-1]
        return new_co

    def create_chain(self, r1, r2, bead_num):
        r1 = np.asarray(r1)
        r2 = np.asarray(r2)
        l_co = r1.copy()
        ll_co = None
        co_list = list()
        s = 0
        while s < bead_num:
            # for s in range(1, bead_num + 1):
            s += 1
            p_max, new_co = self.evaluate(l_co, s, bead_num + 1, r1, r2)
            # print(p_max)
            new_co = self.gen_new_co(l_co, s, bead_num + 1, p_max * 10, r1, r2)
            if len(new_co) == 1:
                print("recal chain")
                s = 0
                co_list = list()
                l_co = r1.copy()
                continue
            # print(new_co)
            co_list.append(new_co)
            ll_co = l_co
            l_co = new_co
        return co_list

    def create_one_end_chain(self, r1, bead_num):
        r1 = np.asarray(r1)
        co_list = list()
        l_co = r1
        s = 0
        while s < bead_num:
            # for s in range(1, bead_num + 1):
            s += 1
            p_max, new_co = self.evaluate_single(l_co, s, bead_num + 1, r1)
            # print(p_max)
            new_co = self.gen_new_co_single(l_co, s, bead_num + 1, p_max * 10, r1)
            if len(new_co) == 1:
                print("recal chain")
                s = 0
                co_list = list()
                continue
            # print(new_co)
            co_list.append(new_co)
            l_co = new_co
        return co_list


class HardBox(csc.Atom3D):
    def __init__(self, chain_model: ChainModel, chain_creator, seg_judge=None, intra_p=1.0):
        super().__init__()
        self.chain_creator = chain_creator
        # self.super_to_mole = dict()
        # self.mole_to_super = dict()
        self.super_to_pair = dict()
        self.super_to_pair_last = None
        self.repeat_num = 0
        self.super_num = 0
        self.chain_model = chain_model
        self.mole_hard_dict = dict()
        self.hard_density = None
        self.grid = 1.0
        self.extend = 30.
        self.seg_judge = seg_judge
        self.intra_p = intra_p

    def random_es_pair(self):
        Es_ids = list()
        mole_es = dict()
        for a_id, atom in self.Atoms.items():
            if atom.type == "Es":
                Es_ids.append(a_id)
                if atom.mole_id in mole_es:
                    mole_es[atom.mole_id].append(atom.id)
                else:
                    mole_es[atom.mole_id] = [atom.id]
        random.shuffle(Es_ids)
        for i in range(self.super_num):
            last_e2 = 0
            for j in range(self.repeat_num):
                e1_id = Es_ids.pop()
                mole_es[self.Atoms[e1_id].mole_id].remove(e1_id)
                e2_id = mole_es[self.Atoms[e1_id].mole_id].pop()
                Es_ids.remove(e2_id)
                if last_e2 > 0:
                    # self.es_pairs.append((e1_id, last_e2))
                    self.super_to_pair[i + 1].append((last_e2, e1_id))
                else:
                    self.super_to_pair[i + 1] = [(0, e1_id)]
                last_e2 = e2_id
            self.super_to_pair[i + 1].append((last_e2, 0))

    def initialize(self, chain_num=204, repeat_num=4):
        self.repeat_num = repeat_num
        self.super_num = chain_num // repeat_num
        fp1 = "../data/20200420_lr=0.05_1blk_40_cg_ex_hard/record/"
        r1 = ibi.Record(fp1, restart=True)
        param1 = r1.get_param(22, 1)
        pattern = ["TO(2)"] * 6 + ["Es", "Ph", "U", "Ph", "Me", "Ph", "U", "Ph", "Es"] + ["TO(2)"] * 7
        pattern[0] = "TO(1)"
        pattern[-1] = "TO(1)"
        self.create_from_pattern_fix2(pattern, chain_num, param=param1, density=1., center=10, soft=False)
        self.random_es_pair()

    def initialize_atom3d(self, atom3d, repeat_num=4):
        atom3d.convert_to_in_cell()
        self.lattice_parameter = atom3d.lattice_parameter.copy()
        self.box_h = atom3d.box_h.copy()
        self.box_l = atom3d.box_l.copy()
        self.repeat_num = repeat_num
        self.super_num = len(atom3d.mole_dict) // repeat_num
        for atom in atom3d.Atoms.values():
            if atom.type[:2] != "TO":
                self.append_element(cp.deepcopy(atom))
        for bond in atom3d.Bonds.values():
            if "TO(1)" in bond.type or "TO(2)" in bond.type:
                continue
            self.append_element(cp.deepcopy(bond))
        self.sort_atom_bond_id()
        self.random_es_pair()
        for a_id, atom in self.Atoms.items():
            if atom.type in self.hard_type or atom.type == "Es":
                if atom.mole_id in self.mole_hard_dict:
                    self.mole_hard_dict[atom.mole_id].append(a_id)
                else:
                    self.mole_hard_dict[atom.mole_id] = [a_id]
        # self.hard_density = self.density_matrix(atom3d.get_fake_copy())
        # self.hard_density /= np.mean(self.hard_density)

    def density_matrix(self, atom3d):
        atom3d.extend_self(extend_range=self.extend)
        n_xa, n_ya, n_za = ((self.lattice_parameter + self.extend * 2) // self.grid + 1).astype(int)
        H_sum = np.zeros((n_xa, n_ya, n_za))
        for atom in atom3d.Atoms.values():
            # exclude soft for better contrast
            atom_type = atom.type
            if atom_type[:2] == "TO" or atom_type[:2] == "Es":
                continue
            x, y, z = ((atom.coordinate - self.box_l + self.extend) / self.grid).astype(int)
            # print(atom_type)
            # print(x, y, z)
            H_sum[x, y, z] += atom3d.type_mass_dict[atom_type]
        # print(H_sum)
        return H_sum

    # def min_r2(self, r1, r2):
    #     dr = r2 - r1
    #     dr = dr % self.lattice_parameter
    #     for i in range(3):
    #         if dr[i] <= 0:
    #             dr[i] = dr[i] + self.lattice_parameter[i] if dr[i] + self.lattice_parameter[i] < -dr[i] else dr[i]
    #         else:
    #             dr[i] = dr[i] - self.lattice_parameter[i] if self.lattice_parameter[i] - dr[i] < dr[i] else dr[i]
    #     dr2 = np.inner(dr, dr)
    #     return dr2, r1 + dr

    def min_r2(self, r1, r2):
        dr = r2 - r1
        dr = dr % self.lattice_parameter
        judge_p = np.where(dr > self.lattice_parameter / 2, 1, 0)
        judge_n = np.where(dr < -self.lattice_parameter / 2, -1, 0)
        judge = judge_n + judge_p
        dr -= judge * self.lattice_parameter
        # if np.inner(dr, dr) > np.inner(self.lattice_parameter/2, self.lattice_parameter/2):
        #     print(dr)

        return np.inner(dr, dr), r1 + dr

    def pair_energy(self, pair, T, N=14, b=5.13):
        if pair[0] == 0 or pair[1] == 0:
            return 0.
        r_2, new_r2 = self.min_r2(self.Atoms[pair[0]].coordinate.copy(), self.Atoms[pair[1]].coordinate.copy())
        E = self.chain_model.energy(N, r_2, T)
        if self.seg_judge is not None:
            image = (new_r2 - self.Atoms[pair[1]].coordinate) // self.lattice_parameter
            judge = self.seg_judge.judge(self.Atoms[pair[0]], self.Atoms[pair[1]], image)
            # print(judge)
            if self.intra_p > 0.9999999 and not judge:
                E += 50000
            elif self.intra_p < 0.0000001 and judge:
                E += 50000
        '''
        r1 = self.Atoms[pair[0]].coordinate.copy()
        gx1, gy1, gz1 = ((r1 - self.box_l + self.extend) / self.grid).astype(int)
        gx2, gy2, gz2 = ((new_r2 - self.box_l + self.extend) / self.grid).astype(int)
        mean_density = np.mean(
            self.hard_density[min(gx1, gx2):max(gx1, gx2)+1, min(gy1, gy2):max(gy1, gy2)+1, min(gz1, gz2):max(gz1, gz2)+1])
        if np.isnan(mean_density):
            mean_density = 1.

        return self.chain_model.energy(N, r_2, T) + (mean_density - 1.) * 1000
        '''
        return E

    def energy_dif(self, op1, op2, np1, np2, T):
        return self.pair_energy(np1, T) + self.pair_energy(np2, T) - self.pair_energy(op1, T) - self.pair_energy(op2, T)

    def move(self, T):
        def new_pair(_p1, _p2):
            _new_p1 = (_p1[0], _p2[1])
            _new_p2 = (_p2[0], _p1[1])
            return _new_p1, _new_p2

        self.super_to_pair_last = cp.copy(self.super_to_pair)
        s1_ran = np.random.randint(1, self.super_num + 1)
        s2_ran = np.random.randint(1, self.super_num + 1)
        p_ran = np.random.randint(1, self.repeat_num)
        p1 = self.super_to_pair[s1_ran][p_ran]
        p2 = self.super_to_pair[s2_ran][p_ran]
        new_p1, new_p2 = new_pair(p1, p2)
        new_s1 = self.super_to_pair[s1_ran][:p_ran] + [new_p1] + self.super_to_pair[s2_ran][p_ran + 1:]
        new_s2 = self.super_to_pair[s2_ran][:p_ran] + [new_p2] + self.super_to_pair[s1_ran][p_ran + 1:]
        self.super_to_pair[s1_ran] = new_s1
        self.super_to_pair[s2_ran] = new_s2
        return self.energy_dif(p1, p2, new_p1, new_p2, T)

    def back(self):
        self.super_to_pair = cp.copy(self.super_to_pair_last)

    def total_energy(self):
        te = 0
        for super_pairs in self.super_to_pair.values():
            for pair in super_pairs:
                te += self.pair_energy(pair, 298)
        return te

    def hard_seg_convert_to_default(self, mole_id):
        for a_id in self.mole_hard_dict[mole_id]:
            atom = self.Atoms[a_id]
            atom.set_default_coordinate(self.lattice_parameter)

    def arrange_hard_seg(self, mole_id, es_id, rn):
        dr = rn - self.Atoms[es_id].coordinate
        for a_id in self.mole_hard_dict[mole_id]:
            atom = self.Atoms[a_id]
            atom.coordinate += dr
        self.in_cell = False
        self.default = False

    def create_soft(self, pair):
        atom_id = len(self.Atoms) + 1
        bond_id = len(self.Bonds) + 1

        to_num = 0
        es_label = 0
        if pair[0] == 0:
            to_num = 6
            es_label = 1
            es2 = self.Atoms[pair[1]]
            self.hard_seg_convert_to_default(es2.mole_id)
        elif pair[1] == 0:
            to_num = 7
            es_label = 0

        if to_num > 0:
            _co = self.Atoms[pair[es_label]].coordinate.copy()
            mole_id = self.Atoms[pair[es_label]].mole_id
            co_list = self.chain_creator.create_one_end_chain(_co, to_num)
            for i, _co in enumerate(co_list):
                if i == 0:
                    self.append_element(
                        bsc.Atom(atom_id, "TO(2)", _co.copy(), mole_id=cp.copy(mole_id)))
                    self.append_element(bsc.Bond(bond_id, pair[es_label], atom_id, _type=("Es", "TO(2)")))
                elif i == to_num - 1:
                    self.append_element(
                        bsc.Atom(atom_id, "TO(1)", _co.copy(), mole_id=cp.copy(mole_id)))
                    self.append_element(bsc.Bond(bond_id, atom_id, atom_id - 1, _type=("TO(1)", "TO(2)")))
                else:
                    self.append_element(
                        bsc.Atom(atom_id, "TO(2)", _co.copy(), mole_id=cp.copy(mole_id)))
                    self.append_element(bsc.Bond(bond_id, atom_id, atom_id - 1, _type=("TO(2)", "TO(2)")))
                atom_id += 1
                bond_id += 1
        else:
            es1 = self.Atoms[pair[0]]
            es2 = self.Atoms[pair[1]]
            # print(es1.coordinate, es2.coordinate)
            self.hard_seg_convert_to_default(es2.mole_id)
            # print(es1.coordinate, es2.coordinate)
            mole_id = es1.mole_id
            r1 = es1.coordinate.copy()
            r2 = es2.coordinate.copy()
            r_2, r2 = self.min_r2(r1, r2)
            self.arrange_hard_seg(es2.mole_id, es2.id, r2)
            co_list = self.chain_creator.create_chain(r1, r2, 13)
            # print(r1, r2)
            while cso.calculate_distance(co_list[-1], r2) > 1.5 * self.chain_creator.b:
                co_list = self.chain_creator.create_chain(r1, r2, 13)
            for i, co in enumerate(co_list):
                if i == 0:
                    self.append_element(
                        bsc.Atom(atom_id, "TO(2)", co.copy(), mole_id=cp.copy(mole_id)))
                    self.append_element(bsc.Bond(bond_id, pair[0], atom_id, _type=("Es", "TO(2)")))
                else:
                    self.append_element(
                        bsc.Atom(atom_id, "TO(2)", co.copy(), mole_id=cp.copy(mole_id)))
                    self.append_element(bsc.Bond(bond_id, atom_id - 1, atom_id, _type=("TO(2)", "TO(2)")))
                atom_id += 1
                bond_id += 1
            self.append_element(bsc.Bond(bond_id, pair[1], atom_id - 1, _type=("Es", "TO(2)")))

    def create_data_file(self, path):
        atom3d = cp.deepcopy(self)
        mole_super = dict()
        mole_dict = dict()
        print("create soft")
        len_block = len(self.super_to_pair.items()) // 20
        for super_id, super_pairs in self.super_to_pair.items():
            if super_id % len_block == 0:
                print("%d%%" % int(super_id // len_block * 5))
            for pair in super_pairs:
                # print(pair)
                atom3d.create_soft(pair)
                if pair[0] > 0:
                    mole_super[atom3d.Atoms[pair[0]].mole_id] = super_id
        print()
        # print(mole_super)
        for a_id in atom3d.Atoms:
            super_id = mole_super[atom3d.Atoms[a_id].mole_id]
            atom3d.Atoms[a_id].mole_id = super_id
            if super_id not in mole_dict:
                mole_dict[super_id] = [a_id]
            else:
                mole_dict[super_id].append(a_id)
        atom3d.mole_dict = mole_dict
        print("sort atom bond id")
        atom3d.sort_atom_bond_id()
        print("cal all joints")
        atom3d.cal_all_joints()
        if "c1" in atom3d.Atom_type_dict:
            LmpParser.create_data_file(atom3d, path, q=False, improper=True)
        else:
            LmpParser.create_data_file(atom3d, path, q=False, improper=False)
        return atom3d

    def out(self):
        if self.seg_judge is not None:
            sum_p = 0
            intra_p = 0
            for super_pairs in self.super_to_pair.values():
                for pair in super_pairs:
                    if pair[0] == 0 or pair[1] == 0:
                        continue
                    r_2, new_r2 = self.min_r2(self.Atoms[pair[0]].coordinate.copy(), self.Atoms[pair[1]].coordinate.copy())
                    image = (new_r2 - self.Atoms[pair[1]].coordinate) // self.lattice_parameter
                    judge = self.seg_judge.judge(self.Atoms[pair[0]], self.Atoms[pair[1]], image)
                    sum_p += 1
                    if judge:
                        intra_p += 1
            print("intra rate: %.4f" % (intra_p / sum_p))


class MonteCarlo:

    def __init__(self, init_state):
        self.state = init_state

    def simulation(self, T, step):
        for i in range(step):
            dE = self.state.move(T)
            if dE > 0:
                p = np.exp(-1 / T * dE)
                if p <= np.random.rand():
                    self.state.back()
            if i % 10000 == 0:
                print("step %d\t, total energy: %.2f" % (i, self.state.total_energy()))
                self.state.out()


if __name__ == '__main__':
    # atom3d = LmpParser.load_data_file("/home/centos/work/SH6S/cgu_1blk_2000chain/init/strain.data")
    # atom3d.renew_coordinate_file("/home/centos/work/SH6S/cgu_1blk_2000chain/1E-6/strain.lammpstrj.s1.0")

    # atom3d = LmpParser.load_data_file("/home/centos/work/1blk_50/cg_1blk_2000chain/init/strain.data")
    # atom3d.renew_coordinate_file("/home/centos/work/1blk_50/cg_1blk_2000chain/init/record.lammpstrj.10300000")

    # atom3d = LmpParser.load_data_file("/home/centos/work/1blk_50/cgu_1blk_3200chain/init/1blk_3200.data")
    # atom3d.renew_coordinate_file(
    #     "/home/centos/work/1blk_50/cgu_1blk_3200chain/init/trj/1blk_3200.lammpstrj.010000000")

    # atom3d = LmpParser.load_data_file("/home/centos/work/1blk_50/cgu_1blk_160chain/init/1blk_160.data")
    # atom3d.renew_coordinate_file(
    #     "/home/centos/work/1blk_50/cgu_1blk_160chain/init/trj/1blk_160.lammpstrj.010000000")
    # chains = 160
    # chain_model = TOChain()
    # chain_creator = ChainCreator(chain_model, 5.13)
    #
    # print(chain_model.energy(1, 25, 298))
    pass
