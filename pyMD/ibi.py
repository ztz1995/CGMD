from pyMD import collective_structure_class as csc
from pyMD import file_parser as fps
from pyMD import functions as f
import numpy as np
from multiprocessing import Pool, Process
import time
import os
import shutil
import pandas as pd
from collections import OrderedDict
import json
import matplotlib.pyplot as plt
import copy as cp
import os


class Record:

    def __init__(self, path, restart=False):
        self.path = path
        self.restart = restart
        if not restart:
            if not os.path.exists(path):
                os.makedirs(path)
            self.current_state = [0] + [False] * 3  # [iter_num, pre_adjust, cal_trajectory, cal_g(r)]
            json.dump(self.current_state, open(self.path + "current_state.json", "w+"))
        else:
            self.current_state = json.load(open(path + "current_state.json"))
        # self.log = CalLog(name, path, restart)
        # self.param_log = CalLog("param.log", path, restart)
        # self.error_log = CalLog("error.log", path, restart)
        # self.run_log = CalLog("run.log", path, restart)

    @staticmethod
    def dump_dict_json(obj, file_path):
        assert isinstance(obj, dict)
        json_dict = dict()
        for key, value in obj.items():
            json_dict[key] = dict()
            for k, v in value.items():
                # print(k)
                # print(v)
                if isinstance(v, (list, np.ndarray, tuple)):
                    json_dict[key]["_".join(k)] = list(v)
                else:
                    json_dict[key]["_".join(k)] = round(v, 3)

        with open(file_path, "w") as file:
            json.dump(json_dict, file)
        pass

    @staticmethod
    def load_dict_json(file_path):
        with open(file_path, "r") as file:
            json_dict = json.load(file)
        obj_dict = dict()
        for key, value in json_dict.items():
            obj_dict[key] = dict()
            for k, v in value.items():
                obj_dict[key][tuple(k.split("_"))] = np.array(v)
                # obj_dict[key][tuple(k.split("_"))] = v if key != "non_bond" else np.array(v)
        return obj_dict

    def renew_state(self, new_state, parameters=None, error=None):
        if new_state == -1:
            # self.param_log.renew_iter(self.current_state[0])
            # self.param_log.renew_dict_of_dict(parameters, "param")
            self.dump_dict_json(parameters, self.path + "param_iter_%02d_0.json" % self.current_state[0])
        elif new_state == 0:
            assert self.current_state[1:] == [True] * 3
            self.current_state = [self.current_state[0] + 1] + [False] * 3
            # self.error_log.renew_dict_of_dict(error, "error")
            # self.error_log.renew_iter(self.current_state[0])
            # if self.current_state[0] != 0:
            # self.param_log.renew_iter(self.current_state[0])
            # self.param_log.renew_dict_of_dict(parameters, "param")
            self.dump_dict_json(error, self.path + "error_iter_%02d.json" % (self.current_state[0] - 1))
            self.dump_dict_json(parameters, self.path + "param_iter_%02d_0.json" % self.current_state[0])
        else:
            assert self.current_state[1:new_state] == [True] * (new_state - 1) and self.current_state[new_state:] == [
                False] * (4 - new_state)
            self.current_state[new_state] = not self.current_state[new_state]
            if new_state == 1:
                # self.param_log.renew_dict_of_dict(parameters, "param")
                self.dump_dict_json(parameters, self.path + "param_iter_%02d_1.json" % self.current_state[0])
        json.dump(self.current_state, open(self.path + "current_state.json", "w+"))
        self.restart = False
        pass

    def get_param(self, iter_num=None, pre=1):
        if iter_num is None:
            if self.current_state[1]:
                return self.load_dict_json(self.path + "param_iter_%02d_1.json" % self.current_state[0])
            else:
                return self.load_dict_json(self.path + "param_iter_%02d_0.json" % self.current_state[0])
        else:
            return self.load_dict_json(self.path + "param_iter_%02d_%d.json" % (iter_num, pre))

    def get_init_param(self):
        return self.load_dict_json(self.path + "param_iter_00_1.json")


class IBI:
    para_list = ["non_bond", 'Bond', 'Angle', 'Dihedral']
    # para_list = ['Bond', 'Angle', 'Dihedral']
    cal_points = {"non_bond": np.arange(1800) * 0.01 + 0.01,
                  # para_list[1]: np.arange(2000) * 0.005 + 0.005,
                  "Bond": np.arange(1000) * 0.01 + 0.01,
                  "Angle": np.arange(1800) * 0.1 + 0.05,
                  "Dihedral": np.arange(1800) * 0.1 + 0.05}  # np.arange(-1800, 1800, 1) * 0.1

    default_style = {"non_bond": "table", 'Bond': "fene/expand", 'Angle': "harmonic", 'Dihedral': "multi/harmonic",
                     "Improper": "class2"}

    style_func = {"harmonic": f.harmonic, "multi/harmonic": f.multi_harmonic, "fene/expand": f.fene_expand,
                  "harmonic_d": f.harmonic_d}

    fit_func = {"non_bond": f.fit_non_bond, 'Bond': f.fit_bond, 'Angle': f.fit_angle, 'Dihedral': f.fit_dihedral}

    # aa_type = {"h1", "n1", "o1", "c3"}
    # new type
    aa_type = {"h1", "n1", "o1", "c1"}

    non_bond_cutoff = 14.0000000001
    cut_idx = 1400

    # # cgu
    # class2_type = {("h1", "n1"): [1.0100, 462.7500, -1053.6355, 1545.7701],  # 13
    #                ("c1", "n1"): [1.3880, 440.6783, -828.3871, 1423.2587],  # 14
    #                ("c1", "o1"): [1.2160, 823.7948, -1878.8288, 2303.4950],  # 15
    #                # angle
    #                ("c1", "n1", "h1"): {" ": [122.9480, 40.4820, -16.2009, 8.3271], "bb": [8.6253, 1.3880, 1.0100],
    #                                     "ba": [34.8312, 15.0778, 1.3880, 1.0100]},  # 21
    #                ("n1", "c1", "o1"): {" ": [125.5320, 101.8765, -41.8101, 7.7222], "bb": [115.4645, 1.3880, 1.2160],
    #                                     "ba": [32.8758, 46.1093, 1.3880, 1.2160]},  # 24
    #                ("n1", "c1", "n1"): {" ": [120.5292, 100.0857, -36.7315, 24.2608], "bb": [84.5263, 1.3880, 1.3880],
    #                                     "ba": [49.0875, 49.0875, 1.3880, 1.3880]},  # 25
    #                # dihedral
    #                ("h1", "n1", "c1", "o1"): {" ": [0.0005, 0.0000, 2.0516, 0.0000, -0.0013, 0.0000],
    #                                           "mbt": [0.0009, 4.4689, -0.0008, 1.3880],
    #                                           "ebt": [-0.0006, -0.0016, -0.0004, 0.0013, 0.0005, 0.0008, 1.0100,
    #                                                   1.2160],
    #                                           "at": [-0.001, -0.001, -0.0008, 0.0015, 0.0015, 0.0015, 122.948,
    #                                                  125.5320],
    #                                           "aat": [-0.0003, 122.9480, 125.5320],
    #                                           "bb13": [-0.0000, 1.0100, 1.2160]},  # 31
    #                ("h1", "n1", "c1", "n1"): {" ": [-1.0630, 0.0000, 1.5631, 0.0000, 0.0000, 0.0000],
    #                                           "mbt": [0.0001, 6.3289, -0.0007, 1.3880],
    #                                           "ebt": [-0.0006, -0.0015, -0.0019, -0.0002, 0.0008, 0.0004, 1.0100,
    #                                                   1.3880],
    #                                           "at": [0.0005, -0.0002, -0.0002, -0.0008, -0.0004, -0.0002, 122.9480,
    #                                                  120.5292],
    #                                           "aat": [0.0006, 122.9480, 120.5292],
    #                                           "bb13": [0.0021, 1.0100, 1.3880]},  # 30
    #                # improper
    #                ("n1", "c1", "n1", "o1"): {" ": [59.37400, 0.00000],  # 16
    #                                           "aa": [-0.0001, -0.0001, -0.0001, 120.5292, 125.5320, 125.5320]},
    #                ("c1", "n1", "Ph", "h1"): {" ": [4.4181, 0.0000],  # 15
    #                                           "aa": [0.0, 0.0, 0.0, 120.07, 122.9480, 116.3230]}}

    # fg
    class2_type = {
        "Bond": {
            ("h1", "n1"): [1.0100, 462.7500, -1053.6355, 1545.7701],  # 13
            ("c1", "n1"): [1.3880, 440.6783, -828.3871, 1423.2587],  # 14
            ("c1", "o1"): [1.2160, 823.7948, -1878.8288, 2303.4950],  # 15
            ("c2", "n1"): [1.3950, 344.0452, -652.1377, 1022.2271],  # 12
            ("c2", "c2"): [1.4170, 470.8361, -627.6245, 1327.6166],  # 10
            ('c2', 'c3'): [1.5010, 321.9021, -521.8355, 572.1488],  # 16
        },
        "Angle": {
            ("c1", "n1", "h1"): {" ": [122.9480, 40.4820, -16.2009, 8.3271], "bb": [8.6253, 1.3880, 1.0100],
                                 "ba": [34.8312, 15.0778, 1.3880, 1.0100]},  # 21
            ("n1", "c1", "o1"): {" ": [125.5320, 101.8765, -41.8101, 7.7222], "bb": [115.4645, 1.3880, 1.2160],
                                 "ba": [32.8758, 46.1093, 1.3880, 1.2160]},  # 24
            ("n1", "c1", "n1"): {" ": [120.5292, 100.0857, -36.7315, 24.2608], "bb": [84.5263, 1.3880, 1.3880],
                                 "ba": [49.0875, 49.0875, 1.3880, 1.3880]},  # 25
            ('c2', 'c2', 'c2'): {" ": [118.9000, 61.0226, -34.9904, 0.0000],
                                 "bb": [68.2856, 1.4170, 1.4170],
                                 "ba": [28.8708, 28.8708, 1.4170, 1.4170]},  # 17
            ('c2', 'c2', 'n1'): {" ": [120.7640, 73.2738, -27.4044, 13.3945],
                                 "bb": [37.8749, 1.4170, 1.3950],
                                 "ba": [35.8865, 53.6977, 1.4170, 1.3950]},  # 20
            ('c2', 'c2', 'c3'): {" ": [120.0500, 44.7148, -22.7330, 0.0000],
                                 "bb": [12.0676, 1.4170, 1.5010],
                                 "ba": [31.0771, 47.0579, 1.4170, 1.5010]},  # 26
            ('c1', 'n1', 'c2'): {" ": [120.0700, 47.1131, -32.5599, 13.1257],
                                 "bb": [41.4233, 1.3880, 1.3950],
                                 "ba": [34.7791, 24.3705, 1.3880, 1.3950]},  # 23
            ('c2', 'n1', 'h1'): {" ": [116.3230, 18.3123, -7.8322, 5.3289],
                                 "bb": [8.2930, 1.3950, 1.0100],
                                 "ba": [10.4568, 12.8217, 1.3950, 1.0100]},  # 22
            # ('c2', 'c3', 'c2'): {" ": [111.0000, 44.3234, -9.4453, 0.0000],
            #                      "bb": [-0.0000, 1.5010, 1.5010],
            #                      "ba": [0.0000, 0.0000, 1.5010, 1.5010]},  # 28
            # ('Es', 'c2', 'c2'): {"": [116.0640, 71.2598, -15.8268, 2.0523],
            #                      "bb": [37.8749, 1.4890, 1.4170],
            #                      "ba": [23.6977, 45.8865, 1.4890, 1.4170]}  # 18
        },
        "Dihedral": {
            ("h1", "n1", "c1", "o1"): {" ": [0.0005, 0.0000, 2.0516, 0.0000, -0.0013, 0.0000],
                                       "mbt": [0.0009, 4.4689, -0.0008, 1.3880],
                                       "ebt": [-0.0006, -0.0016, -0.0004, 0.0013, 0.0005, 0.0008, 1.0100,
                                               1.2160],
                                       "at": [-0.001, -0.001, -0.0008, 0.0015, 0.0015, 0.0015, 122.948,
                                              125.5320],
                                       "aat": [-0.0003, 122.9480, 125.5320],
                                       "bb13": [-0.0000, 1.0100, 1.2160]},  # 31
            ("h1", "n1", "c1", "n1"): {" ": [-1.0630, 0.0000, 1.5631, 0.0000, 0.0000, 0.0000],
                                       "mbt": [0.0001, 6.3289, -0.0007, 1.3880],
                                       "ebt": [-0.0006, -0.0015, -0.0019, -0.0002, 0.0008, 0.0004, 1.0100,
                                               1.3880],
                                       "at": [0.0005, -0.0002, -0.0002, -0.0008, -0.0004, -0.0002, 122.9480,
                                              120.5292],
                                       "aat": [0.0006, 122.9480, 120.5292],
                                       "bb13": [0.0021, 1.0100, 1.3880]},  # 30
            ('c2', 'n1', 'c1', 'o1'): {" ": [-0.0003, 0.0000, 2.0522, 0.0000, 0.0003, 0.0000],
                                       "mbt": [-0.0024, 4.4645, -0.0045, 1.3880],
                                       "ebt": [-0.0002, -0.0010, -0.0005, -0.0006, -0.0011, -0.0020, 1.3950, 1.2160],
                                       "at": [0.0011, 0.0008, 0.0008, 0.0004, 0.0002, -0.0002, 120.0700, 125.5320],
                                       "aat": [0.0008, 120.0700, 125.5320],
                                       "bb13": [0.0009, 1.3950, 1.2160]},  # 33
            # # ('c1', 'n1', 'c2', 'c2'): {" ": [0.0001, 0.0000, 0.5107, 0.0000, 0.0000, 0.0000],
            # ('c1', 'n1', 'c2', 'c2'): {" ": [0.0001, 0.0000, 0.6007, 0.0000, 0.0000, 0.0000],
            #                            "mbt": [-0.0000, 4.9027, -0.0000, 1.3950],
            #                            "ebt": [-0.0000, -0.0000, -0.0000, 0.0000, 0.0000, 0.0000, 1.3880, 1.4170],
            #                            "at": [-0.0000, -0.0000, -0.0000, -0.0005, -0.0002, -0.0001, 120.0700, 120.7640],
            #                            "aat": [0.0000, 120.0700, 120.7640],
            #                            "bb13": [-0.0000, 1.3880, 1.4170]},  # 29
            # # ('c2', 'c2', 'n1', 'h1'): {" ": [0.0000, 0.0000, 0.5107, 0.0000, 0.0000, 0.0000],
            # ('c2', 'c2', 'n1', 'h1'): {" ": [0.0000, 0.0000, 0.6007, 0.0000, 0.0000, 0.0000],
            #                            "mbt": [-0.0000, 2.4730, 0.0000, 1.3950],
            #                            "ebt": [0.0000, -0.0000, -0.0000, -0.0000, -0.0000, -0.0000, 1.4170, 1.0100],
            #                            "at": [0.0000, 0.0000, 0.0000, -0.0003, -0.0001, -0.0000, 120.7640, 116.3230],
            #                            "aat": [-0.0000, 120.7640, 116.3230],
            #                            "bb13": [0.0000, 1.4170, 1.0100]},  # 28
            ('c2', 'c2', 'c2', 'c2'): {" ": [8.3667, 0.0000, 1.2000, 0.0000, 0.0000, 0.0000],
                                       "mbt": [27.5989, -2.3120, 0.0000, 1.4170],
                                       "ebt": [-0.1185, 6.3204, 0.0000, -0.1185, 6.3204, 0.0000, 1.4170, 1.4170],
                                       "at": [1.9767, 1.0239, -0.0000, 1.9767, 1.0239, -0.0000, 118.9000, 118.9000],
                                       "aat": [-0.0000, 118.9000, 118.9000],
                                       "bb13": [53.0000, 1.4170, 1.4170]},  # 22
            ('c2', 'n1', 'c1', 'n1'): {" ": [-1.0631, 0.0000, 1.5631, 0.0000, -0.0002, 0.0000],
                                       "mbt": [0.0009, 6.3274, -0.0001, 1.3880],
                                       "ebt": [0.0003, -0.0005, 0.0004, 0.0001, -0.0001, 0.0001, 1.3950, 1.3880],
                                       "at": [0.0005, 0.0005, 0.0005, 0.0000, 0.0001, 0.0002, 120.0700, 120.5292],
                                       "aat": [-0.0000, 120.0700, 120.5292],
                                       "bb13": [-0.0018, 1.3950, 1.3880]},  # 32
            ('c2', 'c2', 'c2', 'c3'): {" ": [-0.0000, 0.0000, 4.4072, 0.0000, -0.0000, 0.0000],
                                       "mbt": [-0.0000, 9.1792, -0.0000, 1.4170],
                                       "ebt": [0.0000, -0.6918, 0.0000, -0.0000, 0.2421, -0.0000, 1.4170, 1.5010],
                                       "at": [-0.0000, 3.8987, -0.0000, 0.0000, -4.4683, 0.0000, 118.9000, 120.0500],
                                       "aat": [-14.4097, 118.9000, 120.0500],
                                       "bb13": [2.5085, 1.4170, 1.5010]},  # 35
            ('c2', 'c2', 'c3', 'c2'): {" ": [-0.2802, 0.0000, -0.0678, 0.0000, -0.0122, 0.0000],
                                       "mbt": [-0.0000, -0.0000, 0.0000, 1.5010],
                                       "ebt": [0.0000, 0.0000, 0.0000, -0.0000, -0.0000, -0.0000, 1.4170, 1.5010],
                                       "at": [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 120.0500, 111.0000],
                                       "aat": [0.0000, 120.0500, 111.0000],
                                       "bb13": [0.0000, 1.4170, 1.5010]},  # 37
            ('c2', 'c2', 'c2', 'n1'): {" ": [0.0001, 0.0000, 3.4040, 0.0000, 0.0000, 0.0000],
                                       "mbt": [0.0000, 5.2012, 0.0000, 1.4170],
                                       "ebt": [0.0000, -0.0000, -0.0000, -0.0000, 0.0000, 0.0000, 1.4170, 1.3950],
                                       "at": [-0.0000, -0.0000, -0.0000, -0.0005, -0.0002, -0.0001, 118.9000, 120.7640],
                                       "aat": [-0.0000, 118.9000, 120.7640],
                                       "bb13": [-0.0000, 1.4170, 1.3950]},  # 27
        },
        "Improper": {
            ("n1", "c1", "n1", "o1"): {" ": [59.37400, 0.00000],  # 16
                                       "aa": [-0.0001, -0.0001, -0.0001, 120.5292, 125.5320, 125.5320]},
            ("c1", "n1", "c2", "h1"): {" ": [4.4181, 0.0000],  # 15
                                       "aa": [0.0, 0.0, 0.0, 120.07, 122.9480, 116.3230]},
            ('c2', 'c2', 'c2', 'n1'): {" ": [17.0526, 0.0000],
                                       "aa": [0.0000, 0.0000, 0.0000, 118.9000, 120.7640, 120.7640]},  # 14
            ('c2', 'c2', 'c2', 'c3'): {" ": [7.8153, 0.0000],
                                       "aa": [-0.0000, -0.0000, 0.0000, 118.9000, 120.0500, 120.0500]},  # 17
            # ('Es', 'c2', 'c2', 'c2'): {" ": [17.0526, 0.],
            #                            "aa": []}
        }}

    def __init__(self, name, atom3d: csc.Atom3D = None, atom3d_gen=None):
        self.q = False
        self.improper = True
        self.dihedral = True
        self.name = name
        self.record = None
        self.log = None
        self.atom3d_gen = atom3d_gen
        self.atom3d = self.atom3d_gen() if atom3d is None else atom3d
        self.origin_atom3d = cp.deepcopy(self.atom3d)
        self.iter = 0
        self.path = fps.get_project_path() + "data/" + self.name + "/"
        if not os.path.exists(self.path):
            os.makedirs(self.path)
        self.para_style = dict()
        self.coefficient = dict()
        for para in self.para_list:
            self.para_style[para] = dict()
            self.coefficient[para] = dict()
        for atom_type1 in self.atom3d.Atom_type_dict:
            for atom_type2 in self.atom3d.Atom_type_dict:
                type_tuple, _ = csc.sort_types((atom_type1, atom_type2))
                if type_tuple not in self.coefficient["non_bond"]:
                    self.para_style["non_bond"][type_tuple] = "table"
                    self.coefficient["non_bond"][type_tuple] = list()
        for para in self.para_list[1:]:
            for type_tuple in self.atom3d.__dict__[para + "_type_dict"]:
                self.para_style[para][type_tuple] = self.default_style[para]
                # if para == "Bond" and (type_tuple[0][:2] == "TO" or type_tuple[1][:2] == "TO"):
                # if para == "Bond" and (type_tuple[0][0].isupper() or type_tuple[1][0].isupper()):
                #     self.para_style[para][type_tuple] = "fene/expand"
                # if para == "Dihedral" and (
                #         type_tuple[:-1] == ("c2", "c2", "c2") or type_tuple[1:] == ("c2", "c2", "c2")):
                #     self.para_style[para][type_tuple] = "harmonic_d"
                if type_tuple in self.class2_type[para]:
                    self.para_style[para][type_tuple] = "class2"
                self.coefficient[para][type_tuple] = list()

    def cal_gr_serial(self, trj_path, step, aa=None, exclude=None, non_bond=True):
        if aa:
            np.random.seed(step)
            st = np.random.uniform(0, 5)
            # print("sleep: %f" % st)
            time.sleep(st)
            aa.renew_coordinate(trj_path, step=step)
            self.atom3d.renew_coordinate(aa)
        else:
            self.atom3d.renew_coordinate_file(trj_path + ".%d" % step)
        result_dict = dict()
        for para in self.para_list:
            result_dict[para] = dict()
            for type_tuple in self.coefficient[para]:
                if para == "non_bond":
                    if non_bond:
                        result_dict[para][type_tuple] = self.atom3d.cal_inter_para(type_tuple, self.cal_points[para],
                                                                               exclude=exclude)
                    else:
                        result_dict[para][type_tuple] = np.zeros(len(self.cal_points[para])+200)
                else:
                    result_dict[para][type_tuple] = self.atom3d.cal_intra_para(para, type_tuple, self.cal_points[para])
        lp = self.atom3d.lattice_parameter
        result_dict["volume"] = lp[0] * lp[1] * lp[2]
        # print("finish", step)
        return result_dict

    def cal_gr(self, trj_path, out_path, begin_step, frames, each_step=1000, num_parallel=40, cg=False, exclude=None, non_bond=True):
        if cg:
            aa = None
            self.log.info("Cal gr")
        else:
            aa = csc.AA3D()
            aa.get_mass(trj_path[:-9] + "data")
        if not os.path.exists(out_path):
            os.makedirs(out_path)
        result_dict = dict()
        result_dict["volume"] = 0.
        for para_name in self.para_list:
            result_dict[para_name] = dict()
            for type_tuple in self.coefficient[para_name]:
                if para_name == "non_bond" and non_bond:
                    result_dict[para_name][type_tuple] = np.zeros(len(self.cal_points[para_name]) + 200)
                else:
                    result_dict[para_name][type_tuple] = np.zeros(len(self.cal_points[para_name]))

        if num_parallel > 0:
            pools = Pool(num_parallel)
            parallel = list()
            for i in range(frames):
                parallel.append(
                    pools.apply_async(self.cal_gr_serial, args=(trj_path, begin_step + i * each_step, aa, exclude, non_bond)))
                # parallel.append(self.cal_target_serial(aa, aa_trj_path, begin_step + i * each_step))

            pools.close()
            pools.join()
            for result in parallel:
                new_result_dict = result.get()
                # new_result_dict = result
                for para_name in self.para_list:
                    if para_name in self.para_list:
                        for type_tuple in result_dict[para_name]:
                            result_dict[para_name][type_tuple] += new_result_dict[para_name][type_tuple]
                result_dict["volume"] += new_result_dict["volume"]
        else:
            for i in range(1):
                new_result_dict = self.cal_gr_serial(trj_path, begin_step + i * each_step, aa, exclude, non_bond)
                for para_name in self.para_list:
                    if para_name in self.para_list:
                        for type_tuple in result_dict[para_name]:
                            result_dict[para_name][type_tuple] += new_result_dict[para_name][type_tuple]
                result_dict["volume"] += new_result_dict["volume"]

        for para_name in result_dict:
            if para_name in self.para_list:
                for type_tuple in result_dict[para_name]:
                    x = self.cal_points[para_name]
                    y = result_dict[para_name][type_tuple] / frames
                    step = x[1] - x[0]
                    if para_name == "non_bond" or para_name == "Bond":
                        bandwidth = 0.1
                    else:
                        bandwidth = 3.
                    g_points = np.arange(-200 * step, 200 * step, step)
                    s = 2 * np.square(bandwidth)
                    gaussian = 1 / np.sqrt(s * np.pi) * np.exp(-np.square(g_points) / s)
                    if para_name == "non_bond":
                        if not non_bond: continue
                        gr = np.convolve(y, gaussian, mode="same")
                        gr = gr[:len(self.cal_points[para_name])] / 100.
                        gr /= (len(self.atom3d.Atom_type_dict[type_tuple[1]]) / result_dict["volume"] * frames)
                    else:
                        if para_name == "Angle" or para_name == "Dihedral":
                            yf = np.append(y[200:0:-1], y)
                            yf = np.append(yf, y[-1:-201:-1])
                            gr = np.convolve(yf, gaussian, mode="same")
                            gr = gr[200:-200]
                        else:
                            gr = np.convolve(y, gaussian, mode="same")
                        gr = gr / np.sum(gr) / (x[1] - x[0])

                    data = pd.DataFrame(OrderedDict({'x': x, 'y': y[:len(self.cal_points[para_name])], 'gr': gr}))
                    result_dict[para_name][type_tuple] = data
                    if cg:
                        data.to_csv("%s%s_%s.csv" % (out_path, para_name, '_'.join(type_tuple)), index=False)
                    else:
                        data.to_csv("%s%s_%s_target.csv" % (out_path, para_name, '_'.join(type_tuple)), index=False)
        with open(out_path + "log.txt", "w") as file:
            log = list()
            log.append(self.name + "\n")
            log.append(fps.get_date() + "\n")
            log.append(trj_path + "\n")
            log.append("Frames: %d\n" % frames)
            log.append("Average volume: %.3f\n" % (result_dict["volume"] / frames))
            file.writelines(log)
        return result_dict

    def cal_gr_serial_new(self, atom3ds, steps, exclude=None):
        atom3d = atom3ds(steps[0], fake=False)
        result_dict = dict()
        for para in self.para_list:
            result_dict[para] = dict()
            for type_tuple in self.coefficient[para]:
                # print(self.coefficient)
                if para == "non_bond":
                    # print(type_tuple)
                    result_dict[para][type_tuple] = atom3d.cal_inter_para(type_tuple, self.cal_points[para],
                                                                          exclude=exclude)
                else:
                    result_dict[para][type_tuple] = atom3d.cal_intra_para(para, type_tuple, self.cal_points[para])
        lp = atom3d.lattice_parameter
        result_dict["volume"] = lp[0] * lp[1] * lp[2]

        for i in steps[1:]:
            atom3d = atom3ds(i, fake=False)
            for para in self.para_list:
                for type_tuple in self.coefficient[para]:
                    # print(self.coefficient)
                    if para == "non_bond":
                        # print(type_tuple)
                        result_dict[para][type_tuple] += atom3d.cal_inter_para(type_tuple, self.cal_points[para],
                                                                               exclude=exclude)
                    else:
                        result_dict[para][type_tuple] += atom3d.cal_intra_para(para, type_tuple, self.cal_points[para])
            lp = atom3d.lattice_parameter
            result_dict["volume"] += lp[0] * lp[1] * lp[2]
        # print("finish", step)
        return result_dict

    def cal_gr_new(self, atom3d_gen, out_path, steps, num_parallel=40, exclude=None, name="cg"):
        os.makedirs(out_path, exist_ok=True)
        result_dict = dict()
        result_dict["volume"] = 0.
        for para_name in self.para_list:
            result_dict[para_name] = dict()
            for type_tuple in self.coefficient[para_name]:
                if para_name == "non_bond":
                    result_dict[para_name][type_tuple] = np.zeros(len(self.cal_points[para_name]) + 200)
                else:
                    result_dict[para_name][type_tuple] = np.zeros(len(self.cal_points[para_name]))

        pools = Pool(num_parallel)
        parallel = list()
        steps = list(steps)
        s_len = len(steps)
        block = s_len // num_parallel
        for i in range(num_parallel):
            # print("call")
            # parallel.append(pools.apply_async(self.cal_gr_serial_new, args=(cp.deepcopy(atom3d_gen(i, fake=False)), exclude)))
            if i < num_parallel - 1:
                parallel.append(pools.apply_async(self.cal_gr_serial_new,
                                                  args=(atom3d_gen, steps[i * block:(i + 1) * block], exclude)))
            else:
                parallel.append(
                    pools.apply_async(self.cal_gr_serial_new, args=(atom3d_gen, steps[i * block:], exclude)))
        pools.close()
        pools.join()
        for result in parallel:
            new_result_dict = result.get()
            # new_result_dict = result
            for para_name in self.para_list:
                if para_name in self.para_list:
                    for type_tuple in result_dict[para_name]:
                        result_dict[para_name][type_tuple] += new_result_dict[para_name][type_tuple]
            result_dict["volume"] += new_result_dict["volume"]

        for para_name in result_dict:
            if para_name in self.para_list:
                for type_tuple in result_dict[para_name]:
                    x = self.cal_points[para_name]
                    y = result_dict[para_name][type_tuple] / len(list(steps))
                    step = x[1] - x[0]
                    if para_name == "non_bond" or para_name == "Bond":
                        bandwidth = 0.1
                    else:
                        bandwidth = 3.
                    g_points = np.arange(-200 * step, 200 * step, step)
                    s = 2 * np.square(bandwidth)
                    gaussian = 1 / np.sqrt(s * np.pi) * np.exp(-np.square(g_points) / s)
                    if para_name == "non_bond":
                        gr = np.convolve(y, gaussian, mode="same")
                        gr = gr[:len(self.cal_points[para_name])] / 100.
                        gr /= (len(self.atom3d.Atom_type_dict[type_tuple[1]]) / result_dict["volume"] * len(
                            list(steps)))
                    else:
                        if para_name == "Angle" or para_name == "Dihedral":
                            yf = np.append(y[200:0:-1], y)
                            yf = np.append(yf, y[-1:-201:-1])
                            gr = np.convolve(yf, gaussian, mode="same")
                            gr = gr[200:-200]
                        else:
                            gr = np.convolve(y, gaussian, mode="same")
                        gr = gr / np.sum(gr) / (x[1] - x[0])

                    data = pd.DataFrame(OrderedDict({'x': x, 'y': y[:len(self.cal_points[para_name])], 'gr': gr}))
                    result_dict[para_name][type_tuple] = data
                    data.to_csv("%s%s_%s_%s.csv" % (out_path, para_name, '_'.join(type_tuple), name), index=False)
        with open(out_path + "log.txt", "w") as file:
            log = list()
            log.append(self.name + "\n")
            log.append(fps.get_date() + "\n")
            file.writelines(log)
        return result_dict

    def get_target_dict(self, path):
        self.log.info("Get target gr")
        target_dict = dict()
        for para in self.para_list:
            target_dict[para] = dict()
            for type_tuple in self.coefficient[para]:
                target_file_name = "%s%s_%s_target.csv" % (path, para, "_".join(type_tuple))
                target_dict[para][type_tuple] = pd.read_csv(target_file_name)
        return target_dict

    def init_non_bond(self, target_dict, init="same"):
        points = self.cal_points["non_bond"]
        for type_tuple in target_dict["non_bond"]:
            if type_tuple[0] in self.aa_type and type_tuple[1] in self.aa_type:
                self.coefficient["non_bond"][type_tuple] = f.cal_aa_force(points, type_tuple,
                                                                          cut_off=self.non_bond_cutoff, q=True,
                                                                          dielectric=1., dsf=True)
            # elif type_tuple[0] in self.aa_type or type_tuple[1] in self.aa_type:
            #     if type_tuple[0][0].islower():
            #         cg_type = ("U", type_tuple[1])
            #     else:
            #         cg_type = (type_tuple[0], "U")
            #     cg_path = "/home/centos/Projects/CGMD/data/target/1blk_50_600K_cg_ex_hard/aver/"
            #     cg_data = pd.read_csv(cg_path + "non_bond_%s_target.csv" % "_".join(cg_type))
            #     cg_m = f.find_maximum(cg_data["gr"])[0]
            #     yr = np.array(target_dict["non_bond"][type_tuple].gr)
            #     cgu_m = f.find_maximum(yr)[0]
            #     yr = yr * yr[cgu_m] / cg_data["gr"][cg_m]
            #     l = np.argwhere(yr > 1E-7)[0][0]
            #     logy = np.log(yr)
            #     e_new = -8.314 * 300 / 4186 * logy
            #     fit_e = self.fit_func["non_bond"](points, e_new, l) * f.smooth_charmm(points,
            #                                                                           self.non_bond_cutoff - 4.,
            #                                                                           self.non_bond_cutoff)
            #     self.coefficient["non_bond"][type_tuple] = fit_e
            else:
                yr = np.array(target_dict["non_bond"][type_tuple].gr)
                l = np.argwhere(yr > 1E-7)[0][0]
                logy = np.log(yr)
                e_new = -8.314 * 300 / 4186 * logy
                fit_e = self.fit_func["non_bond"](points, e_new, l) * f.smooth_charmm(points,
                                                                                      self.non_bond_cutoff - 4.,
                                                                                      self.non_bond_cutoff)
                self.coefficient["non_bond"][type_tuple] = fit_e * 0.3
        # else:
        #     raise ValueError("init should be same or gr")
        pass

    def init_parameter(self, target_dict, init="same"):
        self.log.info("Init parameters")
        for para in self.para_list:
            for type_tuple in self.coefficient[para]:
                if para != "non_bond":
                    style = self.para_style[para][type_tuple]
                    if style == "class2":
                        continue
                    # print(type_tuple)
                    x = np.array(target_dict[para][type_tuple].x)
                    y = np.array(target_dict[para][type_tuple].gr)
                    # yr = y / y.sum() / (x[1] - x[0])
                    # yr += 1E-5
                    # if para == "Bond":
                    l_id = np.argwhere(y > 1E-5)[0][0]
                    r_id = np.argwhere(y > 1E-5)[-1][0]
                    xr = x[l_id:r_id]
                    yr = y[l_id:r_id]
                    yr += 1E-300
                    # else:
                    #     yr = y + 1E-300
                    #     xr = x
                    # xr = x
                    # print(type_tuple)
                    if para == "Bond":
                        E = -8.314 * 298 / 4186 * np.log(yr / xr / xr)
                    else:
                        E = -8.314 * 298 / 4186 * np.log(yr)
                        # print(E)
                    # E = -8.314 * 298 / 4186 * np.log(yr)
                    E -= E.min()
                    param = self.fit_func[para](xr, E, style)
                    # if para == "Dihedral":
                    #     plt.plot(xr, E)
                    #     sf = self.style_func[self.para_style[para][type_tuple]]
                    #     #     if self.para_style[para][type_tuple] == "harmonic":
                    #     #         plt.plot(xr, sf(xr, param[0]/180/180*np.pi*np.pi, param[1]))
                    #     #     else:
                    #     plt.plot(xr, sf(xr, *param))
                    #     plt.title("_".join(type_tuple))
                    #     plt.show()
                    #     print(type_tuple)
                    #     print(param)
                    # if para == "Angle":
                    #     plt.plot(xr, sf(xr / 180 * np.pi, param[0], param[1] / 180 * np.pi), label="fit")
                    # if para == "Dihedral":
                    #     print(param)
                    #     plt.plot(xr, E, label="origin")
                    #     plt.plot(xr, sf(xr, *param), label="fit")
                    #     plt.title("_".join(type_tuple))
                    #     plt.legend()
                    #     plt.show()
                    #     plt.close()
                    # import time
                    # time.sleep(0.5)
                    self.coefficient[para][type_tuple] = param
        # c2 angle
        # self.coefficient["Angle"][("Es", "c2", "c2")] = cp.copy(self.coefficient["Angle"][("c2", "c2", "c2")])
        if ("Ph", "U", "Ph") in self.coefficient["Angle"].keys():
            self.coefficient["Angle"][("Ph", "U", "Ph")] = [63.52, 117.4]
        self.init_non_bond(target_dict, init)
        self.record.renew_state(-1, parameters=self.coefficient)

    def create_in_file(self, out_path, in_name):
        def format_param(pms):
            pm_str = ""
            for pm in pms:
                pm_str += "{:.5f} ".format(pm)
            return pm_str

        def potential_param():
            bond_string = "# potentials for %s, iter %d\n" % (self.name, self.iter)
            bond_string += "bond_style      hybrid harmonic fene/expand class2\n"
            bond_type_list = list(self.atom3d.Bond_type_dict.keys())
            bond_type_list.sort()
            # print(self.para_style)
            for joint_id, joint_type in enumerate(bond_type_list):
                style = self.para_style["Bond"][joint_type]
                if style == "class2":
                    p = cp.copy(self.class2_type["Bond"][joint_type])
                else:
                    p = cp.copy(self.coefficient["Bond"][joint_type])
                if style == "fene/expand":
                    p = [p[0], p[2], 0., 0., p[1]]
                bond_string += "bond_coeff     {:<2} {:<12}  ".format(joint_id + 1, style)
                bond_string += format_param(p)
                #
                # if style == "harmonic":
                #     p = self.coefficient["Bond"][joint_type]
                #     string += "bond_coeff     {:<2} {:<12}  {:.5f} {:.5f} ".format(joint_id + 1, style, *p)
                # elif style == "fene/expand":
                #     p = self.coefficient["Bond"][joint_type]
                #     string += "bond_coeff     {:<2} {:<12}  {:.5f} {:.5f} 0.0 0.0 {:.5f}".format(joint_id + 1, style, p[0],
                #                                                                              p[2], p[1])
                # elif style == "class2":
                #     string += "bond_coeff     {:<2} {:<12}  {:.5f} {:.5f} {:.5f} {:.5f} ".format(joint_id + 1, style, *p)
                # else:
                #     raise KeyError("No such bond style!")

                bond_string += " # %s \n" % "_".join(joint_type)

            bond_string += "\nangle_style     hybrid harmonic class2\n"
            # string += "\nangle_style     hybrid table linear 3601 harmonic\n"
            # string += "\nangle_style     table linear\n"
            angle_type_list = list(self.atom3d.Angle_type_dict.keys())
            angle_type_list.sort()
            for joint_id, joint_type in enumerate(angle_type_list):
                style = self.para_style["Angle"][joint_type]
                if style == "class2":
                    p = cp.copy(self.class2_type["Angle"][joint_type])
                    for k, v in p.items():
                        bond_string += "angle_coeff     {:<2} {:<12} {:<4} ".format(joint_id + 1, style, k)
                        bond_string += format_param(v) + "\n"
                else:
                    p = cp.copy(self.coefficient["Angle"][joint_type])
                    bond_string += "angle_coeff     {:<2} {:<12}  ".format(joint_id + 1, style)
                    bond_string += format_param(p)
                # if style == "harmonic":
                # string += "angle_coeff    {:>2} {:.5f} {:.5f} ".format(joint_id + 1,
                #                                                        *self.coefficient["Angle"][joint_type])
                bond_string += " # %s \n" % "_".join(joint_type)

            if self.dihedral:
                bond_string += "\ndihedral_style     hybrid multi/harmonic harmonic class2\n"
                dihedral_type_list = list(self.atom3d.Dihedral_type_dict.keys())
                dihedral_type_list.sort()
                for joint_id, joint_type in enumerate(dihedral_type_list):
                    # if self.para_style["Dihedral"][joint_type] == "harmonic_d":
                    #     string += "dihedral_coeff    {:>2} harmonic {:.5f} {:.0f} 1" \
                    #         .format(joint_id + 1, *self.coefficient["Dihedral"][joint_type])
                    # else:
                    #     string += "dihedral_coeff    {:>2} multi/harmonic {:.5f} {:.5f} {:.5f} {:.5f} {:.5f}" \
                    #         .format(joint_id + 1, *self.coefficient["Dihedral"][joint_type])
                    style = self.para_style["Dihedral"][joint_type]
                    if style == "class2":
                        p = cp.copy(self.class2_type["Dihedral"][joint_type])
                        for k, v in p.items():
                            bond_string += "dihedral_coeff     {:<2} {:<12} {:<4} ".format(joint_id + 1, style, k)
                            bond_string += format_param(v) + "\n"
                    else:
                        p = cp.copy(self.coefficient["Dihedral"][joint_type])
                        bond_string += "dihedral_coeff     {:<2} {:<12}  ".format(joint_id + 1, style)
                        bond_string += format_param(p)
                    bond_string += " # %s \n" % "_".join(joint_type)
            # string += "\nimproper_style     class2\n"
            #             # improper_type_list = list(self.atom3d.Improper_type_dict.keys())
            #             # improper_type_list.sort()
            #             # for joint_id, joint_type in enumerate(improper_type_list):
            #             #     string += "improper_coeff    {:>2}  {:.5f} {:.5f}" \
            #             #         .format(joint_id + 1, *self.improper_coefficient[joint_type])
            #             #     string += " # %s \n" % "_".join(joint_type)
            #             #         string += '''
            #             # improper_style     class2
            #             # improper_coeff     1  17.05260 0.00000 # Es_c2_c2_c2
            #             # improper_coeff     * aa 0.0000     5.9863     0.0000   116.0640   116.0640   118.9000
            #             # improper_coeff     2  7.81530 0.00000 # c2_c2_c2_c1
            #             # improper_coeff     * aa -0.0000    -0.0000     0.0000   118.9000   120.0500   120.0500
            #             # improper_coeff     3  17.05260 0.00000 # c2_c2_c2_n1
            #             # improper_coeff     * aa 0.0000     0.0000     0.0000   118.9000   120.7640   120.7640
            #             # improper_coeff     4  59.37400 0.00000 # n1_c3_n1_o1
            #             # improper_coeff     * aa -0.0001    -0.0001    -0.0001   120.5292   125.5320   125.5320
            #             # '''
            #             # new type
            if self.improper:
                bond_string += "\nimproper_style     class2\n"
                improper_type_list = list(self.atom3d.Improper_type_dict.keys())
                improper_type_list.sort()
                for joint_id, joint_type in enumerate(improper_type_list):
                    p = cp.copy(self.class2_type["Improper"][joint_type])
                    for k, v in p.items():
                        bond_string += "improper_coeff     {:<2} {:<4} ".format(joint_id + 1, k)
                        bond_string += format_param(v) + "\n"
                    bond_string += " # %s \n" % "_".join(joint_type)
            #         string += '''
            # improper_style     class2
            # improper_coeff     1  59.37400 0.00000 # n1_c1_n1_o1
            # improper_coeff     1  aa -0.0001    -0.0001    -0.0001   120.5292   125.5320   125.5320
            # '''

            # string += "\npair_style      table linear %d\n" % len(self.cal_points["non_bond"] + 1)
            if in_name != "dpd.in":
                if self.q:
                    # string += "\npair_style      table linear %d pppm\n" % int(self.param/non_bond_cutoff * 100)
                    bond_string += "\npair_style      hybrid/overlay table linear %d coul/dsf 0.01 12.\n" % self.cut_idx
                else:
                    if self.improper:
                        bond_string += "\npair_style      hybrid table linear %d table linear %d\n" % (
                            self.cut_idx, self.cut_idx)
                    else:
                        bond_string += "\npair_style      table linear %d\n" % self.cut_idx
                atom_type_list = list(self.atom3d.Atom_type_dict.keys())
                atom_type_list.sort()
                for type1, type2 in self.coefficient["non_bond"].keys():
                    if type1 not in atom_type_list or type2 not in atom_type_list:
                        continue
                    type_list = [atom_type_list.index(type1) + 1, atom_type_list.index(type2) + 1]
                    if self.q:
                        bond_string += "pair_coeff     {:>2} {:>2} table  {:<35} {}\n" \
                            .format(type_list[0], type_list[1], "param/non_bond_%s_%s.param" % (type1, type2),
                                    "%s_%s" % (type1, type2))
                    else:
                        if self.improper:
                            flag = 2 if type1[0].islower() and type2[0].islower() else 1
                            bond_string += "pair_coeff     {:>2} {:>2} table {} {:<35} {}\n" \
                                .format(type_list[0], type_list[1], flag, "param/non_bond_%s_%s.param" % (type1, type2),
                                        "%s_%s" % (type1, type2))
                        else:
                            bond_string += "pair_coeff     {:>2} {:>2} {:<35} {}\n" \
                                .format(type_list[0], type_list[1], "param/non_bond_%s_%s.param" % (type1, type2),
                                        "%s_%s" % (type1, type2))
                if self.q:
                    bond_string += "pair_coeff      * * coul/dsf\n"
                    with open(self.path + "record/dielectric_iter_%02d.txt" % self.iter, "r") as file:
                        di = float(file.read().split()[0])
                    bond_string += "dielectric      %.5f\n\n" % di
            if self.improper:
                bond_string += '''special_bonds   lj/coul 0 0.9999999999 1
pair_modify		pair table 1 special lj/coul 0 0.9999999999 1
pair_modify     pair table 2 special lj/coul 0 1e-10 1
            '''
            else:
                bond_string += "special_bonds   lj/coul 0 1 1\n"
            with open(out_path + "potential.param", "w") as file:
                file.write(bond_string)

        string = "# Initialization\n"
        string += "variable        n string %s\n" % "1blk_40"
        string += "log             $n.log\n"
        string += "units           real\n"
        if self.q:
            string += "atom_style      full\n"
        else:
            string += "atom_style      molecular\n"
        string += "boundary        p p p \n\n" \
                  "##############################################\n"
        string += "# Forcefield parameters\n" \
                  "read_data       $n.data\n" \
                  "include         potential.param\n\n"
        potential_param()
        with open(fps.get_project_path() + "template/%s" % in_name) as template_file:
            template_string = template_file.read()
        string += template_string
        with open(out_path + in_name, "w") as in_file:
            in_file.write(string)
        pass

    def create_param_file(self, para_name, type_tuple, out_path):
        def out_put(_x, _y):
            dy = f.cal_derivative(_x, _y)
            string = "_".join(type_tuple) + "\n"
            string += "N %d\n\n" % len(_x)
            for i in range(len(_x)):
                string += "%d  %.6f  %.12f  %.10f\n" % (i, _x[i], _y[i], -dy[i])
            with open(out_path + "param/" + para_name + "_" + "_".join(type_tuple) + ".param", "w") as param_file:
                param_file.write(string)

        # print(type_tuple)
        os.makedirs(out_path + "param/", exist_ok=True)
        if para_name == "non_bond":
            x = cp.copy(self.cal_points[para_name])
            y = cp.copy(self.coefficient[para_name][type_tuple])
            x = np.append(1e-5, x)
            y = np.append(2 * y[0] + y[3] - 2 * y[2], y)
            out_put(x[:self.cut_idx + 1], y[:self.cut_idx + 1])
        pass

    def prepare_simulation_file(self, file_path, in_name):
        if self.atom3d_gen is not None:
            self.origin_atom3d = self.atom3d_gen()
        fps.LmpParser.create_data_file(self.origin_atom3d, file_path + "1blk_40.data", q=self.q, improper=self.improper, dihedral=self.dihedral)
        os.makedirs(file_path + "trj/", exist_ok=True)
        self.create_in_file(file_path, in_name)
        for para in self.para_list:
            for type_tuple in self.coefficient[para]:
                self.create_param_file(para, type_tuple, file_path)

    def invoke_lmp_process(self, file_path, in_name, num_parallel, screen=False):
        success = False
        time_out = 120 * 60
        while not success:
            if os.path.exists(file_path):
                shutil.rmtree(file_path)
            os.mkdir(file_path)
            self.prepare_simulation_file(file_path, in_name)
            p_lmp = Process(target=fps.invoke_lmp, args=(file_path, in_name, num_parallel, screen))
            p_lmp.start()
            start = time.time()
            while time.time() - start <= time_out:
                if p_lmp.is_alive():
                    time.sleep(1)
                else:
                    success = True
                    if not os.path.exists(file_path + "finish.restart"):
                        success = False
                        self.log.error("error, recalculate trajectory")
                        p_lmp.terminate()
                        p_lmp.join()
                        time.sleep(10)
                    break
            else:
                print("time out, recalculate trajectory")
                p_lmp.terminate()
                p_lmp.join()
                time.sleep(10)

    def pre_adjust(self, num_parallel, with_pre=True):
        self.log.info("Pre Adjust")

        for type_tuple in self.coefficient["non_bond"].keys():
            # first order correction
            # if type_tuple[0] not in self.aa_type and type_tuple[1] not in self.aa_type:
            self.coefficient["non_bond"][type_tuple] -= self.coefficient["non_bond"][type_tuple][self.cut_idx - 1]
            # second order correction
            # self.coefficient["non_bond"][type_tuple] *= f.smooth_step(self.cal_points["non_bond"], 10., 12.)
            # k = (self.coefficient["non_bond"][type_tuple][cut_idx - 1] - self.coefficient["non_bond"][type_tuple][
            #     cut_idx]) / (self.cal_points["non_bond"][1] - self.cal_points["non_bond"][0])
            # self.coefficient["non_bond"][type_tuple] -= k * self.non_bond_cutoff - k * self.cal_points["non_bond"]

            self.coefficient["non_bond"][type_tuple] = np.append(
                self.coefficient["non_bond"][type_tuple][:self.cut_idx],
                np.zeros(len(self.cal_points["non_bond"]) - self.cut_idx))
        if with_pre:
            file_path = self.path + "iter_%02d/pre/" % self.iter
            converge = False
            for i in range(1):
                # if self.iter > 0 or i > 0:
                self.invoke_lmp_process(file_path, "pre.in", num_parallel, screen=False)
                with open(file_path + "1blk_40.log") as log_f:
                    log_lines = log_f.readlines()
                begin_id = 0
                recreate_flag = 0
                for log_id, log_line in enumerate(log_lines):
                    sp = log_line.split()
                    if sp:
                        if sp[0] == "Step":
                            begin_id = log_id + 1
                        elif sp[0] == "WARNING:":
                            recreate_flag += 1
                pressure = 0.
                for log_line in log_lines[begin_id: begin_id + 200]:
                    pressure += float(log_line.split()[2])
                pressure /= 200.
                self.log.info("pressure: %.2f" % pressure)
                if abs(pressure) > 50:
                    for type_tuple in self.coefficient["non_bond"].keys():
                        if type_tuple[0] not in self.aa_type and type_tuple[1] not in self.aa_type:
                            self.coefficient["non_bond"][type_tuple] -= (-self.cal_points[
                                "non_bond"] / self.non_bond_cutoff + 1) * pressure / 3E4
                else:
                    converge = True
                    break
            # if not converge:
            #     raise RuntimeError("Can not adjust pressure!")

    def cg_simulation(self, num_parallel, ensemble="nvt"):
        self.log.info("CG simulation")
        file_path = self.path + "iter_%02d/cal/" % self.iter
        if not os.path.exists(file_path):
            os.makedirs(file_path)
        self.invoke_lmp_process(file_path, "%s.in" % ensemble, num_parallel=num_parallel, screen=False)

    def cg_gr(self, num_parallel=40, lmp_parallel=12, with_pre=True, ensemble="nvt", non_bond=True):
        def prepare_file_path(fp):
            if not os.path.exists(fp):
                os.makedirs(fp)

        file_path = self.path + "iter_%02d/" % self.iter
        result_dict = None
        if self.record.restart:
            self.coefficient = self.record.get_param()
        while self.record.current_state[1:] != [True] * 3:
            if not self.record.current_state[1]:
                prepare_file_path(file_path)
                self.pre_adjust(lmp_parallel, with_pre)
                self.record.renew_state(1, self.coefficient)
            elif not self.record.current_state[2]:
                prepare_file_path(file_path)
                self.cg_simulation(lmp_parallel, ensemble)
                self.record.renew_state(2)
            elif not self.record.current_state[3]:
                self.atom3d = fps.LmpParser.load_data_file(file_path + "cal/1blk_40.data")
                # result_dict = self.cal_gr(file_path + "cal/trj/300K.lammpstrj", file_path + "csv/", 0, 1000, each_step=500,
                result_dict = self.cal_gr(file_path + "cal/trj/300K.lammpstrj", file_path + "csv/", 0, 1000,
                                          each_step=200,
                                          num_parallel=num_parallel, cg=True, non_bond=non_bond)
                self.record.renew_state(3)
        if result_dict is None:
            result_dict = dict()
            for para in self.coefficient:
                result_dict[para] = dict()
                for type_tuple in self.coefficient[para]:
                    data_file = "%scsv/%s_%s.csv" % (file_path, para, '_'.join(type_tuple))
                    result_dict[para][type_tuple] = pd.read_csv(data_file)
        return result_dict

    def fit_all_param(self, target, fit, learning_rate, num_parallel=48, non_bond=True):
        self.log.info("Fit parameters")
        error_dict = dict()
        # param_20 = self.record.get_param(self.iter - 20)
        type_tuple_list = list()
        error_list = list()
        for para_name in self.para_list:
            error_dict[para_name] = dict()
            for type_tuple in self.coefficient[para_name]:
                # _error = float(
                #     f.error_l2(np.array(fit[para_name][type_tuple].gr), np.array(target[para_name][type_tuple].gr)))
                _error = float(
                    f.error_l1(np.array(fit[para_name][type_tuple].gr), np.array(target[para_name][type_tuple].gr)))
                error_dict[para_name][type_tuple] = _error
                if para_name == "non_bond":
                    if not (type_tuple[0] in self.aa_type and type_tuple[1] in self.aa_type):
                        error_list.append(_error)
                        type_tuple_list.append(type_tuple)
                else:
                    self.log.info("%s, %s: %.3f" % (para_name, "_".join(type_tuple), _error))
        error_list = np.array(error_list)
        mean_error = float(np.sqrt(np.square(error_list).mean()))
        self.log.info("non_bond msq error: %.3f" % mean_error)
        min_idx = int(np.argmin(error_list))
        max_idx = int(np.argmax(error_list))
        self.log.info("non_bond min error: %.3f, %s" % (error_list[min_idx], "_".join(type_tuple_list[min_idx])))
        self.log.info("non_bond max error: %.3f, %s" % (error_list[max_idx], "_".join(type_tuple_list[max_idx])))
        # print(type(error_dict["non_bond"][("c2", "c2")]))
        new_di = 1.
        # h_dif = np.array(fit["non_bond"][("o1", "o1")].gr).max() / np.array(target["non_bond"][("o1", "o1")].gr).max()
        # # h_dif = np.argmax(np.array(target["non_bond"][("h1", "o1")].gr)) - np.argmax(
        # #     np.array(fit["non_bond"][("h1", "o1")].gr))
        # with open(self.path + "record/dielectric_iter_%02d.txt" % self.iter, "r") as file:
        #     di = float(file.read().split()[0])
        # new_di = di * (1 + (h_dif - 1) * 0.05)
        # # if self.iter == 1:
        # #     new_di = 1.
        # with open(self.path + "record/dielectric_iter_%02d.txt" % (self.iter + 1), "w") as file:
        #     file.write(str(new_di))
        if not os.path.exists(self.path + "iter_%02d/png/" % self.iter):
            os.makedirs(self.path + "iter_%02d/png/" % self.iter)
        pool = Pool(num_parallel)
        parallel = list()
        for para_name in self.para_list:
            if (not non_bond) and para_name == "non_bond":
                continue
            for type_tuple in self.coefficient[para_name]:
                # if para_name == "Angle" and type_tuple == ("Es", "c2", "c2"):
                #     continue
                # if para_name == "Angle" and type_tuple == ("Ph", "U", "Ph"):
                #     pass
                # else:
                if type_tuple == type_tuple_list[max_idx]:
                    max_bool = False  # True
                else:
                    max_bool = False
                # if error_dict[para_name][type_tuple] > np.median(error_list):
                #     fit_bool = True
                # else:
                #     fit_bool = False
                # self.fit_parameter(para_name, type_tuple, fit, target, learning_rate, fit_bool=fit_bool,
                #                    new_di=new_di)
                parallel.append(pool.apply_async(self.fit_parameter,
                                                 args=(para_name, type_tuple, fit, target, learning_rate, max_bool,
                                                       new_di)))
        # c2 angle
        # self.coefficient["Angle"][("Es", "c2", "c2")] = cp.copy(self.coefficient["Angle"][("c2", "c2", "c2")])
        pool.close()
        pool.join()
        for result in parallel:
            para_name, type_tuple, param = result.get()
            self.coefficient[para_name][type_tuple] = param
        if ("Ph", "U", "Ph") in self.coefficient["Angle"].keys():
            self.coefficient["Angle"][("Ph", "U", "Ph")] = [63.52, 117.4]
        if ("Es", "c2", "c2") in self.coefficient["Angle"].keys():
            self.coefficient["Angle"][("Es", "c2", "c2")] = [1.737485325897969, 94.74749999999999]
        return error_dict

    def fit_parameter(self, para_name, type_tuple, result_dict, target_dict, learning_rate=0.024, max_bool=False,
                      new_di=1.):
        def fit_non_bond(_fit_data, _target_data, cupa):
            points = self.cal_points["non_bond"]
            e_present = cp.copy(cupa)
            target = _target_data
            fit = _fit_data
            # _error = f.error_l2(fit, target)
            # print(type_tuple)

            judge_t = np.where(target > 1E-60, 1, 0)
            judge_f = np.where(fit > 1E-60, 1, 0)
            judge = judge_t & judge_f
            # print(np.argwhere(judge_t > 1E-5)[0][0])
            # print(np.argwhere(judge_f > 1E-5)[0][0])
            # if type_tuple == ("c3","c3"):
            #     l = np.argwhere(target > 0.5)[0][0]
            #     lf = np.argwhere(fit > 0.5)[0][0]
            # elif type_tuple == ("TO(2)", "TO(2)"):
            #     l = np.argwhere(target > 1E-2)[0][0]
            #     lf = np.argwhere(fit > 1E-2)[0][0]
            # else:
            l = np.argwhere(target > 1E-5)[0][0]
            lf = np.argwhere(fit > 1E-5)[0][0]
            if lf - l < 400:
                if type_tuple[0] in self.aa_type and type_tuple[1] in self.aa_type:
                    ylim = None
                    aa = True
                else:
                    ylim = (-20, 20)
                    aa = False

                target += np.ones(judge.shape) - judge
                fit += np.ones(judge.shape) - judge
                dif = 8.314 * 300 / 4186 * (np.log(fit) - np.log(target)) * judge
                if type_tuple == ("TO(2)", "h1") or type_tuple == ("TO(1)", "h1"):
                    fit_m = f.find_maximum(fit)[0]
                    tar_m = f.find_maximum(target)[0]
                    if fit[fit_m] < target[tar_m]:
                        dif *= 5
                # dif *= f.linear_func(points, -2 / self.non_bond_cutoff, 2)
                dif *= f.smooth_charmm(points, self.non_bond_cutoff - 4., self.non_bond_cutoff)
                # dif *= np.sqrt(_erdif)
                # dif *= abs(min(_l_param) - min(cupa)) / 0.02

                # if type_tuple[0] != type_tuple[1]:
                #     dif *= 0.5

                # dif = dif * abs(dif) * 2 / max(abs(dif[l+50:]))
                # correction
                # dif *= (target - _l_gr) / (fit - _l_gr)
                if max_bool:
                    dif *= 1
                e_new = e_present + dif * learning_rate  # * (1 + target / target.max() * judge)

                e_new = f.fit_non_bond(points, e_new, l, aa=aa)
                # print(type_tuple, l)

                # if type_tuple == ("c3", "c3"):
                #     plt.plot(points[l-50:], e_new[l-50:])
                #     plt.show()

                # fps.plot_figs([points[l + 10::6]] * 3,
                #               [e_present[l + 10::6], e_new[l + 10::6], 10 * dif[l + 10::6] * learning_rate], 'distance',
                fps.plot_figs([points[10::6]] * 3,
                              [e_present[10::6], e_new[10::6], 10 * dif[10::6] * learning_rate], 'distance',
                              'potential',
                              'iter_%02d-%s-%s-Potential' % (self.iter, para_name, '_'.join(type_tuple)),
                              save_dictionary=self.path + "iter_%02d/png/" % self.iter,
                              label=['E_present', 'E_next', 'difX10'], show=False, ylim=ylim)
            else:
                print(type_tuple)
                print("not cal")
                e_new = e_present

            return e_new  # , _error

        def fit_bond(_fit_data, _target_data, cupa):
            points = self.cal_points[para_name]
            target = np.array(_target_data)
            fit = np.array(_fit_data)
            _error = f.error_l2(fit, target)
            # self.log.info("%s, %s: %.3f" % (para_name, "_".join(type_tuple), _error))
            if _error < 0.04:
                return cupa
            judge_min = 1E-5 if para_name == "Angle" else 1E-10
            judge_target = np.where(target > judge_min, 1, 0)
            judge_fit = np.where(fit > judge_min, 1, 0)
            judge = judge_fit & judge_target
            if para_name == "Bond" or para_name == "Angle":
                l_id = np.where(judge > 0)[0][0]
                r_id = np.where(judge > 0)[0][-1]
                points = points[l_id:r_id]
                fit = fit[l_id:r_id]
                target = target[l_id:r_id]
                judge = np.ones(points.shape)
            target += 1e-300
            fit += 1e-300
            sf = self.style_func[self.para_style[para_name][type_tuple]]
            # if para_name == "Angle" or (para_name == "Dihedral" and sf == f.harmonic):
            if para_name == "Angle":
                e_present = sf(points, cupa[0] / 180 / 180 * 3.14 * 3.14, cupa[1])
            else:
                e_present = sf(points, *cupa)
            dif = 8.314 * 300 / 4186 * (np.log(fit) - np.log(target)) * judge
            e_next = e_present + dif * 0.4
            # plt.plot(points, e_present)
            # plt.plot(points, e_next)
            # plt.title("_".join(type_tuple))
            # plt.show()
            d = points[np.argmax(target)] - points[np.argmax(fit)]
            if para_name == "Dihedral":
                d = max(target) / max(fit)
            p_new = self.fit_func[para_name](points, e_next, self.para_style[para_name][type_tuple], p0=cupa, d=d)
            # if para_name == "Angle" or (para_name == "Dihedral" and sf == f.harmonic):
            if para_name == "Angle":
                e_new = sf(points, p_new[0] / 180 / 180 * 3.14 * 3.14, p_new[1])
            else:
                e_new = sf(points, *p_new)
            fps.plot_figs([points] * 2, [e_present, e_new], 'distance', 'potential',
                          'iter_%02d-%s-%s-Potential' % (self.iter, para_name, '_'.join(type_tuple)),
                          save_dictionary=self.path + "iter_%02d/png/" % self.iter,
                          label=['fit_current', 'fit_next'], show=False)
            return p_new

        def fit_aa(di):
            if not self.q:
                E = f.cal_aa_force(point, type_tuple, cut_off=self.non_bond_cutoff, q=True, dielectric=di, dsf=True)
            else:
                E = f.cal_aa_force(point, type_tuple, cut_off=self.non_bond_cutoff, q=False)
            E = np.append(E[:self.cut_idx], np.zeros(len(point) - self.cut_idx))
            return E

        point = cp.copy(self.cal_points[para_name])
        parameters = cp.copy(self.coefficient[para_name][type_tuple])
        if para_name == "non_bond":
            target_data = np.array(target_dict[para_name][type_tuple].gr)
            fit_data = np.array(result_dict[para_name][type_tuple].gr)
            # last_param = np.array(param_20[para_name][type_tuple])
            # last_gr = np.array(pd.read_csv(self.path+"iter_%02d/csv/non_bond_%s.csv"%(self.iter-1, '_'.join(type_tuple))).gr)
            fps.plot_figs([point] * 2,
                          [target_data, fit_data], 'distance', 'g(r)',
                          'iter_%02d-%s-%s-RDF' % (self.iter, para_name, '_'.join(type_tuple)),
                          save_dictionary=self.path + "iter_%02d/png/" % self.iter, label=['target', 'fit'])
            # if para_name == "non_bond" and type_tuple[0] in self.aa_type and type_tuple[1] in self.aa_type:
            #     pass
            # else:
            if type_tuple[0] in self.aa_type and type_tuple[1] in self.aa_type:
                # self.coefficient[para_name][type_tuple] = fit_aa(new_di)
                ret_param = fit_aa(new_di)
            # elif fit_bool:
            else:
                # self.coefficient[para_name][type_tuple] = fit_non_bond(fit_data, target_data, parameters)  # , last_gr
                ret_param = fit_non_bond(fit_data, target_data, parameters)  # , last_gr

        else:
            target_data = np.array(target_dict[para_name][type_tuple].gr)
            fit_data = np.array(result_dict[para_name][type_tuple].gr)
            fps.plot_figs([point] * 2, [target_data, fit_data], 'distance', 'g(r)',
                          'iter_%02d-%s-%s-RDF' % (self.iter, para_name, '_'.join(type_tuple)),
                          save_dictionary=self.path + "iter_%02d/png/" % self.iter, label=['target', 'fit'])
            # self.coefficient[para_name][type_tuple] = fit_bond(fit_data, target_data, parameters)
            if self.para_style[para_name][type_tuple] == "class2":
                ret_param = []
            else:
                ret_param = fit_bond(fit_data, target_data, parameters)
        return para_name, type_tuple, ret_param
        # return error_sum

    def iterative_fit_params(self, max_iter, target_path, learning_rate=0.04, num_parallel=20, lmp_parallel=18,
                             restart=False, ensemble="nvt", q=False, improper=True, with_pre=True, non_bond=True):
        self.q = q
        self.improper = improper
        self.record = Record(self.path + "record/", restart=restart)
        self.log = fps.create_log(self.path + "record/", restart=restart)
        if not restart:
            self.log.info("Project: %s" % self.name)
        else:
            self.log.info("Project restart")
        self.log.info("Learning rate: %.3f" % learning_rate)
        target = self.get_target_dict(target_path)
        self.iter = self.record.current_state[0]
        if not restart:
            self.init_parameter(target)
        while self.iter < max_iter:
            self.log.info("Iter %d: " % self.iter)
            fit = self.cg_gr(num_parallel, lmp_parallel, ensemble=ensemble, with_pre=with_pre, non_bond=non_bond)
            error = self.fit_all_param(target, fit, learning_rate, num_parallel=num_parallel, non_bond=non_bond)

            self.iter += 1
            self.record.renew_state(0, self.coefficient, error)


if __name__ == '__main__':
    fp = "/home/centos/model/aa/1blk_40/0/"
    atom3d = csc.create_atom3d(fp, 40)

    ibi = IBI("test", atom3d)
    ibi.iterative_fit_params(100, "../1blk_40_target_exclude_hard_convolve/aver/", restart=False, num_parallel=48,
                             lmp_parallel=24)
