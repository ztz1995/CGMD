from pyMD import cartesian_operation as cso
from pyMD import file_parser as fps
import numpy as np
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt


class U:
    def __init__(self, c3_id):
        self.o1 = 0
        self.c3 = c3_id
        self.h1 = list()
        self.id_co = dict()

    @property
    def c3_co(self):
        return self.id_co[self.c3]

    def get_copy(self):
        nu = U(self.c3)
        nu.o1 = self.o1
        nu.h1 = self.h1[:]
        for k, v in self.id_co.items():
            nu.id_co[k] = v.copy()
        return nu


def get_u_list(a3d):
    _u_list = list()
    for atom_id, atom in a3d.Atoms.items():
        if atom.type == "c3":
            u = U(atom_id)
            u.id_co[atom_id] = atom.coordinate.copy()
            connects = a3d.get_connected_atom_id(atom_id)
            for a_id in connects:
                a = a3d.Atoms[a_id]
                if a.type == "o1":
                    u.o1 = a_id
                    u.id_co[a_id] = a.coordinate.copy()
                elif a.type == "n1":
                    con = a3d.get_connected_atom_id(a_id)
                    for n_id in con:
                        nc = a3d.Atoms[n_id]
                        if nc.type == "h1":
                            u.h1.append(n_id)
                            u.id_co[n_id] = nc.coordinate.copy()
            _u_list.append(u)
    return _u_list


def get_extended_u_in_range(a3d, c_u_list, extend_range):
    def get_moved_u(_u, x, y, z):
        new_u = _u.get_copy()
        d = np.asarray([x, y, z], dtype=int)
        for k in new_u.id_co.keys():
            new_u.id_co[k] -= d * a3d.lattice_parameter
        return new_u

    extended_list = list()
    upper = a3d.lattice_parameter + a3d.box_h - extend_range
    lower = a3d.box_l + extend_range
    for u in c_u_list:
        record = (u.c3_co > upper).astype(int) + (u.c3_co <= lower).astype(int) * -1
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


def cal_h_bond(u1, u2, h_cut=2.5):
    h = 0
    for h_id in u1.h1:
        d = cso.calculate_distance(u1.id_co[h_id], u2.id_co[u2.o1])
        if d <= h_cut:
            h += 1
    for h_id in u2.h1:
        d = cso.calculate_distance(u2.id_co[h_id], u1.id_co[u1.o1])
        if d <= h_cut:
            h += 1
    return h


if __name__ == '__main__':
    fp = "/home/centos/work/fgh_7blk_200chain/1E-7/"
    atom3d = fps.LmpParser.load_data_file(fp + "strain.data")
    h1_list = list()
    h2_list = list()
    h_sum_list = list()
    strain = list()
    for i in range(0, 9500000, 500000):
        atom3d.renew_coordinate_file(fp + "strain.lammpstrj.s1.%d"%i)
        atom3d.convert_to_in_cell()
        cen_u_list = get_u_list(atom3d)
        ext_u_list = get_extended_u_in_range(atom3d, cen_u_list, 10.)
        # for u in ext_u_list:
        #     print(u.id_co)
        cu_co_list = list(map(lambda u:u.c3_co, cen_u_list))
        eu_co_list = list(map(lambda u:u.c3_co, ext_u_list))
        all_u_list = cen_u_list + ext_u_list
        au_co_list = cu_co_list + eu_co_list
        tree = cKDTree(au_co_list)
        cut_off = 8.
        h1 = 0
        h2 = 0
        for cu_id, cu in enumerate(cen_u_list):
            idx_list = tree.query_ball_point(cu.c3_co, cut_off)
            for idx in idx_list:
                if idx != cu_id:
                    h_num = cal_h_bond(cu, all_u_list[idx], h_cut=2.5)
                    if h_num == 1:
                        h1 += 1
                    elif h_num == 2:
                        h2 += 1
                    elif h_num > 2:
                        raise ValueError("More than 2 h-bond!")
            # print(h1, h2)
        strain.append(i * 2E-7)
        h1_list.append(h1/2)
        h2_list.append(h2/2)
        h_sum_list.append(h1/2 + h2)
        print(i)

    plt.plot(strain, h1_list, label="h1")
    plt.plot(strain, h2_list, label="h2")
    plt.plot(strain, h_sum_list, label="h_sum")
    plt.legend()
    plt.show()
