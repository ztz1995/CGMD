import os
import numpy as np
import matplotlib.pyplot as plt
from pyMD import functions as func
from scipy.optimize import curve_fit


def load_ff(path):

    with open(path, "r") as file:
        lines = file.readlines()

    index = list()
    r = list()
    e = list()
    f = list()

    for l in lines[3:]:
        args = l.split()
        if args:
            index.append(args[0])
            r.append(float(args[1]))
            e.append(float(args[2]))
            f.append(float(args[3]))

    r = np.asarray(r,dtype=np.float64)
    e = np.asarray(e,dtype=float)
    f = np.asarray(f,dtype=float)
    return index, r, e, f


def write_ff(path, index, r, e, f, append_name=""):
    with open(path, "r") as file:
        lines = file.readlines()
    new_lines = lines[:3]
    for j in range(len(index)):
        new_lines.append("%s %.20f %f %f\n" % (index[j], r[j], e[j], f[j]))
    with open(path+append_name, "w") as file:
        file.writelines(new_lines)


if __name__ == '__main__':

    # path = "/home/centos/ztz/change_h_bond/8blk_50_BD_2/param/"
    # files = os.listdir(path)
    # for fn in files:
    #     if fn[:4] == "non_":
    #         _index, _r, _e, _f = load_ff(path + fn)
    #         for i in range(1, 11):
    #             ne = _e * i / 10
    #             nf = -func.cal_derivative(_r, ne)
    #             write_ff(path+fn, _index, _r, ne, nf, "_%d" % i)

    path = "/home/centos/ztz/change_h_bond/change_TO/soft/param/"
    files = os.listdir(path)
    for fn in files:
        if fn == "non_bond_TO(2)_TO(2).param":
            _index, _r, _e, _f = load_ff(path + fn)
            from pyMD.functions import pow2
            ne = pow2(_r, 5., 5.)
            nf = -func.cal_derivative(_r, ne)
            write_ff(path+fn, _index, _r, ne, nf, "_soft")