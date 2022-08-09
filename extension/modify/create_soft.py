import numpy as np
import os
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


def write_ff(path, index, r, e, f):
    with open(path, "r") as file:
        lines = file.readlines()
    new_lines = lines[:3]
    for j in range(len(index)):
        new_lines.append("%s %.20f %f %f\n" % (index[j], r[j], e[j], f[j]))
    with open(path, "w") as file:
        file.writelines(new_lines)



def harmonic(x, *parameters):
    k, r, d = parameters
    return k * np.square(x - r) + d


def gradual(length, l, r):
    x1 = np.ones(l)
    x2 = np.linspace(1, 0, r-l)
    x3 = np.zeros(length-r)
    x = np.append(np.append(x1, x2), x3)
    y = np.ones(length) - x
    return x, y


def soften(index, r, e, f):
    im = e.argmin()
    print(im)
    p0 = [1, r[im], e[im]]
    p, _ = curve_fit(harmonic, r[im - 300:im-100], e[im - 300:im-100], p0=p0, maxfev=int(1E5))
    print(p)
    # plt.plot(r[im-50:im], harmonic(r[im-50:im], *p))


    g1, g2 = gradual(len(r), im - 300, im-100)
    enew = harmonic(r, *p) * g1 + e * g2
    forc_new = -func.cal_derivative(r, enew)

    plt.plot(r[im-300:], enew[im-300:])
    plt.plot(r[im-300:], e[im-300:])
    plt.show()

    return index, r, enew, forc_new



if __name__ == '__main__':

    fp = "/home/centos/ztz/change_h_bond/8blk_200_soft/"

    files = os.listdir(fp)

    gradual(100, 40, 50)


    for fn in files:
        if fn[:4] == "non_":
            tt = fn.split(".")[0].split("_")[-2:]
            if tt[0][0].isupper() and tt[1][0].isupper():
            # if tt[0][0].isupper() and tt[1][0].islower() or tt[0][0].islower() and tt[1][0].isupper():
            # if tt[0][:2] == "TO" or tt[1][:2] == "TO":
            #     if tt[0][0].islower() or tt[1][0].islower(): continue
                # ind, dist, eng, forc = load_ff(fp + fn)
                print(fn)
                ind, dist, eng, forc = soften(*load_ff(fp + fn))


                dist[0] = 1E-20


                write_ff(fp + fn, ind, dist, eng, forc)



