import numpy as np
from pyMD import functions as func
import matplotlib.pyplot as plt
if __name__ == '__main__':

    # fp = "/home/centos/ztz/solvent/test2_3/non_bond_h1_sv1.param"
    # fp = "/home/centos/ztz/solvent/test2_3/non_bond_o1_sv1.param"
    # fp = "/home/centos/ztz/solvent/test2_3/non_bond_c1_sv1.param"
    # fp = "/home/centos/ztz/solvent/test2_3/non_bond_n1_sv1.param"

    # fp = "/home/centos/ztz/solvent/test2_3/non_bond_c1_sv2.param"
    # fp = "/home/centos/ztz/solvent/test2_3/non_bond_n1_sv2.param"
    # fp = "/home/centos/ztz/solvent/test2_3/non_bond_o1_sv2.param"
    # fp = "/home/centos/ztz/solvent/test2_3/non_bond_h1_sv2.param"

    # fp = "/home/centos/ztz/solvent/test2_2/non_bond_sv2_sv2.param"
    # fp = "/home/centos/ztz/solvent/test2_2/non_bond_TO(2)_sv2.param"
    # fp = "/home/centos/ztz/solvent/test2_2/non_bond_TO(1)_sv2.param"

    # path = "/home/centos/ztz/change_h_bond/8blk_50_sov_8000/"
    # path = "/home/centos/ztz/dpd/cgu/8blk_200/dpd_hybrid/param/"
    path = "/home/centos/ztz/dpd/cgu/8blk_200/1E-6_hb2/param/"
    types = {"h1":7, "o1":9, "n1":8, "c1":6}
    tp = dict()
    for i in types:
        for j in types:
            tp["_".join(sorted([i, j]))] = sorted([types[i], types[j]])

    # for i in range(1, 6):
    for i in range(1):
        for t in tp:

            fp = path + "non_bond_%s.param" % t

            with open(fp, "r") as file:
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

            r = np.asarray(r)
            # e = func.cal_aa_force(r, tuple(t.split("_")), cut_off=14., q=True, dielectric=i, dsf=True)
            e = func.cal_aa_force(r, tuple(t.split("_")), cut_off=14., q=True, dielectric=0.5, dsf=True)
            f = -func.cal_derivative(r, e)
            # # e *= 0.5
            # # f *= 0.5
            # e *= 10
            # f *= 10
            new_lines = lines[:3]
            for j in range(len(index)):
                new_lines.append("%s %s %f %f\n" % (index[j], r[j], e[j], f[j]))

            # fn = fp+"_%d"%i
            fn = fp+"_0.5"
            with open(fn, "w") as file:
                file.writelines(new_lines)
            print("pair_coeff   %d  %d  table 2 %s   %s" % (tp[t][0], tp[t][1], "param/"+fn.split("/")[-1], t))

    # e_list = list()
    # for i in range(1):
    #     for t in tp:
    #
    #         fp = path + "non_bond_%s.param" % t
    #
    #         with open(fp, "r") as file:
    #             lines = file.readlines()
    #
    #         index = list()
    #         r = list()
    #         e = list()
    #         f = list()
    #
    #         for l in lines[3:]:
    #             args = l.split()
    #             if args:
    #                 index.append(args[0])
    #                 r.append(float(args[1]))
    #                 e.append(float(args[2]))
    #                 f.append(float(args[3]))
    #
    #         bidx = 210
    #         # plt.plot(r[bidx:], e[bidx:])
    #         # plt.title(t)
    #         # plt.show()
    #         print(t, e[bidx])
    #         e_list.append(e[bidx])
    #         # print()
    # num_list = [4, 4, 8, 4, 1, 4, 2, 4, 4, 1]
    # esum = 0.
    # for i in range(10):
    #     esum += num_list[i] * e_list[i]
    #
    # e_dict = dict()
    # for i, t in enumerate(tp):
    #     e_dict[t] = 25. / esum * e_list[i]
    # print(e_dict.__repr__())
