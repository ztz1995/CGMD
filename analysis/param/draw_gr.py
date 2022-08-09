import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
import os
import time
from pyMD.functions import error_l1
import numpy as np

if __name__ == '__main__':
    plt.style.use(["science", "ieee", "no-latex"])
    plt.rc("font", family='Times New Roman')
    # plt.style.use(["science", "ieee"])

    # # label = "cg"
    # fig, ax = plt.subplots(6, 4, dpi=1000, figsize=(8, 9))
    # cg_path = "data/cg/"
    # cgu_path = "data/cgu2cg/"
    # aa_path = "data/aa2cg/"
    # files = os.listdir(cg_path)
    # data_dict = dict()
    # count = 0
    # for fn in files:
    #     if fn[-3:] == "csv":
    #         if fn[:5] == "non_b":
    #         # if fn[:5] != "non_b":
    #             cg_data = pd.read_csv(cg_path + fn)
    #             aa_name = fn[:-4-len("cg")] + "aa2cg.csv"
    #             aa_data = pd.read_csv(aa_path + aa_name)
    #             cgu_name = fn[:-4-len("cg")] + "cgu2cg.csv"
    #             cgu_data = pd.read_csv(cgu_path + cgu_name)
    #             x = count // 4
    #             y = count % 4
    #             ax[x][y].plot(aa_data.x, aa_data.gr, label="AA")
    #             ax[x][y].plot(cgu_data.x, cgu_data.gr, label="CGU")
    #             ax[x][y].plot(cg_data.x, cg_data.gr, label="CG")
    #             fn = fn.split("_")
    #             # title = "-".join(fn[1:-1])
    #             title = "-".join(fn[2:4])
    #             print(title)
    #             ax[x][y].set_title(title, fontsize=10)
    #             # if fn[0] == "Angle":
    #             #     xlabel = "Bond Angle (degrees)"
    #             # elif fn[0] == "Bond":
    #             #     xlabel = "Bond Length ($\mathrm{\AA}$)"
    #             # else:
    #             #     xlabel = "Dihedral Angle (degrees)"
    #             xlabel = "Distance ($\mathrm{\AA}$)"
    #             ax[x][y].set_xlabel(xlabel, fontsize=10)
    #             count += 1
    # plt.show()

    # label = "cgu"
    # fig, ax = plt.subplots(6, 4, dpi=1000, figsize=(8, 9))
    # # fig, ax = plt.subplots(6, 4, dpi=1000, figsize=(8, 9))
    #
    # cgu_path = "data/cgu/"
    # aa_path = "data/aa2cgu/"
    # files = os.listdir(cgu_path)
    # data_dict = dict()
    # count = 0
    # for fn in files:
    #     if fn[-3:] == "csv":
    #         # if fn[:5] == "non_b":
    #         if fn[:5] != "non_b":
    #             fnl = fn.split("_")
    #             # title = "-".join(fnl[2:-1])
    #             title = "-".join(fnl[1:-1])
    #             # tl = [_[0].islower() for _ in fnl[2:-1]]
    #             tl = [_[0].islower() for _ in fnl[1:-1]]
    #             if not any(tl):
    #                 continue
    #             cgu_data = pd.read_csv(cgu_path + fn)
    #             aa_name = fn[:-4-len("cgu")] + "aa2cgu.csv"
    #             aa_data = pd.read_csv(aa_path + aa_name)
    #             x = count // 4
    #             y = count % 4
    #             ax[x][y].plot(aa_data.x, aa_data.gr, label="AA")
    #             ax[x][y].plot(cgu_data.x, cgu_data.gr, label="CGU")
    #             ax[x][y].set_title(title, fontsize=10)
    #             if fn[0] == "Angle":
    #                 xlabel = "Bond Angle (degrees)"
    #             elif fn[0] == "Bond":
    #                 xlabel = "Bond Length ($\mathrm{\AA}$)"
    #             else:
    #                 xlabel = "Dihedral Angle (degrees)"
    #             # xlabel = "Distance ($\mathrm{\AA}$)"
    #             ax[x][y].set_xlabel(xlabel, fontsize=10)
    #             print(title)
    #             count += 1
    # plt.show()



    label = "cgu"
    # fig, ax = plt.subplots(8, 4, dpi=100, figsize=(8, 12))
    exclude = []
    cgu_path = "data/10ns/cgu/"
    aa_path = "data/10ns/aa2cgu/"
    files = os.listdir(aa_path)
    print(files)
    count = 0
    # _type = 0
    _type = 1
    if _type == 0:
        fig, ax = plt.subplots(8, 4, dpi=1000, figsize=(8, 12))
    else:
        fig, ax = plt.subplots(6, 4, dpi=1000, figsize=(8, 9))
    for fn in files:
        if fn[-3:] == "csv":
            if _type == 0:
                if fn[:5] != "non_b":
                    continue
            else:
                if fn[:5] == "non_b":
                    continue
            fnl = fn.split("_")
            start_idx = 2 if _type == 0 else 1
            title = "-".join(fnl[start_idx:-1])
            tl = [_[0].islower() for _ in fnl[start_idx:-1]]
            if not any(tl):
                continue
            print(title)
            aa_name = fn
            aa_data = pd.read_csv(aa_path + aa_name)
            x = count // 4
            y = count % 4
            ax[x][y].plot(aa_data.x, aa_data.gr, label="AA")
            ax[x][y].set_title(title, fontsize=10)
            cgu_data = pd.read_csv(cgu_path + fn[:-4 - len("aa2cgu")] + "cgu.csv")
            ax[x][y].plot(cgu_data.x, cgu_data.gr, label="CGU")
            # for i in range(10):
            #     if i in exclude:
            #         continue
            #     cgu_data = pd.read_csv(cgu_path + "%d/" % i + fn[:-4-len("aa2cgu")] + "cgu.csv")
            #     ax[x][y].plot(cgu_data.x, cgu_data.gr, label="CGU_%d" % i)
            if _type == 0:
                xlabel = "Distance ($\mathrm{\AA}$)"
            else:
                if fn[0] == "Angle":
                    xlabel = "Bond Angle (degrees)"
                elif fn[0] == "Bond":
                    xlabel = "Bond Length ($\mathrm{\AA}$)"
                else:
                    xlabel = "Dihedral Angle (degrees)"
            ax[x][y].set_xlabel(xlabel, fontsize=10)

            count += 1
    plt.show()


    # label = "cg"
    # # fig, ax = plt.subplots(8, 4, dpi=100, figsize=(8, 12))
    # cgu_path = "data/10ns/cgu2cg/"
    # cg_path = "data/10ns/cg/"
    # aa_path = "data/10ns/aa2cg/"
    # files = os.listdir(aa_path)
    # print(files)
    # count = 0
    # # _type = 0
    # _type = 1
    # cg_error = list()
    # cgu_error = list()
    # if _type == 0:
    #     fig, ax = plt.subplots(6, 4, dpi=1000, figsize=(8, 9))
    # else:
    #     fig, ax = plt.subplots(6, 4, dpi=1000, figsize=(8, 9))
    # for fn in files:
    #     if fn[-3:] == "csv":
    #         if _type == 0:
    #             if fn[:5] != "non_b":
    #                 continue
    #         else:
    #             if fn[:5] == "non_b":
    #                 continue
    #         fnl = fn.split("_")
    #         start_idx = 2 if _type == 0 else 1
    #         title = "-".join(fnl[start_idx:-1])
    #         print(title)
    #         aa_name = fn
    #         aa_data = pd.read_csv(aa_path + aa_name)
    #         x = count // 4
    #         y = count % 4
    #         ax[x][y].plot(aa_data.x, aa_data.gr, label="AA")
    #         ax[x][y].set_title(title, fontsize=10)
    #         cgu_data = pd.read_csv(cgu_path + fn[:-4 - len("aa2cg")] + "cgu2cg.csv")
    #         ax[x][y].plot(cgu_data.x, cgu_data.gr, label="CGU")
    #         cg_data = pd.read_csv(cg_path + fn[:-4 - len("aa2cg")] + "cg.csv")
    #         ax[x][y].plot(cg_data.x, cg_data.gr, label="CG")
    #         cg_error.append(error_l1(cg_data.gr, aa_data.gr))
    #         cgu_error.append(error_l1(cgu_data.gr, aa_data.gr))
    #         # for i in range(10):
    #         #     if i in exclude:
    #         #         continue
    #         #     cgu_data = pd.read_csv(cgu_path + "%d/" % i + fn[:-4-len("aa2cgu")] + "cgu.csv")
    #         #     ax[x][y].plot(cgu_data.x, cgu_data.gr, label="CGU_%d" % i)
    #         if _type == 0:
    #             xlabel = "Distance ($\mathrm{\AA}$)"
    #         else:
    #             if fn[0] == "Angle":
    #                 xlabel = "Bond Angle (degrees)"
    #             elif fn[0] == "Bond":
    #                 xlabel = "Bond Length ($\mathrm{\AA}$)"
    #             else:
    #                 xlabel = "Dihedral Angle (degrees)"
    #         ax[x][y].set_xlabel(xlabel, fontsize=10)
    #
    #         count += 1
    # plt.show()
    # print(max(cg_error))
    # print(max(cgu_error))
    # print(np.average(cg_error))
    # print(np.average(cgu_error))