import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from matplotlib import cm


# plt.style.use(["science", "ieee"])


def load_data(fp):
    data = pd.read_csv(fp, sep="\s+", skiprows=1, header=None)
    # data.columns = ["step", "stress_x", "stress_y", "stress_z", "l_x", "l_y", "l_z", "temp"]
    # data.columns = ["step", "stress_x", "stress_y", "stress_z", "ke_hard", "ke_soft", "virial_hard", "virial_soft", "pair_hard", "pair_soft"]
    data.columns = ["step", "stress_x", "stress_y", "stress_z", "ke_hard_x", "ke_soft_x", "virial_hard_x", "virial_soft_x", "pair_hard_x",
                    "pair_soft_x", "ke_hard_y", "ke_soft_y", "virial_hard_y", "virial_soft_y", "pair_hard_y", "pair_soft_y", "ke_hard_z", "ke_soft_z",
                    "virial_hard_z", "virial_soft_z", "pair_hard_z", "pair_soft_z"]
    return data


if __name__ == '__main__':

    file_path1 = "/home/centos/work/1blk_3200/cgu_8blk_400chain_new/no_relax/8blk_400.relax.txt"
    file_path2 = "/home/centos/work/1blk_3200/cg_8blk_400chain_new/no_relax/8blk_400.relax.txt"

    # file_path1 = "/home/centos/work/1blk_3200/cgu_8blk_400chain/no_relax/8blk_400.relax.txt"
    # file_path2 = "/home/centos/work/1blk_3200/cg_8blk_400chain/no_relax/8blk_400.relax.txt"

    data1 = load_data(file_path2)
    title = "cg"
    # print(data1)
    value_dict = {"%s_hard" % title: dict(), "%s_soft" % title: dict(), "%s_x"%title:0., "%s_y"%title:0., "%s_z"%title:0.}
    for i in ["x", "y", "z"]:
        hard = data1["ke_hard_%s" % i] + data1["virial_hard_%s" % i]
        soft = data1["ke_soft_%s" % i] + data1["virial_soft_%s" % i]
        value_dict["%s_hard" % title][i] = hard.mean()
        value_dict["%s_soft" % title][i] = soft.mean()
        value_dict["%s_%s"%(title,i)] = data1["stress_%s" % i].mean()

    print(value_dict)

    # {'cgu_hard': {'x': 50.765709531785014, 'y': 49.135818963697254, 'z': 47.020965637464606}, 'cgu_soft': {'x': -52.473627764230905, 'y': -50.60273012692953, 'z': -49.99395140978971}}
    # {'cg_hard': {'x': 51.74498946434563, 'y': 51.506702981299746, 'z': 50.08626441759919}, 'cg_soft': {'x': -48.62960914596994, 'y': -47.78661370375951, 'z': -49.4993237109538}}

    # {'cg_hard': {'x': 48.291250050677114, 'y': 49.81049457523296, 'z': 48.366125924518556}, 'cg_soft': {'x': -53.85276241753623, 'y': -54.353069423872206, 'z': -54.63832790574639}}
    # {'cgu_hard': {'x': 46.71195235304606, 'y': 46.47041905655383, 'z': 45.69100164015891}, 'cgu_soft': {'x': -58.04333093565993, 'y': -57.27057988524728, 'z': -56.01408430923586}}
