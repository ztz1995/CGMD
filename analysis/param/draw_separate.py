import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
import os
import time

if __name__ == '__main__':
    plt.style.use(["science", "ieee", "no-latex"])
    plt.rc("font", family='Times New Roman')
    # plt.style.use(["science", "ieee"])

    # label = "cgu2cg"
    label = "cg"
    path = "data/%s/" % label
    aa_path = "data/aa/"
    data_dict = dict()
    for i in range(8):
        cg_data = pd.read_csv(path + "%d/non_bond_U_U_%s.csv" % (i, label))
        plt.plot(cg_data.x, cg_data.gr)
    plt.show()
