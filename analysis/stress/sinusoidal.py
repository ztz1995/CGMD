import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from matplotlib import cm


def main():

    # plt.figure(figsize=(7,5))
    process = "s1"
    # dir_path1 = "/home/centos/work/1blk_50/cg_7blk_300chain/sin_test/"   # cg, sin test
    dir_path1 = "/home/centos/work/1blk_50/cg_7blk_300chain/sin_test_1/"   # cg, sin test

    # print(colors)
    # file_path_1 = dir_path1 + "strain.%s.txt" % process
    # file_path_1 = dir_path1 + "strain.%d.txt" % 100000
    file_path_1 = dir_path1 + "strain.%d.txt" % 1000000
    data = pd.read_csv(file_path_1, skiprows=1, header=None, sep=" ",
                       names=["strain", "stress_x", "stress_y", "stress_z", "l_x", "l_y", "l_z", "temp"])
    # print(data)
    plt.plot(data.index[:100000:100], data.strain[:100000:100], label="strain")
    plt.plot(data.index[:100000:100], data.stress_x[:100000:100], label="stress")
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    main()


