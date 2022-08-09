import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':

    # Step Temp Pxx Pyy Pzz Lx Ly Lz E_pair E_bond E_angle E_dihed
    # fp = "/home/centos/work/200-125_1E-7_new/strain.log"
    # fp = "/home/centos/work/7blk_2000chain/200-125_1E-7/strain.log"
    # fp = "/home/centos/work/7blk_500chain/1E-6_fene/strain.log"
    # fp = "/home/centos/work/7blk_500chain/5E-7_df/strain.log"
    # fp = "/home/centos/work/fgh_7blk_200chain/1E-7/strain.log"
    # fp = "/home/centos/work/cg_7blk_200chain/1E-6/strain.log"
    # fp = "/home/centos/work/cg_7blk_200chain/1E-5/strain.log"
    fp = "/home/centos/work/cgu_7blk_200chain/1E-7/strain.log"
    with open(fp, "r") as file:
        lines = file.readlines()
    begin = 0
    end = len(lines)-1
    flag = 0
    for idx, line in enumerate(lines):
        args = line.split()
        if len(args) > 2:
            if flag == 0:
                if args[1] == "first":
                    flag = 1
            if flag == 1:
                if args[0] == "Step":
                    begin = idx
                elif args[0] == "Loop":
                    end = idx
                    flag = 2
    print(begin, end)
    data = pd.read_csv(fp, sep="\\s+", skiprows=begin, nrows=end-begin-1)
    # print(data.head())

    # data = pd.read_csv("relax.csv", sep="\\s+")
    # data = pd.read_csv("strain.csv", sep="\\s+")
    strain = np.asarray(data.Step) / 1E7
    # strain = np.asarray(data.Step) / 2E6
    # strain = np.asarray(data.Step) * 1E-6 * 5
    # strain = np.asarray(data.Step) * 5E-7 * 5
    # k = 2.5
    k = 1
    av = 10
    # plt.plot(data.Temp, data.E_pair)
    # plt.show()
    plt.plot(strain, (data.E_pair - data.E_pair[:av].mean())/k, label="E_pair")
    plt.plot(strain, (data.E_bond - data.E_bond[:av].mean())/k, label="E_bond")
    plt.plot(strain, (data.E_angle - data.E_angle[:av].mean())/k, label="E_angle")
    plt.plot(strain, (data.E_dihed - data.E_dihed[:av].mean())/k, label="E_dihed")
    # plt.plot(strain, data.E_dihed - data.E_dihed[:10].mean(), label="E_dihed")
    plt.legend()
    plt.show()

