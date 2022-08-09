import pandas as pd
from pyMD import analysis as ans
import matplotlib.pyplot as plt
import numpy as np
import tables
import pickle as pkl


if __name__ == '__main__':
    # fp = "/home/centos/work/1blk_50/cg_1blk_50chain/long_npt/nohup.out"
    # data = pd.read_csv(fp, sep="\s+", skiprows=295)

    # fp = "/home/centos/work/1blk_50/cg_1blk_50chain/long_nvt/stress.dat"
    # data = pd.read_csv(fp, sep=",", skiprows=1, header=None)
    # data.columns = ["Step", "Pxy", "Pxz", "Pyz"]
    # data.to_hdf("/home/centos/fast/temp/cg_nvt.h5", 'df')
    data = pd.read_hdf("/home/centos/fast/temp/cg_nvt.h5", 'df')

    # print(data.head())

    frames = 1000000
    steps = 10
    nums = 100000
    # start_delta = 5
    starts = np.random.randint(1, 100000, nums)

    labels = ["Pxy", "Pxz", "Pyz"]
    G = np.zeros(frames)
    for label in labels:
        stress_cal = ans.StressCalculator(np.asarray(data[label]), step=1)
        auto = ans.AutoCorrelation(stress_cal, frames, steps, starts)
        g = auto.cal()
        print(g)
        G += g
    with open("data/G_cg_nvt_%d_%d_%d.pkl" % (frames, steps, nums), "wb") as f:
        pkl.dump(G/3, f)

    t = np.arange(frames) * steps
    v = np.asarray(G/G[0])
    neg = np.where(v<0)
    v = np.delete(v, neg)
    t = np.delete(t, neg)
    plt.plot(t, v)

    # plt.plot(np.arange(frames) * steps, G)
    plt.xscale("log")
    plt.yscale("log")
    plt.show()
    print(G/3)