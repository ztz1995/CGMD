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
    print("finish load")

    frames = 1000000
    steps = 100
    nums = 1
    # start_delta = 5
    starts = np.random.randint(1, 100000, nums)

    labels = ["Pxy", "Pxz", "Pyz"]
    G_sum = np.zeros(frames)
    for i, start in enumerate(starts):
        print(i)
        G = np.zeros(frames)
        for label in labels:
            v = np.asarray(data[label])
            vs = v[start:start + frames * steps:steps]
            g = np.correlate(vs, vs, mode="full")
            g = g[len(g) // 2:]
            g /= np.arange(len(g), 0, -1)
            G += g
        G_sum += G

        # stress_cal = ans.StressCalculator(None, step=1)
        # auto = ans.AutoCorrelation(stress_cal, frames, steps, starts)
        # print("finish 1")
        # g = auto.cal_parallel_shared(list(data[label]), parallel=5)
        # print(g)
        # G += g
    G_sum /= nums
    with open("data/G_cg_nvt_%d_%d_%d.pkl" % (frames, steps, nums), "wb") as f:
        pkl.dump(G_sum, f)

    # with open("data/G_cg_nvt_%d_%d_%d.pkl" % (frames, steps, nums), "rb") as f:
    #     G = pkl.load(f)

    t = np.arange(frames) * steps
    v = np.asarray(G / G[0])

    v_aver = list()
    for i in range(int(len(v) / 1.1)):
        l = round(i * 0.9)
        r = round(i * 1.1)
        v_aver.append(v[l:r + 1].mean())

    t = t[:len(v_aver)]
    v = np.array(v_aver)
    neg = np.where(v < 0)
    v = np.delete(v, neg)
    t = np.delete(t, neg)
    plt.plot(t, v)

    # plt.plot(np.arange(frames) * steps, G)
    plt.xscale("log")
    plt.yscale("log")
    plt.show()
