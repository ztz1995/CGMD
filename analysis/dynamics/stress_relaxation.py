import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

if __name__ == '__main__':
    # plt.figure(figsize=(5, 3))
    # plt.axes(xscale="log", yscale="log")
    # # plt.xlabel("log(t/s)")
    # plt.xlabel("t/s")
    # # plt.ylabel("log(G/G$^0$)")
    # plt.ylabel("G/G$^0$")
    #
    # fp2 = "cgu_1blk_40.csv"
    # data2 = pd.read_csv(fp2, sep="\\s+", skiprows=0)
    # v2 = (data2.v_pxy + data2.v_pxz + data2.v_pyz) / 3
    # t2 = np.asarray(data2.Time)
    # v2 = np.asarray(v2/v2[0])
    # neg2 = np.where(v2<0)
    # v2 = np.delete(v2, neg2)
    # t2 = np.delete(t2, neg2)
    # l = 45
    # r2 = np.argmin(v2)
    # # v -= v.min()-1E-10
    # # print(v)
    # # plt.axes()
    # # plt.axes(xscale="log", yscale="log")
    # plt.plot(t2[l:r2]*1E-15, v2[l:r2], label="cgu")
    #
    #
    # # fp = "cgu_1blk_40.csv"
    # fp = "cg_1blk_40_2.csv"
    # data = pd.read_csv(fp, sep="\\s+", skiprows=0)
    # v = (data.v_pxy + data.v_pxz + data.v_pyz) / 3
    # t = np.asarray(data.Time)
    # v = np.asarray(v/v[0])
    # neg = np.where(v<0)
    # v = np.delete(v, neg)
    # l1 = 50
    # r1 = np.argmin(v)
    # t = np.delete(t, neg)
    # # v -= v.min()-1E-10
    # # print(v)
    # # plt.axes()
    #
    # plt.plot(t[l:r1]*1E-15, v[l:r1], label="cg")
    #
    # plt.legend()
    # plt.show()

    fp1 = "/home/centos/work/1blk_50/cg_1blk_50chain/long_npt/S0St.dat"
    # fp = "/home/centos/work/1blk_50/cg_1blk_50chain/long_nvt/S0St.dat"
    # data = pd.read_csv(fp, sep="\s+", skiprows=3486, header=None)
    fp2 = "/home/centos/work/1blk_50/cgu_1blk_50chain/long_npt_1/S0St_long.dat.bak"
    fp3 = "/home/centos/model/aa/1blk_40/long/1/S0St.dat"

    fps = [fp1, fp2, fp3]
    skips = [1101, 3585049, 391]
    label = ["cg", "cgu", "aa"]
    plt.axes(xscale="log", yscale="log")

    for i in range(3):
        fp = fps[i]
        skip = skips[i]
        data = pd.read_csv(fp, sep="\s+", skiprows=skip, header=None)
        data.columns = ["time","pxy", "pxz", "pyz"]
        v = (data.pxy + data.pxz + data.pyz) / 3
        t = np.asarray(data.time)
        v = np.asarray(v/v[0])
        neg = np.where(v<0)
        v = np.delete(v, neg)
        l = 50
        r1 = np.argmin(v)
        t = np.delete(t, neg)
        print(v)
        # plt.xlabel("log(t/s)")
        # plt.xlabel("t/s")
        # plt.ylabel("log(G/G$^0$)")
        # plt.ylabel("G/G$^0$")

        # plt.plot(t[l:r1]*1E-15, v[l:r1], label="cg")
        plt.plot(t[l:r1], v[l:r1], label=label[i])


    plt.legend()
    plt.show()
