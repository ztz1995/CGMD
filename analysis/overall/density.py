import pandas as pd
from pyMD import analysis as ans
import matplotlib.pyplot as plt
import numpy as np


def properties(_v, _blk, _chain):
    print("a = %f" % np.power(v, 1 / 3))
    print("density = %.4f" % (1515.98 * _blk * _chain / 6.02 / v * 10))


if __name__ == '__main__':
    # # CG
    # fp = "/home/centos/work/1blk_50/cg_1blk_50chain/long_npt/nohup.out"
    # data = pd.read_csv(fp, sep="\s+", skiprows=295, nrows=1000000)
    # v = data.Volume.mean()
    # print("V = %.4f" % v)
    # properties(v)
    # V = 117510.2091
    # a = 48.980724
    # density = 1.0715

    # # CGU
    # fp = "/home/centos/work/1blk_50/cgu_1blk_50chain/long_npt_1/nohup.out"
    # data = pd.read_csv(fp, sep="\s+", skiprows=193)
    # print(data)
    # v = data.Volume.mean()
    # print("V = %.4f" % v)
    # properties(v)
    # V = 117426.0284
    # a = 48.969025
    # density = 1.0723

    # # AA
    # v = 117401.955
    # print("V = %.4f" % v)
    # properties(v)
    # V = 117401.9550
    # a = 48.965678
    # density = 1.0725

    # # CG 1blk_3200
    # title = "cg"
    # blk = 1
    # chain = 3200
    # start = 0
    # fp = "/home/centos/work/1blk_3200/%s_%dblk_%dchain/init/output.17043438" % (title, blk, chain)
    # data = pd.read_csv(fp, sep="\s+", skiprows=99, nrows=10000)
    # print(data)
    # v = data.Volume[start].mean()
    # print("V = %.4f" % v)
    # properties(v, blk, chain)
    # V = 7529197.0000
    # a = 195.997058
    # density = 1.0703

    # # CG 8blk_400
    # title = "cg"
    # blk = 8
    # chain = 400
    # start = 2000
    # fp = "/home/centos/work/1blk_3200/%s_%dblk_%dchain/init/8blk_400.log" % (title, blk, chain)
    # data = pd.read_csv(fp, sep="\s+", skiprows=132, nrows=10000)
    # print(data)
    # v = data.Volume[start].mean()
    # print("V = %.4f" % v)
    # properties(v, blk, chain)



    # # CGU 8blk_400
    # title = "cgu"
    # blk = 8
    # chain = 400
    # start = 1000
    # fp = "/home/centos/work/1blk_3200/%s_%dblk_%dchain/init/output.17083589" % (title, blk, chain)
    # data = pd.read_csv(fp, sep="\s+", skiprows=115, nrows=10000)
    # print(data)
    # v = data.Volume[start].mean()
    # print("V = %.4f" % v)
    # properties(v, blk, chain)
    # V = 7242587.2000
    # a = 193.477846
    # density = 1.1126

    # # CGU 1blk_160
    # title = "cgu"
    # blk = 1
    # chain = 160
    # start = 2000
    # fp = "/home/centos/work/1blk_160/%s_%dblk_%dchain/equi/output.17042366" % (title, blk, chain)
    # data = pd.read_csv(fp, sep="\s+", skiprows=272, nrows=10000)
    # print(data)
    # v = data.Volume[start].mean()
    # print("V = %.4f" % v)
    # properties(v, blk, chain)
    # V = 376286.6200
    # a = 72.194857
    # density = 1.0708

    v = 328.84314 * 328.73317 * 327.55211
    properties(v, 8, 2000)