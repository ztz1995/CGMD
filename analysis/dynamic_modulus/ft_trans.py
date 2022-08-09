import numpy as np
from numpy import fft
import pandas as pd
import matplotlib.pyplot as plt


def time2index(_t):
    _record = {16: 16, 32: 24, 64: 32, 128: 40, 256: 48, 512: 56, 1024: 64, 2048: 72, 4096: 80, 8192: 88, 16384: 96,
               32768: 104, 65536: 112, 131072: 120, 262144: 128, 524288: 136, 1048576: 144, 2097152: 152, 4194304: 160,
               8388608: 168, 16777216: 176, 33554432: 184, 67108864: 192}
    if _t == 0.:
        return 0
    log2t = int(np.log2(_t))
    if log2t < 4:
        return int(_t)
    flag = np.power(2, log2t)
    return int((_t - flag)/(flag/8) + _record[flag])


if __name__ == '__main__':

    # t = 983040
    # t_idx = time2index(t)

    data = pd.read_csv("cgu_Gt.csv", sep=" ")
    data["Gaver"] = (data.Gx + data.Gy + data.Gz) / 3
    # print(data.time[t_idx])

    f_list = list()
    g_list = list()
    for i in range(2, 20):
        d = np.power(2, i)
        time = range(0, 16 * d, d)
        time_idx = list(map(time2index, time))
        # print(data.loc[time_idx, :])
        # print(np.asarray(data.loc[time_idx, ["Gaver"]]).flatten())
        fG = fft.fft(np.asarray(data.loc[time_idx, ["Gaver"]]).flatten())
        freq = fft.fftfreq(16, d=d)
        f_list.append(freq[1])
        g_list.append(fG[1])

        # plt.scatter(freq[2], fG.real[2])
    plt.loglog(f_list, g_list)
    plt.show()

    # t = list(data.time)
    # G = np.asarray((data.Gx + data.Gy + data.Gz) / 3)

    # print(t[:16+1])
    # print(t[:16:2] + t[16:24+1])
    # print(t[:16:4] + t[16:24:2] + t[24:32+1])
    # print(t[:16:8] + t[16:24:4] + t[24:32:2] + t[32:40+1])

    # t = 16
    # l_idx = 16
    # record_dict = dict()
    # for i in range(24):
    #     record_dict[t] = l_idx
    #     t *= 2
    #     l_idx += 8
    # print(record_dict)

    # print(t[:64:4])
    # fG = fft.fft(G[:16])
    # freq = fft.fftfreq(16)
    # plt.scatter(freq.real[1:], fG[1:])
    # print(freq)
    # print(fG)
    # plt.show()

    # print()
    # print(fft.fftfreq(16, d=2))
    # print(fft.fftfreq(16, d=4))
    # print(fft.fftfreq(16, d=8))
