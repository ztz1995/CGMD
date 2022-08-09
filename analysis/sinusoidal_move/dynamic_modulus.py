import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pyMD import functions as f
from scipy.optimize import curve_fit


def load_data(fp):
    data = pd.read_csv(fp, skiprows=1, header=None, sep=" ",
                       names=["strain", "stress_x", "stress_y", "stress_z", "l_x", "l_y", "l_z", "temp"])
    return data


def draw(data):
    plt.plot(data.index[:100000:100], data.strain[:100000:100], label="strain")
    plt.plot(data.index[:100000:100], data.stress_x[:100000:100], label="stress")
    plt.legend()
    plt.tight_layout()
    plt.show()


def cal_modulus(data, A_strain, _tp, dump=10):
    def sin_stress(x, A_stress, fi):
        return f.sin_shift(x, A_stress, 2 * np.pi / _tp * 1000, fi)

    samples = 10
    begin_idx = -_tp * samples // dump
    x = np.asarray(data.index[begin_idx:]) / 1000 * dump  # use ps unit
    p, _ = curve_fit(sin_stress, x, data.stress_x[begin_idx:], p0=[data.stress_x.max(), 0.])
    draw_begin = int(-_tp * samples / 2 // dump)
    plt.plot(x[draw_begin:], data.strain[draw_begin:])
    plt.plot(x[draw_begin:], data.stress_x[draw_begin:])
    plt.plot(x[draw_begin:], sin_stress(x[draw_begin:], *p))
    plt.show()

    return p


if __name__ == '__main__':
    # Tp = 1000000
    # Tp = 100000
    # Tp = 10000
    # Tp = 1000
    Tp = 10000
    # dir_path1 = "/home/centos/work/1blk_50/cg_7blk_300chain/sin_test/"  # cg, sin test
    # dir_path1 = "/home/centos/work/1blk_50/cg_7blk_300chain/sin_test_2/"  # cg, sin test
    # dir_path1 = "/home/centos/work/1blk_50/cg_7blk_300chain/sin_test_100/"  # cg, sin test
    # dir_path1 = "/home/centos/work/1blk_50/cg_7blk_300chain/sin_test_1/"  # cg, sin test
    # file_path_1 = dir_path1 + "strain.%d.txt" % Tp

    dir_path1 = "/home/centos/work/1blk_50/cgu_7blk_300chain/sin_%d/" % Tp  # cg, sin test
    file_path_1 = dir_path1 + "strain.s1.txt"
    data1 = load_data(file_path_1)
    # p = cal_modulus(data1, 0.03, Tp)
    p = cal_modulus(data1, 0.03, Tp, dump=10)
    # p = cal_modulus(data1, 0.1, Tp, dump=100)
    print(p)
    print(np.sin(p[1]) * p[0])
    print(np.cos(p[1]) * p[0])

    # 1E5
    # [0.00964584 0.53212377]
    # 0.004893955818370446
    # 0.008312120649259139

    # 1E4
    # [0.01094599 0.44727178]
    # 0.004734222166896639
    # 0.0098692410300611
