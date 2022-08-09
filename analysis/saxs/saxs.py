from pyMD import analysis as ans
import matplotlib.pyplot as plt
import numpy as np
from pyMD.file_parser import LmpParser


def plot_2d(fig_name, _axis, _v, _save_path=None):
    # print(_axis)
    # plt.figure(figsize=(20, 20))
    plt.figure(figsize=(6, 5))
    # plt.imshow(_v, interpolation='bilinear', v max=1.6, vmin=0.6)  # cmap="gray",
    # plt.imshow(_v, interpolation='bilinear', vmax=1000, vmin=0)  # cmap="gray",
    # plt.imshow(_v, interpolation='bilinear', vmax=_v.max(), vmin=0)  # cmap="gray",
    # plt.imshow(_v, interpolation='bilinear', vmax=3E8, vmin=0)  # cmap="gray",
    plt.pcolormesh(_v, vmax=5E7, vmin=0)  # cmap="gray",

    # plt.contourf(np.arange(-rc, rc, 0.1), np.arange(-rc, rc, 0.1), _v, levels=levels)
    plt.colorbar()
    plt.title(fig_name)
    # plt.xticks(np.arange(0, 25, 5), np.arange(-0.5, 0.75, 0.25))
    # plt.yticks(np.arange(0, 25, 5), np.arange(-0.5, 0.75, 0.25))
    # plt.yticks(np.arange(0, 110, 10), np.arange(-2.5, 3.0, 0.5))
    plt.xlabel(_axis[1])
    plt.ylabel(_axis[0])
    if _save_path:
        plt.savefig("%s/%s.png" % (_save_path, fig_name), dpi=300)
    else:
        plt.show()
    plt.close()


if __name__ == '__main__':

    # fp = "/home/centos/work/cgu_7blk_500chain/init/"
    # fp = "/home/centos/work/cgu_7blk_500chain/init/"
    # fp = "/home/centos/work/cgu_2blk_800chain_2/init/"
    # fp = "/home/centos/work/cgu_7blk_200chain/anneal_480/"
    # fp = "/home/centos/work/cgu_7blk_200chain/init/"
    # fp = "/home/centos/work/cg_7blk_200chain/init/"
    # fp = "/home/centos/work/cg_1blk_1000chain/init/"
    # fp = "/home/centos/work/cgu_1blk_1000chain/init/"
    # start = 10000000

    # fp = "/home/centos/work/SH6S/cg_7blk_300chain/init/"
    # start = 8000000
    # start = 200000
    # fp = "/home/centos/work/cgu_7blk_200chain/init/"
    # start = 15400000
    # fp = "/home/centos/work/cgu_7blk_500chain/init/"
    # start = 11000000
    # fp = "/home/centos/work/cgu_1blk_1000chain/init/"
    # start = 10300000

    # fp = "/home/centos/work/SH6S/cgu_7blk_300chain/init/"
    # start = 9500000
    # fp = "/home/centos/work/SH6S/cgu_7blk_300chain/anneal_400/"
    # start = 17200000
    # start = 7200000
    # start = 200000

    # fp = "/home/centos/work/SH6S/cg_7blk_300chain/anneal_400/"
    # start = 11500000

    fp = "/home/centos/work/SH6S/cgu_1blk_2000chain/init/"
    start = 10300000


    fp_data = fp + "strain.data"
    fp_trj = fp + "record.lammpstrj."

    # start = 6900000
    max_len = 120
    # atom3ds = ans.atom3d_generator(fp_data, fp_trj, start=7300000, step=100000, frames=1)
    # atom3ds = ans.atom3d_generator(fp_data, fp_trj, start=3600000, step=100000, frames=1)
    # atom3ds = ans.atom3d_generator(fp_data, fp_trj, start=15400000, step=100000, frames=1)
    # atom3ds = ans.atom3d_generator(fp_data, fp_trj, start=15400000, step=100000, frames=1)
    atom3ds = ans.atom3d_generator(fp_data, fp_trj, start=start, step=100000, frames=1)
    a = next(atom3ds)
    # a = LmpParser.load_data_file(fp_data)
    # a = LmpParser.load_data_file("../pure.data")
    grids = 2 * np.pi / max_len / np.linspace(0.5, 2.5, 25, endpoint=False) * 10  # A -> nm
    # fs = np.asarray([0.06, 0.12, 0.24])
    # grids = 2 * np.pi / max_len / np.linspace(0.08, 0.16, 2, endpoint=True)
    # grids = 2 * np.pi / max_len / fs
    print(grids)
    plt.figure(figsize=(5, 4), dpi=300)
    plt.xlabel("q/nm$^{-1}$")
    freq_list = list()
    ints_list = list()
    ints_list1 = list()
    ints_list2 = list()
    for grid in grids:
        ds = np.zeros((10, 10))
        for i in range(3):
            d = ans.StaticAnalyzer.saxs(a, grid=grid, bin_num=max_len, iter_num=100, drop_axis=i)[:10, :10]
            # plot_2d(str(grid), "yz", d)
            ds += d
        # ds = ans.StaticAnalyzer.saxs(a, grid=grid, bin_num=120, iter_num=200, drop_axis=1)
        # plot_2d(str(grid), "yz", ds/3)
        # freq = np.arange(1, 2) / max_len / grid * 2 * np.pi
        freq = 1. / max_len / grid * 2 * np.pi * 10
        # freq = 2. / max_len / grid * 2 * np.pi
        freqs = np.arange(1,10) / max_len / grid * 2 * np.pi * 10
        # plt.plot(freqs, ds[1:10, 0])
        # plt.show()
        print(freqs)
        print(ds[1:10, 0])
        # plt.show()
        # intensity = ds[1:8, 0]
        intensity = (ds[1, 0] + ds[0, 1])/2/3  * (grids[0]/grid) *(grids[0]/grid)
        intensity1 = (ds[1, 0] + ds[0, 1])/2/3
        # intensity = (ds[2, 0] + ds[0, 2])/2/3
        # intensity1 = ds[1, 0]/3
        # intensity2 = ds[0, 1]/3
        freq_list.append(freq)
        ints_list.append(intensity)
        ints_list1.append(intensity1)
        # ints_list2.append(intensity2)
    # plt.show()
    freqs = np.asarray(freq_list)
    freqs = freqs.flatten(order="F")
    intss = np.asarray(ints_list)
    intss = intss.flatten(order="F")
    print(freqs, intss)
    print(freqs.__repr__(), intss.__repr__())
    plt.plot(freqs, intss)
    plt.show()
    # freqs = np.asarray(freqs)
    # intss = np.asarray(intss)
    # plt.plot(freqs, intss*freqs*freqs)
    # plt.show()
    # # ds = np.fft.fftshift(ds)
    # # plot_2d("test", "yz", ds[:2, :20])
    # plt.plot(freq_list, ints_list1)
    # plt.show()
    # freqs = np.asarray(freq_list)
    # intss = np.asarray(ints_list1)
    # plt.plot(freqs, intss*freqs*freqs)
    # plt.show()

    # plt.plot(freq_list, ints_list2)
    # plt.show()
