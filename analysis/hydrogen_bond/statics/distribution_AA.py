from pyMD import analysis as ans
import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl


if __name__ == '__main__':

    # AA
    # for i in [1, 2, 3, 5, 7, 8, 9]:
    # hb = list()
    # for j in [0]:
    #     h_bond_sum = list()
    #
    #     fp = "/home/centos/model/aa/1blk_50/%d/" % j
    #     title = "aa"
    #     fp_data = fp + "1blk_50.data"
    #     fp_trj = fp + "1blk_50.lammpstrj"
    #     fp_info = fp + "cgu_struct_info.pkl"
    #     frames = 1000
    #     atom3ds = ans.Atom3dGenerator(fp_data, fp_trj, one_file=True, info_path=fp + "cgu_struct_info.pkl", pad=0)
    #
    #     start = 15000000
    #     steps = 1000
    #     h_cal = ans.HydrogenBond(atom3ds())
    #
    #     for i in range(frames):
    #         print(i)
    #         atom3d = atom3ds(start + i * steps, only_aa=True)
    #         hbonds = h_cal.distribution(atom3d, 4.0)
    #         h_bond_sum += hbonds
    #
    #     with open("data/distribution_%s_500K_%d_%d.pkl" % (title, frames, j), "wb") as f:
    #         pkl.dump(h_bond_sum, f)
    #     hb += h_bond_sum
    #     h_lens = [c[0] for c in h_bond_sum]
    #     h_angs = [c[1] for c in h_bond_sum]
    #     plt.hist2d(h_lens, h_angs, range=[[1., 4], [90, 180]], bins=[100,100], density=True)
    #     plt.colorbar()
    #     plt.show()

    plt.style.use(["science", "ieee", "no-latex"])
    plt.rc("font", family='Times New Roman')
    font_size = 14

    fig = plt.figure(dpi=1000, figsize=(3.,2.7))
    ax = fig.add_subplot(111)
    ax.xaxis.set_tick_params(top=False,bottom=False,left=False,right=False)
    ax.yaxis.set_tick_params(top=False,bottom=False,left=False,right=False)

    frames = 1000
    title = "aa"
    hb = list()
    for j in [0, 1, 2, 3, 5, 7, 8, 9]:
        with open("data/distribution_%s_500K_%d_%d.pkl" % (title, frames, j), "rb") as f:
            hb_l = pkl.load(f)
        hb += hb_l
    h_lens = [c[0] for c in hb]
    h_angs = [c[1] for c in hb]
    plt.hist2d(h_lens, h_angs, range=[[1., 4], [90, 180]], bins=[100, 100], density=True, vmin=0., vmax=0.064)
    y_major_locator = plt.MultipleLocator(20)
    plt.gca().yaxis.set_major_locator(y_major_locator)
    # plt.yticks(fontname="cmr10", fontsize=font_size)
    # plt.xticks(fontname="cmr10", fontsize=font_size)
    plt.ylabel("Angle (degrees)", fontsize=font_size)
    plt.xlabel("Distance ($\mathrm{\AA}$)", fontsize=font_size)
    plt.show()