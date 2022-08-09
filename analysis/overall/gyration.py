import matplotlib.pyplot as plt
import pickle as pkl


def get_gyration(fp, chains):
    with open(fp, "r") as file:
        lines = file.readlines()
    frames = len(lines) // (chains + 1)
    gyrations = list()
    for i in range(frames):
        for j in range(chains):
            g = lines[3 + i * (chains + 1) + j + 1].split()[-1]
            gyrations.append(float(g))
            # print(g)
    return gyrations


if __name__ == '__main__':

    # # AA
    # gyration_list = list()
    # for i in [0, 1, 2, 3, 5, 7, 8, 9]:
    #     fp = "/home/centos/model/aa/1blk_50/%d/gyration.dat" % i
    #     gyration_list += get_gyration(fp, 50)
    # plt.hist(gyration_list, bins=40, range=[min(gyration_list), max(gyration_list)])
    # plt.show()
    # with open("data/gyration_aa.pkl", "wb") as file:
    #     pkl.dump(gyration_list, file)

    # # CG
    # fp = "/home/centos/work/1blk_50/cg_1blk_160chain/equi/gyration.dat"
    # gyration_list = get_gyration(fp, 160)
    # plt.hist(gyration_list, bins=40, range=[min(gyration_list), max(gyration_list)])
    # plt.show()
    # with open("data/gyration_cg.pkl", "wb") as file:
    #     pkl.dump(gyration_list, file)

    # CGU
    fp = "/home/centos/work/1blk_50/cgu_1blk_160chain/equi/gyration.dat"
    gyration_list = get_gyration(fp, 160)
    plt.hist(gyration_list, bins=40, range=[min(gyration_list), max(gyration_list)])
    plt.show()
    with open("data/gyration_cgu.pkl", "wb") as file:
        pkl.dump(gyration_list, file)