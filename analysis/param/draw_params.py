from pyMD import ibi
import matplotlib.pyplot as plt


if __name__ == '__main__':
    plt.style.use(["science", "ieee", "no-latex"])
    plt.rc("font", family='Times New Roman')

    # fp0 = "../../data/20200719_lr=0.05_1blk_40_cg_ex_hard/record/"
    # r0 = ibi.Record(fp0, restart=True)
    # param0 = r0.get_param(0, 0)

    fp1 = "../../data/20200720_lr=0.05_1blk_50_cg_ex_hard/record/"
    r1 = ibi.Record(fp1, restart=True)
    param1 = r1.get_param(52, 1)
    # fig, ax = plt.subplots(6, 4, dpi=100, figsize=(8, 9))

    # fp1 = "../../data/20200722_lr=0.05_1blk_50_cgu_ex_hard/record/"
    # r1 = ibi.Record(fp1, restart=True)
    # param1 = r1.get_param(353, 1)
    # fig, ax = plt.subplots(9, 4, dpi=1000, figsize=(8, 13.5))
    #
    # draw_x = ibi.IBI.cal_points["non_bond"]
    # count = 0
    # for k, v in param1["non_bond"].items():
    #     if k[0][0].islower() and k[1][0].islower():
    #         continue
    #     title = "-".join(k)
    #     x, y = count//4, count%4
    #     # ax[x][y].plot(draw_x, param0["non_bond"][k])
    #     ax[x][y].plot(draw_x, v, "r")
    #     ax[x][y].set_title(title)
    #     count += 1
    #     ax[x][y].set_ylim(-10, 10)
    #     ax[x][y].set_xlim(0, 15)
    # plt.show()
    for k, v in param1.items():
        if k != "non_bond":
            print(k)
            for tt, pa in v.items():
                if len(pa)>0:
                    print("-".join(tt), end="\t")
                    for value in pa:
                        print("%.2f" % value, end="\t")
                    print()

