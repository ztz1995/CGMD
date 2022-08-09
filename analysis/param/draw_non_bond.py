import matplotlib.pyplot as plt
import pandas as pd
import matplotlib

if __name__ == '__main__':
    plt.style.use(["science", "ieee", "no-latex"])
    plt.rc("font", family='Times New Roman')
    # plt.style.use(["science", "ieee"])

    labels = ["AA", "CGU-initial", "CGU-final"]
    plt.figure(figsize=(4, 3), dpi=1000)
    type_tuple = ("h1", "o1")
    # type_tuple = ("o1", "o1")
    # fps = ["/home/centos/Projects/CGMD/data/target/1blk_50_cgu_ex_hard/aver/non_bond_%s_target.csv" % "_".join(type_tuple),
    #        "/home/centos/Projects/CGMD/data/20200722_lr=0.05_1blk_50_cgu_ex_hard/iter_00/csv/non_bond_%s.csv" % "_".join(type_tuple),
    #        "/home/centos/Projects/CGMD/data/20200722_lr=0.05_1blk_50_cgu_ex_hard/iter_353/csv/non_bond_%s.csv" % "_".join(type_tuple)]

    fps = ["/home/centos/Projects/CGMD/data/target/1blk_50_cgu_ex_hard/aver/non_bond_%s_target.csv" % "_".join(type_tuple),
           "/home/centos/Projects/CGMD/data/20200722_lr=0.05_1blk_50_cgu_ex_hard/iter_00/csv/non_bond_%s.csv" % "_".join(type_tuple),
           "data/10ns/cgu/non_bond_%s_cgu.csv" % "_".join(type_tuple)]
    for i in range(3):
        data = pd.read_csv(fps[i])
        x = data.x
        y = data.gr

        # y /= y.sum() * (x[1] - x[0])
        # x = x[::5]
        # y = y[::5]
        plt.plot(x, y, label=labels[i])
    plt.ylabel("g(r)", fontsize=14)
    plt.xlabel("Distance ($\mathrm{\AA}$)", fontsize=14)
    plt.legend(fontsize=12)
    plt.show()