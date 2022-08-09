import json
import numpy as np
import matplotlib.pyplot as plt
import os

if __name__ == '__main__':

    # project_name = "20200329_lr=0.1"
    # project_name = "20200329_lr=0.05_SH6S"
    # project_name = "20200408_lr=0.05_hard_aa"
    # project_name = "20200410_lr=0.05_SH6S_fgh"
    # project_name = "20200418_lr=0.05_SH6S_cg_ex_hard"
    # project_name = "20200419_lr=0.05_SH6S_cgu_ex_hard"
    # project_name = "20200422_lr=0.05_1blk_40_cgu_ex_hard"
    # project_name = "20200429_lr=0.05_SH6S_cg_ex_hard"
    # project_name = "20200502_SH6S_cg"
    # project_name = "20200719_lr=0.05_1blk_40_cg_ex_hard"
    # project_name = "20200720_lr=0.05_1blk_40_cg_ex_hard_random_init"
    # project_name = "20200720_lr=0.05_1blk_50_cg_ex_hard"
    # project_name = "20200722_lr=0.05_1blk_50_cgu_ex_hard"
    # project_name = "20200728_lr=0.05_1blk_50_cgu_ex_hard"
    # project_name = "20200914_lr=0.05_1blk_50_600K_cg_ex_hard"
    # project_name = "20200917_lr=0.05_1blk_50_600K_cgu_ex_hard"
    project_name = "20210401_lr=0.05_1blk_50_fg_ex_hard"

    para_name = "non_bond"

    type_tuple = "o1_o1"

    # plot_type = True
    plot_type = False
    start = 1
    max_iter = 500
    error_dict_list = list()
    cm_iter = 0
    for iter_num in range(start, max_iter):
        file_name = "/home/centos/Projects/CGMD/data/%s/record/error_iter_%02d.json" % (project_name, iter_num)
        if os.path.isfile(file_name):
            with open(file_name, "r") as f:
                error_dict_list.append(json.load(f))
            cm_iter = iter_num
        else:
            break
    error_list = list()
    min_error = 1.
    for error_dict in error_dict_list:
        if plot_type:
            error = error_dict[para_name][type_tuple]
        else:
            error = np.array(list(error_dict[para_name].values()))
            # error = np.sqrt(np.square(error).mean())
            error = error.max()
        error_list.append(error)
        min_error = error if error < min_error else min_error

    fig = plt.figure()
    # plt.subplot(2, 1, 1)
    plt.axes([0.08, 0.12, 0.9, 0.8])
    # plt.plot(list(range(cm_iter + 1))[1:], error_list[1:])
    plt.plot(list(range(start, cm_iter + 1)), error_list[:])
    plt.xlim([start, cm_iter + 1])
    plt.ylim(bottom=0.)
    # plt.hlines(min_error, -1, cm_iter + 1, colors="r", linestyles="dashed")
    plt.hlines(0.05, -1, cm_iter + 1, colors="r", linestyles="dashed")
    if plot_type:
        plt.title("%s Error" % type_tuple)
    else:
        plt.title("M. S. Error vs Iteration")
    # ax = fig.add_axes([0, 0.8, 1., 1.0])
    ax = plt.axes([0.06, 0.02, 0.9, 0.05])
    ax.text(0.01, 0.00, project_name, fontdict={'size': 10, 'color': 'black'})
    ax.set_axis_off()
    print("save fig")
    plt.savefig("extension/fit/error.png")
    plt.close()

    parameter_dict_list = list()
    for i in range(start, cm_iter + 1):
        with open("../../data/%s/record/param_iter_%02d_1.json" % (project_name, i), "r") as f:
            parameter_dict_list.append(json.load(f))

    # type_tuple_list = ["B(0)_TO(1)"]
    # type_tuple_list = ["Ph_Ph", "U_U", "Ph_U", "TO(2)_TO(2)", "TO(2)_U","Me_Me"]
    # type_tuple_list = ["U_U"]
    # type_tuple_list = ["Ph_Ph", "Ph_h1", "Ph_o1", "TO(2)_TO(2)", "TO(2)_h1","Me_Me"]
    # type_tuple_list = ["TO(2)_n1", "TO(2)_c1", "TO(2)_o1", "TO(2)_h1","Me_Me"]
    # type_tuple_list = ["TO(1)_n1", "TO(1)_c1", "TO(1)_o1", "TO(1)_h1","c3_c3"]
    type_tuple_list = ["c1_c2", "c2_c2", "c2_c3", "c2_h1", "c2_o1", "c2_n1"]
    # type_tuple_list = ["Ph_Ph"]
    # type_tuple_list = ["Ph_h1"]
    # type_tuple_list = ["c1_c1", "TO(1)_c1", "TO(1)_n1", "TO(1)_o1", "Es_Es"]
    # type_tuple_list = ["TO(1)_TO(1)", "B(0)_TO(1)", "TO(1)_U(0)", "Ph(0)_TO(1)", "Es(0)_U(0)"]
    # type_tuple_list = ["Ph(0)_TO(1)", "B(0)_Ph(0)", "Ph(0)_U(0)", "Ph(0)_Ph(0)"]
    for type_tuple in type_tuple_list:
        # for type_tuple in parameter_dict_list[0][para_name]:
        param_list = list()
        for param_dict in parameter_dict_list:
            param_list.append(np.array(param_dict[para_name][type_tuple]).min())
        plt.plot(list(range(start, cm_iter + 1)), param_list)
        plt.title(type_tuple)
        plt.xlim([start-1., cm_iter + 1])
        plt.show()
