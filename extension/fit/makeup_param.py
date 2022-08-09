from pyMD import ibi

if __name__ == '__main__':

    # fp1 = "../data/20200410_lr=0.05_SH6S_fgh/record/"
    # fp2 = "../data/20200408_lr=0.05_hard_aa/record/"
    # # fp2 = "../data/20200329_lr=0.1/record/"
    #
    # r1 = ibi.Record(fp1)
    # param1 = r1.get_param(0, 0)
    #
    # r2 = ibi.Record(fp2)
    # param2 = r2.get_param(80, 1)
    #
    # new_param = dict()
    # for para in param1:
    #     new_param[para] = dict()
    #     for type_tuple in param1[para]:
    #         if type_tuple in param2[para]:
    #             new_param[para][type_tuple] = param2[para][type_tuple]
    #         else:
    #             new_tt = list()
    #             for _type in type_tuple:
    #                 if _type == "TO(3)":
    #                     new_tt.append("TO(1)")
    #                 else:
    #                     new_tt.append(_type)
    #             new_tt = tuple(new_tt)
    #             new_param[para][type_tuple] = param2[para][new_tt]
    #
    #
    #             # new_param[para][type_tuple] = param1[para][type_tuple]
    #         # if para == "non_bond":
    #         #     if type_tuple[0] not in ibi.IBI.aa_type and type_tuple[1] not in ibi.IBI.aa_type:
    #         #         new_param[para][type_tuple] = param1[para][type_tuple]
    #         #     else:
    #         #         new_param[para][type_tuple] = param2[para][type_tuple]
    #         # else:
    #         #     # if type_tuple in param2[para]:
    #         #     #     new_param[para][type_tuple] = param2[para][type_tuple]
    #         #     # else:
    #         #     new_param[para][type_tuple] = param1[para][type_tuple]
    # ibi.Record.dump_dict_json(new_param, fp1+"param_iter_01_1.json")

    # fp1 = "../../data/20200418_lr=0.05_SH6S_cg_ex_hard/record/"
    # fp1 = "../../data/20200418_lr=0.05_SH6S_cg_ex_hard/record/"
    # fp2 = "../../data/20200420_lr=0.05_1blk_40_cg_ex_hard/record/"
    # # fp2 = "../data/20200329_lr=0.1/record/"
    #
    # r1 = ibi.Record(fp1, restart=True)
    # # param1 = r1.get_param(184, 1)
    # param1 = r1.get_param(289, 1)
    #
    # r2 = ibi.Record(fp2, restart=True)
    # param2 = r2.get_param(22, 0)
    #
    # new_param = dict()
    # for para in param2:
    #     new_param[para] = dict()
    #     for type_tuple in param2[para]:
    #         # if para == "non_bond" and ((type_tuple[0][0].islower() and type_tuple[1][0].isupper()) or (type_tuple[1][0].islower() and type_tuple[0][0].isupper())):
    #         #     new_param[para][type_tuple] = param1[para][type_tuple] * 0.2
    #         # else:
    #         #     new_param[para][type_tuple] = param2[para][type_tuple]
    #         if type_tuple in param1[para]:
    #             new_param[para][type_tuple] = param1[para][type_tuple]
    #         else:
    #             new_param[para][type_tuple] = param2[para][type_tuple]
    #         # # elif para == "non_bond" and not (type_tuple[0][0].islower() and type_tuple[1][0].islower()):
    #         # #     if type_tuple[0][0].islower():
    #         # #         new_param[para][type_tuple] = param1[para][("U", type_tuple[1])]
    #         # #     else:
    #         # #         new_param[para][type_tuple] = param1[para][(type_tuple[0], "U")]
    #         # else:
    #         #     # new_tt = list()
    #         #     # for _type in type_tuple:
    #         #     #     if _type == "TO(3)":
    #         #     #         new_tt.append("TO(1)")
    #         #     #     else:
    #         #     #         new_tt.append(_type)
    #         #     # new_tt = tuple(new_tt)
    #         #     if para == "non_bond":
    #         #         new_param[para][type_tuple] = param2[para][type_tuple] * 0.2
    #         #     else:
    #         #         new_param[para][type_tuple] = param2[para][type_tuple]
    #
    # ibi.Record.dump_dict_json(new_param, fp2+"param_iter_23_0.json")

    # fp1 = "../../data/20200419_lr=0.05_SH6S_cgu_ex_hard/record/"
    # fp3 = "../../data/20200420_lr=0.05_1blk_40_cg_ex_hard/record/"
    # fp2 = "../../data/20200422_lr=0.05_1blk_40_cgu_ex_hard/record/"
    # fp2 = "../data/20200329_lr=0.1/record/"

    # fp1 = "../../data/20200427_lr=0.05_SH6S_cgu_ex_hard/record/"
    # fp1 = "../../data/20200418_lr=0.05_SH6S_cg_ex_hard/record/"
    # fp1 = "../../data/20200429_lr=0.05_SH6S_cg_ex_hard/record/"
    # fp1 = "../../data/20200427_lr=0.05_SH6S_cgu_ex_hard/record/"
    # fp2 = "../../data/20200427_lr=0.05_1blk_40_cgu_ex_hard/record/"
    # fp2 = "../../data/20200429_lr=0.05_SH6S_cg_ex_hard/record/"
    # fp2 = "../../data/20200429_lr=0.05_1blk_40_cg_ex_hard/record/"
    # fp1 = "../../data/20200502_SH6S_cgu/record/"
    # fp2 = "../../data/20200427_lr=0.05_1blk_40_cgu_ex_hard/record/"
    # fp3 = "../../data/20200506_1blk_40_cgu/record/"

    # fp1 = "/home/centos/Projects/CGMD/data/20200914_lr=0.05_1blk_50_600K_cg_ex_hard/record/"
    # fp2 = "/home/centos/Projects/CGMD/data/20200917_lr=0.05_1blk_50_600K_cgu_ex_hard/record/"
    # fp3 = "/home/centos/Projects/CGMD/data/20200917_lr=0.05_1blk_50_600K_cgu_ex_hard/record/"
    #
    # r1 = ibi.Record(fp1, restart=True)
    # param1 = r1.get_param(40, 1)
    #
    # r2 = ibi.Record(fp2, restart=True)
    # param2 = r2.get_param(0, 0)
    #
    # # r3 = ibi.Record(fp3, restart=True)
    # # param3 = r3.get_param(32, 1)
    # import pandas as pd
    # import numpy as np
    # from pyMD import functions as f
    # import matplotlib.pyplot as plt
    #
    # cgu_path = "/home/centos/Projects/CGMD/data/target/1blk_50_600K_cgu_ex_hard/aver/"
    # cg_path = "/home/centos/Projects/CGMD/data/target/1blk_50_600K_cg_ex_hard/aver/"
    #
    # new_param = dict()
    # for para in param2:
    #     new_param[para] = dict()
    #     for type_tuple in param2[para]:
    #         if para == "non_bond":
    #             if type_tuple in param1[para]:
    #                 new_param[para][type_tuple] = param1[para][type_tuple]
    #             elif type_tuple[0][0].islower() and type_tuple[1][0].islower():
    #                 new_param[para][type_tuple] = param2[para][type_tuple]
    #             else:
    #                 if type_tuple[0][0].islower():
    #                     cg_type = ("U", type_tuple[1])
    #                 else:
    #                     cg_type = (type_tuple[0], "U")
    #                 cg_data = pd.read_csv(cg_path + "non_bond_%s_target.csv" % "_".join(cg_type))
    #                 cgu_data = pd.read_csv(cgu_path + "non_bond_%s_target.csv" % "_".join(type_tuple))
    #                 cg_m = f.find_maximum(cg_data["gr"])[0]
    #                 cgu_m = f.find_maximum(cgu_data["gr"])[0]
    #                 if cgu_m < cg_m:
    #                     param = np.append(param1[para][cg_type][cg_m - cgu_m:], np.zeros(cg_m - cgu_m))
    #                     new_param[para][type_tuple] = param * 0.8 + param2[para][type_tuple] * 0.2
    #                     plt.plot(cg_data.x, param2[para][type_tuple])
    #                     plt.plot(cg_data.x, new_param[para][type_tuple])
    #                     plt.title("_".join(type_tuple))
    #                     plt.ylim(-10, 10)
    #                     plt.show()
    #                 else:
    #                     new_param[para][type_tuple] = param1[para][cg_type]
    #         else:
    #             if type_tuple in param1[para]:
    #                 new_param[para][type_tuple] = param1[para][type_tuple]
    #             else:
    #                 new_param[para][type_tuple] = param2[para][type_tuple]
    #
    #         # if para == "non_bond" and type_tuple[0][0].islower() and type_tuple[1][0].islower():
    #         #     new_param[para][type_tuple] = param1[para][type_tuple]
    #         # else:
    #         #     new_param[para][type_tuple] = param2[para][type_tuple]
    #
    # #         # elif type_tuple in param3[para]:
    # #         #     new_param[para][type_tuple] = param3[para][type_tuple]
    # #         # else:
    # #         #     new_param[para][type_tuple] = param2[para][type_tuple]
    # #         #     print(para, type_tuple)
    # ibi.Record.dump_dict_json(new_param, fp3 + "param_iter_01_0.json")
    # # import numpy as np
    # # x = np.arange(1800) * 0.01 + 0.01
    # # for type_tuple in param2["non_bond"]:
    # #     if type_tuple[0] == type_tuple[1]:
    # #         print(type_tuple[0])
    # #         print(x[np.argmin(param2["non_bond"][type_tuple])])


    fp1 = "/home/centos/Projects/CGMD/data/20210401_lr=0.05_1blk_50_fg_ex_hard/record/"
    fp2 = "/home/centos/Projects/CGMD/data/20200410_lr=0.05_SH6S_fgh/record/"

    r1 = ibi.Record(fp1)
    param1 = r1.get_param(0, 0)

    r2 = ibi.Record(fp2)
    param2 = r2.get_param(135, 1)

    new_param = dict()
    for para in param1:
        new_param[para] = dict()
        for type_tuple in param1[para]:
            new_tt = list()
            for _type in type_tuple:
                if _type == "TO(2)":
                    new_tt.append("TO(1)")
                elif _type == "TO(1)":
                    new_tt.append("TO(0)")
                elif _type == "c1":
                    new_tt.append("c3")
                elif _type == "c3":
                    new_tt.append("c1")
                else:
                    new_tt.append(_type)
            new_tt = tuple(new_tt)
            if new_tt in param2[para]:
                new_param[para][type_tuple] = param2[para][new_tt]
            else:
                new_param[para][type_tuple] = param1[para][type_tuple]
    ibi.Record.dump_dict_json(new_param, fp1+"param_iter_01_1.json")