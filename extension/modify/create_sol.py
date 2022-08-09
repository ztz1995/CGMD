import os

if __name__ == '__main__':

    # hard_type = {"Ph":3, "Me":2, "c1":6, "o1":9, "h1":7, "n1":8, "TO(1)":4, "TO(2)":5, "Es":1}
    hard_type = {"Ph":3, "Me":2, "c1":6, "o1":9, "h1":7, "n1":8}
    # soft_type = {}

    fp = "/home/centos/ztz/change_h_bond/8blk_50_sv/"
    # files = os.listdir(fp)
    # for fn in files:
    #     if fn[:3] == "non":
    #         args = fn.split(".")
    #         args = args[0].split("_")
    #         t1 = args[2]
    #         t2 = args[3]
    #

    t1 = "sv"
    # t1 = "sv1"
    for t2 in hard_type:
        if t2 == "TO(1)":
            tt2 = "TO\(1\)"
        elif t2 == "TO(2)":
            tt2 = "TO\(2\)"
        else:
            tt2 = t2
        # origin_tuple = tuple(sorted(["TO\(2\)", tt2]))
        origin_tuple = tuple(sorted(["Ph", tt2]))
        # origin_tuple = tuple(sorted(["o1", tt2]))
        new_tuple = tuple(sorted([t1, tt2]))
        old_name = "non_bond_" + "_".join(origin_tuple) + ".param"
        new_name = "non_bond_" + "_".join(new_tuple) + ".param"

        os.system("cp %s%s %s%s" % (fp, old_name, fp, new_name))
        print("pair_coeff      %d  10 table 1 non_bond_%s_%s.param         %s_%s" % (hard_type[t2], t2, t1, t2, t1))
        new_name = "non_bond_" + "_".join([t2, t1]) + ".param"
        fn = "%s%s" % (fp, new_name)
        with open(fn, "r") as file:
            lines = file.readlines()
        lines[0] = "%s_%s\n" % (t2, t1)
        with open(fn, "w") as file:
            file.writelines(lines)