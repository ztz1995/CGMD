import numpy as np

if __name__ == '__main__':

    sigma = 15.
    dt = 983.91 / 5.92 * sigma
    f_unit = 0.59584 / sigma
    s0 = f_unit * np.sqrt(dt) * 3

    # type_dict = {1: "Es", 2: "Me", 3: "Ph", 4: "TO(1)", 5: "TO(2)", 6: "U"}
    # type_type_dict = {1: 0, 2: 1, 3: 1, 4: 0, 5: 0, 6: 2}

    type_dict = {1: "Es", 2: "Ph", 3: "TO(1)", 4: "TO(2)", 5: "U"}
    type_type_dict = {1: 0, 2: 1, 3: 0, 4: 0, 5: 2}

    type_mass_dict = {"TO(1)": 73.1142, "TO(2)": 72.1062, "TO(3)": 73.1142, "Es": 72.0593, "Me": 14.0267,
                      "Ph": 76.093640, "U": 58.039890}

    print("pair_style      dpd/fdt 300.0 %.3f 245455" % sigma)
    # print("pair_style      hybrid dpd/fdt 300.0 15.000 245455 dpd/fdt 300.0 15.000 245455 zero 15.000 nocoeff")
    for i in range(1, len(type_dict) + 1):
        for j in range(i, len(type_dict) + 1):
            cutoff = 15.0
            a = 0.
            tti = type_type_dict[i]
            ttj = type_type_dict[j]
            s = np.power(type_mass_dict[type_dict[i]] * type_mass_dict[type_dict[j]] / 68.9 / 68.9, 1 / 8) * s0

            if tti == ttj:
                # UU
                if tti == 2:
                    # a = 20.
                    a = 25.

                # other self
                else:
                    a = 25.
            else:
                if tti * ttj == 0:
                    # soft U
                    if tti == 2 or ttj == 2:
                        a = 91.15
                        # a = (a-25)/2 + 25.
                    # soft Ph
                    else:
                        a = 29.48
                        # a = (a-25)/2 +25.

                else:
                    # U Ph
                    a = 55.13
                    # a = (a - 25) / 2 + 25.

            if a == 0:
                print("pair_coeff      %d  %d  zero # %s %s" % (i, j, type_dict[i], type_dict[j]))
            else:

                print("pair_coeff      %d  %d  %f  %f  %.1f # %s %s" % (
                    i, j, a * f_unit, s, cutoff, type_dict[i], type_dict[j]))