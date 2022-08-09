import numpy as np

if __name__ == '__main__':

    sigma = 15.
    dt = 983.91 / 5.92 * sigma
    f_unit = 0.59584 / sigma
    s0 = f_unit * np.sqrt(dt) * 3

    # type_dict = {1: "Es", 2: "Me", 3: "Ph", 4: "TO(1)", 5: "TO(2)", 6: "U"}
    type_dict = {1: "Es", 2: "Me", 3: "Ph", 4: "TO(1)", 5: "TO(2)", 6: "c1", 7: "h1", 8: "n1", 9: "o1"}

    type_type_dict = {1: 0, 2: 1, 3: 1, 4: 0, 5: 0, 6: 2, 7: 2, 8: 2, 9: 2}

    type_mass_dict = {"TO(1)": 73.1142, "TO(2)": 72.1062, "TO(3)": 73.1142, "Es": 72.0593, "Me": 14.0267,
                      "Ph": 76.093640, "U": 58.039890,
                      "c1": 12.0107, "o1": 15.9994, "n1": 14.0067, "h1": 1.008, "sv1": 16, "sv2": 56, "sv": 75}

    # u_dict = {'h1_h1': 2.4783786187661985, 'h1_o1': -3.527486993250609, 'h1_n1': -3.2523203479424376, 'c1_h1': 3.766733774623837, 'o1_o1': 5.651929475155983, 'n1_o1': 5.701033701243296, 'c1_o1': -4.499307466178653, 'n1_n1': 5.937447895319012, 'c1_n1': -3.6759951904238393, 'c1_c1': 6.644801015629236}
    u_dict = {'h1_h1': 0.6226941764320029, 'h1_o1': -0.8783465820371735, 'h1_n1': -0.802844472770812, 'c1_h1': 0.9582944730989507, 'o1_o1': 2.1854777298477583, 'n1_o1': 2.5912701810604997, 'c1_o1': -0.25572449853446116, 'n1_n1': 3.180465563969187, 'c1_n1': 0.3202033424731127, 'c1_c1': 2.570402429401341}

    # print("pair_style      dpd/fdt 300.0 %.3f 245455" % sigma)
    print("pair_style      hybrid dpd/fdt 300.0 15.000 245455 dpd/fdt 300.0 15.000 245455 zero 15.000 nocoeff")
    for i in range(1, len(type_dict) + 1):
        for j in range(i, len(type_dict) + 1):
            cutoff = 15.0
            a = 0.
            tti = type_type_dict[i]
            ttj = type_type_dict[j]
            hybrid = 1
            # s = np.power(type_mass_dict[type_dict[i]] * type_mass_dict[type_dict[j]] / 68.9 / 68.9, 1 / 8) * s0
            s = np.power(type_mass_dict[type_dict[i]] * type_mass_dict[type_dict[j]] / 68.9 / 68.9, 1 / 4) * s0
            # s =  s0

            if tti == ttj:
                # UU
                if tti == 2:
                    # hybrid = 2
                    # # cutoff = 8.
                    # cutoff = 15.
                    # # if i == 7 and j == 9:
                    # #     a = 0.01
                    # # elif i == j and i in [7, 8, 9]:
                    # #     a = 1.
                    # # else:
                    # #     a = 0.3
                    # a = u_dict["_".join([type_dict[i], type_dict[j]])] * 4

                    hybrid = 2
                    cutoff = 15.
                    if i == j and i == 6:
                        a = 25.
                    else:
                        a = 0.0

                # other self
                else:
                    # Ph Ph
                    if i == j == 3:
                        # a = 25
                        a = 25.
                    else:
                        a = 25.
            else:
                if tti * ttj == 0:
                    # soft U
                    if tti == 2 or ttj == 2:
                        if i == 6 or j == 6:
                            # a = 91.15
                            a = 45.
                        else:
                            a = 0
                    # soft Ph
                    else:
                        # a = 29.48
                        a = 40.
                else:
                    # U Ph
                    if i == 6 or j == 6:
                        # a = 55.13
                        a = 30.
                    else:
                        a = 0
                    # a = 35.
            if a == 0:
                print("pair_coeff      %d  %d  zero # %s %s" % (i, j, type_dict[i], type_dict[j]))
            else:
                # print("pair_coeff      %d  %d  dpd/fdt  %d  %f  %f  %.1f # %s %s" % (
                #     i, j, hybrid, a * f_unit, s, cutoff, type_dict[i], type_dict[j]))
                print("pair_coeff      %d  %d  dpd/fdt  %f  %f  %.1f # %s %s" % (
                    i, j, a * f_unit, s, cutoff, type_dict[i], type_dict[j]))