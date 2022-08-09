if __name__ == '__main__':
    type_mass_dict = {"TO(1)": 73.1142, "TO(2)": 72.1062, "TO(3)": 73.1142, "Es": 72.0593, "Me": 14.0267,
                      "Ph": 76.093640, "U": 58.039890}

    aa = {'TO(1)': 0.010474380679788413, 'TO(2)': 0.009723223696421531, 'Es': 0.006734062852657406, 'Ph': 0.00747707931138862, 'U': 0.004479780979015864, 'Me': 0.0007546265333446588}
    cg = {'TO(1)': 0.0007673239615788553, 'TO(2)': 0.000765766553625944, 'Es': 0.0007718419424313741, 'Ph': 0.0007545026051536193, 'U': 0.0007626133015589796, 'Me': 0.0007546265333446589}

    for _t in aa:
        new_mass = type_mass_dict[_t] * aa[_t] / cg[_t]
        print(_t, new_mass)