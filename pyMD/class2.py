from collections import OrderedDict

class COMPASS:
    atom_types = {'c3a', 'o2s', 'h1', "c3'", 'o2e', 'h1n', 'c4', 'c3"', 'o1=', 'n3mh', 'c4o'}
    bond_increment = {('c3a', 'c3a'): 0.0, ('c3a', 'c4'): 0.0, ('c3a', 'h1'): -0.1268, ('c4', 'c4'): 0.0,
                    ('c4', 'h1'): -0.053,
                    ('o2z', 'si4'): -0.2225, ('c3a', 'o2e'): 0.042, ('c3a', 'o2h'): 0.042, ('c4', 'o2e'): 0.16,
                    ('c4', 'o2h'): 0.16,
                    ('h1', 'o2'): 0.42, ('h1', 'o2h'): 0.41, ('c3a', 'n2='): 0.199, ('c3a', 'p4='): -0.06,
                    ('c4', 'n2='): 0.345,
                    ('c4', 'p4='): -0.05, ('cl1p', 'p4='): -0.12, ('f1p', 'p4='): -0.18, ('h1', 'n2='): 0.328,
                    ('h1', 'p4='): -0.05,
                    ('n2=', 'n2='): 0.0, ('n2=', 'n3'): 0.025, ('n2=', 'o2'): -0.043, ('n2=', 'p4='): -0.35,
                    ('n3', 'p4='): -0.12,
                    ('o2', 'p4='): -0.14, ('c1o', 'o1c'): -0.0203, ('c2=', 'o1='): 0.4, ('c2=', 's1='): 0.0258,
                    ('n2o', 'o1='): 0.073,
                    ('h1h', 'h1h'): 0.0, ('n1n', 'n1n'): 0.0, ('n1o', 'o1n'): 0.0288, ('o1=', 's2='): -0.2351,
                    ('o1o', 'o1o'): 0.0,
                    ('c3a', 'n3o'): 0.239, ('c4', 'n3o'): 0.21, ('c4', 'o2n'): 0.317, ('h1', 'n3o'): 0.188,
                    ('n3o', 'o1='): 0.428,
                    ('n3o', 'o2n'): 0.001, ("c3'", 'o2e'): 0.112, ("c3'", 'c4'): 0.0, ("c3'", 'o1='): 0.45,
                    ("c3'", 'c3a'): 0.035,
                    ("c3'", 'n3m'): 0.0, ('c3a', 'n3m'): 0.095, ('n2z', 'c4'): -0.311, ('n2z', 'h1'): -0.335,
                    ('n2t', 'n1t'): 0.386,
                    ('n2t', 'n2z'): 0.247, ('c3a', 'si4'): -0.117, ('c4', 'si4'): -0.135, ('h1', 'si4'): -0.126,
                    ('o1=', 'c3"'): -0.5,
                    ('n3mh', 'c3"'): -0.016, ('n3mh', 'h1n'): -0.351, ('n3mh', 'c3a'): -0.095, ('c4', 'c4o'): 0.,
                    ('c4o', 'o2e'): 0.16,
                    ('c4o', 'h1'): -0.053, ('c4o', 'o2s'): 0.16, ('o2s', "c3'"): -0.112}
    bond_increment_no_h = {('c3a', 'c3a'): 0.0, ('c3a', 'c4'): 0.0, ('c3a', 'h1'): -0.1268, ('c4', 'c4'): 0.0,
                    ('c4', 'h1'): -0.053,
                    ('o2z', 'si4'): -0.2225, ('c3a', 'o2e'): 0.042, ('c3a', 'o2h'): 0.042, ('c4', 'o2e'): 0.16,
                    ('c4', 'o2h'): 0.16,
                    ('h1', 'o2'): 0.42, ('h1', 'o2h'): 0.41, ('c3a', 'n2='): 0.199, ('c3a', 'p4='): -0.06,
                    ('c4', 'n2='): 0.345,
                    ('c4', 'p4='): -0.05, ('cl1p', 'p4='): -0.12, ('f1p', 'p4='): -0.18, ('h1', 'n2='): 0.328,
                    ('h1', 'p4='): -0.05,
                    ('n2=', 'n2='): 0.0, ('n2=', 'n3'): 0.025, ('n2=', 'o2'): -0.043, ('n2=', 'p4='): -0.35,
                    ('n3', 'p4='): -0.12,
                    ('o2', 'p4='): -0.14, ('c1o', 'o1c'): -0.0203, ('c2=', 'o1='): 0.4, ('c2=', 's1='): 0.0258,
                    ('n2o', 'o1='): 0.073,
                    ('h1h', 'h1h'): 0.0, ('n1n', 'n1n'): 0.0, ('n1o', 'o1n'): 0.0288, ('o1=', 's2='): -0.2351,
                    ('o1o', 'o1o'): 0.0,
                    ('c3a', 'n3o'): 0.239, ('c4', 'n3o'): 0.21, ('c4', 'o2n'): 0.317, ('h1', 'n3o'): 0.188,
                    ('n3o', 'o1='): 0.428,
                    ('n3o', 'o2n'): 0.001, ("c3'", 'o2e'): 0.112, ("c3'", 'c4'): 0.0, ("c3'", 'o1='): 0.45,
                    ("c3'", 'c3a'): 0.035,
                    ("c3'", 'n3m'): 0.0, ('c3a', 'n3m'): 0.095, ('n2z', 'c4'): -0.311, ('n2z', 'h1'): -0.335,
                    ('n2t', 'n1t'): 0.386,
                    ('n2t', 'n2z'): 0.247, ('c3a', 'si4'): -0.117, ('c4', 'si4'): -0.135, ('h1', 'si4'): -0.126,
                    ('o1=', 'c3"'): -0.05,
                    ('n3mh', 'c3"'): -0.016, ('n3mh', 'h1n'): -0.0351, ('n3mh', 'c3a'): -0.095, ('c4', 'c4o'): 0.,
                    ('c4o', 'o2e'): 0.16,
                    ('c4o', 'h1'): -0.053, ('c4o', 'o2s'): 0.16, ('o2s', "c3'"): -0.112}
    type_map = {
        'Atoms': {'c4': 1, 'h1': 2, 'c4o': 3, 'o2e': 4, 'o2s': 5, "c3'": 6, 'c3a': 7,
                'n3mh': 8, 'h1n': 9, 'c3\"': 10, 'o1=': 11},
        'Bonds': {('c4', 'c4'): 1, ('c4', 'h1'): 2, ('c4', 'c4o'): 3, ('c4o', 'o2e'): 4, ('c4o', 'h1'): 5, ('c4o', 'o2s'): 6, 
                ("c3'", 'o2s'): 7, ("c3'", 'o1='): 8, ("c3'", 'c3a'): 9, ('c3a', 'c3a'): 10, ('c3a', 'h1'): 11, 
                ('c3a', 'n3mh'): 12, ('h1n', 'n3mh'): 13, ('c3"', 'n3mh'): 14, ('c3"', 'o1='): 15, ('c3a', 'c4'): 16}, 
        'Angles': {('c4', 'c4', 'h1'): 1, ('h1', 'c4', 'h1'): 2, ('c4', 'c4', 'c4'): 3, ('c4', 'c4', 'c4o'): 4, 
                ('c4o', 'c4', 'h1'): 5, ('c4', 'c4o', 'o2e'): 6, ('c4', 'c4o', 'h1'): 7, ('h1', 'c4o', 'o2e'): 8,
                ('h1', 'c4o', 'h1'): 9, ('c4o', 'o2e', 'c4o'): 10, ('c4', 'c4o', 'o2s'): 11, ('h1', 'c4o', 'o2s'): 12, 
                ("c3'", 'o2s', 'c4o'): 13, ('c3a', "c3'", 'o1='): 14, ('o1=', "c3'", 'o2s'): 15, ('c3a', "c3'", 'o2s'): 16, 
                ('c3a', 'c3a', 'c3a'): 17, ("c3'", 'c3a', 'c3a'): 18, ('c3a', 'c3a', 'h1'): 19, ('c3a', 'c3a', 'n3mh'): 20, 
                ('c3"', 'n3mh', 'h1n'): 21, ('c3a', 'n3mh', 'h1n'): 22, ('c3"', 'n3mh', 'c3a'): 23, ('n3mh', 'c3"', 'o1='): 24,
                ('n3mh', 'c3"', 'n3mh'): 25, ('c3a', 'c3a', 'c4'): 26, ('c3a', 'c4', 'h1'): 27, ('c3a', 'c4', 'c3a'): 28}, 
        'Dihedrals': {('c4', 'c4', 'c4', 'h1'): 1, ('h1', 'c4', 'c4', 'h1'): 2, ('c4', 'c4', 'c4', 'c4'): 3,
                ('c4', 'c4', 'c4', 'c4o'): 4, ('c4o', 'c4', 'c4', 'h1'): 5, ('c4', 'c4', 'c4o', 'o2e'): 6,
                ('c4', 'c4', 'c4o', 'h1'): 7, ('h1', 'c4', 'c4o', 'o2e'): 8, ('h1', 'c4', 'c4o', 'h1'): 9, 
                ('c4', 'c4o', 'o2e', 'c4o'): 10, ('c4o', 'o2e', 'c4o', 'h1'): 11, ('c4o', 'c4', 'c4', 'c4o'): 12, 
                ('c4', 'c4', 'c4o', 'o2s'): 13, ('h1', 'c4', 'c4o', 'o2s'): 14, ("c3'", 'o2s', 'c4o', 'c4'): 15, 
                ("c3'", 'o2s', 'c4o', 'h1'): 16, ('c4o', 'o2s', "c3'", 'o1='): 17, ('c3a', "c3'", 'o2s', 'c4o'): 18, 
                ('c3a', 'c3a', "c3'", 'o1='): 19, ('c3a', 'c3a', "c3'", 'o2s'): 20, ('c3a', 'c3a', 'c3a', 'h1'): 21, 
                ('c3a', 'c3a', 'c3a', 'c3a'): 22, ("c3'", 'c3a', 'c3a', 'h1'): 23, ("c3'", 'c3a', 'c3a', 'c3a'): 24, 
                ('h1', 'c3a', 'c3a', 'h1'): 25, ('h1', 'c3a', 'c3a', 'n3mh'): 26, ('c3a', 'c3a', 'c3a', 'n3mh'): 27, 
                ('c3a', 'c3a', 'n3mh', 'h1n'): 28, ('c3"', 'n3mh', 'c3a', 'c3a'): 29, ('h1n', 'n3mh', 'c3"', 'n3mh'): 30, 
                ('h1n', 'n3mh', 'c3"', 'o1='): 31, ('c3a', 'n3mh', 'c3"', 'n3mh'): 32, ('c3a', 'n3mh', 'c3"', 'o1='): 33, 
                ('c4', 'c3a', 'c3a', 'h1'): 34, ('c3a', 'c3a', 'c3a', 'c4'): 35, ('c3a', 'c3a', 'c4', 'h1'): 36, 
                ('c3a', 'c3a', 'c4', 'c3a'): 37}, 
        'Impropers': {('c4', 'c4', 'h1', 'h1'): 1, ('h1', 'c4', 'h1', 'h1'): 2, ('c4', 'c4', 'c4', 'h1'): 3, 
                ('c4', 'c4', 'c4o', 'h1'): 4, ('c4o', 'c4', 'h1', 'h1'): 5, ('c4', 'c4o', 'h1', 'o2e'): 6, 
                ('c4', 'c4o', 'h1', 'h1'): 7, ('h1', 'c4o', 'h1', 'o2e'): 8, ('c4', 'c4o', 'h1', 'o2s'): 9, 
                ('h1', 'c4o', 'h1', 'o2s'): 10, ('c3a', "c3'", 'o1=', 'o2s'): 11, ("c3'", 'c3a', 'c3a', 'c3a'): 12, 
                ('c3a', 'c3a', 'c3a', 'h1'): 13, ('c3a', 'c3a', 'c3a', 'n3mh'): 14, ('c3"', 'n3mh', 'c3a', 'h1n'): 15, 
                ('n3mh', 'c3"', 'n3mh', 'o1='): 16, ('c3a', 'c3a', 'c3a', 'c4'): 17, ('c3a', 'c4', 'h1', 'h1'): 18, 
                ('c3a', 'c4', 'c3a', 'h1'): 19, ('h1', 'c4o', 'h1', 'h1'): 20}}
    param_types = OrderedDict([
                # ("Masses", "Atoms"),
                ("Pair Coeffs", "Atoms"),
                ("Bond Coeffs", "Bonds"),
                ("Angle Coeffs", "Angles"),
                ("Dihedral Coeffs", "Dihedrals"),
                ("Improper Coeffs", "Impropers"),
                ("BondBond Coeffs", "Angles"),
                ("BondAngle Coeffs", "Angles"),
                ("AngleAngle Coeffs", "Impropers"),
                ("AngleAngleTorsion Coeffs", "Dihedrals"),
                ("EndBondTorsion Coeffs", "Dihedrals"),
                ("MiddleBondTorsion Coeffs", "Dihedrals"),
                ("BondBond13 Coeffs", "Dihedrals"),
                ("AngleTorsion Coeffs", "Dihedrals")])
    params = {
        'Masses': {1: '12.011150', 2: '1.007970', 3: '12.011150', 4: '15.999400', 5: '15.999400', 6: '12.011150', 7: '12.011150', 8: '14.006700', 9: '1.007970', 10: '12.011150', 11: '15.999400'}, 
        'Pair Coeffs': {1: '0.0620000000    3.8540000000', 2: '0.0230000000    2.8780000000', 3: '0.0748000000    3.8700000000', 4: '0.1200000000    3.3000000000', 5: '0.0960000000    3.3000000000', 6: '0.0640000000    3.9000000000', 7: '0.0680000000    3.9150000000', 8: '0.2000000000    3.7000000000', 9: '0.0100000000    1.4500000000', 10: '0.0640000000    3.9000000000', 11: '0.1920000000    3.4300000000'}, 
        'Bond Coeffs': {1: '1.5300    299.6700    -501.7675    679.8014', 2: '1.1010    345.0000    -691.8975    844.5945', 3: '1.5300    299.6700    -501.7675    679.8014', 4: '1.4200    400.3954    -835.1848    1313.0166', 5: '1.1010    345.0000    -691.8975    844.5945', 6: '1.4200    400.3954    -835.1848    1313.0166', 7: '1.3750    368.7309    -832.4837    1274.0390', 8: '1.2160    823.7948    -1878.8288    2303.4950', 9: '1.4890    339.3574    -655.7403    670.2309', 10: '1.4170    470.8361    -627.6245    1327.6166', 11: '1.0982    372.8251    -803.4381    894.3328', 12: '1.3950    344.0452    -652.1377    1022.2271', 13: '1.0100    462.7500    -1053.6355    1545.7701', 14: '1.3880    440.6783    -828.3871    1423.2587', 15: '1.2160    823.7948    -1878.8288    2303.4950', 16: '1.5010    321.9021    -521.8355    572.1488'}, 
        'Angle Coeffs': {1: '110.7700    41.4530    -10.6037    5.1278', 2: '107.6600    39.6410    -12.9230    -2.4300', 3: '112.6700    39.5160    -7.4449    -9.5589', 4: '112.6700    39.5160    -7.4449    -9.5589', 5: '110.7700    41.4530    -10.6037    5.1278', 6: '111.2700    54.5381    -8.3662    -13.0836', 7: '110.7700    41.4530    -10.6037    5.1278', 8: '108.7280    58.5446    -10.8074    -12.3997', 9: '107.6600    39.6410    -12.9230    -2.4300', 10: '104.5000    35.7454    -10.0051    -6.2733', 11: '111.2700    54.5381    -8.3662    -13.0836', 12: '108.7280    58.5446    -10.8074    -12.3997', 13: '109.0000    38.9739    -6.2592    -8.1729', 14: '125.5320    72.3167    -16.0615    2.0827', 15: '118.9855    98.6813    -22.2526    10.3714', 16: '108.4400    84.8377    -19.9623    2.7402', 17: '118.9000    61.0226    -34.9904    0.0000', 18: '116.0640    71.2598    -15.8268    2.0523', 19: '117.9400    35.1558    -12.4698    -0.0000', 20: '120.7640    73.2738    -27.4044    13.3945', 21: '122.9480    40.4820    -16.2009    8.3271', 22: '116.3230    18.3123    -7.8322    5.3289', 23: '120.0700    47.1131    -32.5599    13.1257', 24: '125.5320    101.8765    -41.8101    7.7222', 25: '120.5292    100.0857    -36.7315    24.2608', 26: '120.0500    44.7148    -22.7330    0.0000', 27: '111.0000    44.3234    -9.4453    0.0000', 28: '111.0000    44.3234    -9.4453    0.0000'}, 
        'Dihedral Coeffs': {1: '0.0000    0.0000    0.0316    0.0000    -0.1681    0.0000', 2: '-0.1432    0.0000    0.0617    0.0000    -0.1531    0.0000', 3: '0.0001    0.0000    0.0514    0.0000    -0.1430    0.0000', 4: '0.0001    0.0000    0.0514    0.0000    -0.1430    0.0000', 5: '0.0000    0.0000    0.0316    0.0000    -0.1681    0.0000', 6: '0.7138    0.0000    0.2660    0.0000    -0.2545    0.0000', 7: '0.0000    0.0000    0.0316    0.0000    -0.1681    0.0000', 8: '-0.1432    0.0000    0.2531    0.0000    -0.0908    0.0000', 9: '-0.1432    0.0000    0.0617    0.0000    -0.1531    0.0000', 10: '-0.4000    0.0000    -0.4028    0.0000    -0.2450    0.0000', 11: '0.5299    0.0000    0.0002    0.0000    -0.3960    0.0000', 12: '0.0001    0.0000    0.0514    0.0000    -0.1430    0.0000', 13: '0.7138    0.0000    0.2660    0.0000    -0.2545    0.0000', 14: '-0.1432    0.0000    0.2531    0.0000    -0.0908    0.0000', 15: '0.1303    0.0000    -0.3250    0.0000    0.1134    0.0000', 16: '0.9511    0.0000    0.1158    0.0000    0.0003    0.0000', 17: '0.8907    0.0000    3.2639    0.0000    0.2647    0.0000', 18: '-2.5595    0.0000    2.2013    0.0000    0.0327    0.0000', 19: '-0.0001    0.0000    0.7799    0.0000    -0.0000    0.0000', 20: '-0.0000    0.0000    0.7800    0.0000    -0.0000    0.0000', 21: '0.0000    0.0000    3.9661    0.0000    0.0000    0.0000', 22: '8.3667    0.0000    1.2000    0.0000    0.0000    0.0000', 23: '-0.0002    0.0000    2.1672    0.0000    0.0002    0.0000', 24: '0.0000    0.0000    4.6282    0.0000    0.0000    0.0000', 25: '0.0005    0.0000    2.3497    0.0000    -0.0005    0.0000', 26: '0.0003    0.0000    3.4041    0.0000    -0.0002    0.0000', 27: '0.0001    0.0000    3.4040    0.0000    0.0000    0.0000', 28: '0.0000    0.0000    0.5107    0.0000    0.0000    0.0000', 29: '0.0001    0.0000    0.5107    0.0000    0.0000    0.0000', 30: '-1.0630    0.0000    1.5631    0.0000    0.0000    0.0000', 31: '0.0005    0.0000    2.0516    0.0000    -0.0013    0.0000', 32: '-1.0631    0.0000    1.5631    0.0000    -0.0002    0.0000', 33: '-0.0003    0.0000    2.0522    0.0000    0.0003    0.0000', 34: '0.0000    0.0000    1.5590    0.0000    0.0000    0.0000', 35: '-0.0000    0.0000    4.4072    0.0000    -0.0000    0.0000', 36: '-0.2802    0.0000    -0.0678    0.0000    -0.0122    0.0000', 37: '-0.2802    0.0000    -0.0678    0.0000    -0.0122    0.0000'}, 
        'Improper Coeffs': {1: '-0.0000    0.0000', 2: '-0.0000    0.0000', 3: '-0.0000    0.0000', 4: '-0.0000    0.0000', 5: '-0.0000    0.0000', 6: '-0.0000    0.0000', 7: '-0.0000    0.0000', 8: '-0.0000    0.0000', 9: '-0.0000    0.0000', 10: '-0.0000    0.0000', 11: '49.3740    0.0000', 12: '17.0526    0.0000', 13: '4.8912    0.0000', 14: '17.0526    0.0000', 15: '4.4181    0.0000', 16: '59.3740    0.0000', 17: '7.8153    0.0000', 18: '-0.0000    0.0000', 19: '0.0000    0.0000', 20: '-0.0000    0.0000'}, 
        'BondBond Coeffs': {1: '3.3872    1.5300    1.1010', 2: '5.3316    1.1010    1.1010', 3: '0.0000    1.5300    1.5300', 4: '0.0000    1.5300    1.5300', 5: '3.3872    1.5300    1.1010', 6: '11.4318    1.5300    1.4200', 7: '3.3872    1.5300    1.1010', 8: '23.1979    1.1010    1.4200', 9: '5.3316    1.1010    1.1010', 10: '-7.1131    1.4200    1.4200', 11: '11.4318    1.5300    1.4200', 12: '23.1979    1.1010    1.4200', 13: '26.1360    1.3750    1.4200', 14: '116.9445    1.4890    1.2160', 15: '210.1813    1.2160    1.3750', 16: '69.9445    1.4890    1.3750', 17: '68.2856    1.4170    1.4170', 18: '37.8749    1.4890    1.4170', 19: '1.0795    1.4170    1.0982', 20: '37.8749    1.4170    1.3950', 21: '8.6253    1.3880    1.0100', 22: '8.2930    1.3950    1.0100', 23: '41.4233    1.3880    1.3950', 24: '115.4645    1.3880    1.2160', 25: '84.5263    1.3880    1.3880', 26: '12.0676    1.4170    1.5010', 27: '2.9168    1.5010    1.1010', 28: '-0.0000    1.5010    1.5010'}, 
        'BondAngle Coeffs': {1: '20.7540    11.4210    1.5300    1.1010', 2: '18.1030    18.1030    1.1010    1.1010', 3: '8.0160    8.0160    1.5300    1.5300', 4: '8.0160    8.0160    1.5300    1.5300', 5: '20.7540    11.4210    1.5300    1.1010', 6: '2.6868    20.4033    1.5300    1.4200', 7: '20.7540    11.4210    1.5300    1.1010', 8: '4.6189    55.3270    1.1010    1.4200', 9: '18.1030    18.1030    1.1010    1.1010', 10: '-2.8112    -2.8112    1.4200    1.4200', 11: '2.6868    20.4033    1.5300    1.4200', 12: '4.6189    55.3270    1.1010    1.4200', 13: '21.5366    -16.6748    1.3750    1.4200', 14: '72.8758    76.1093    1.4890    1.2160', 15: '79.4497    57.0987    1.2160    1.3750', 16: '72.8758    76.1093    1.4890    1.3750', 17: '28.8708    28.8708    1.4170    1.4170', 18: '23.6977    45.8865    1.4890    1.4170', 19: '20.0033    24.2183    1.4170    1.0982', 20: '35.8865    53.6977    1.4170    1.3950', 21: '34.8312    15.0778    1.3880    1.0100', 22: '10.4568    12.8217    1.3950    1.0100', 23: '34.7791    24.3705    1.3880    1.3950', 24: '32.8758    46.1093    1.3880    1.2160', 25: '49.0875    49.0875    1.3880    1.3880', 26: '31.0771    47.0579    1.4170    1.5010', 27: '26.4608    11.7717    1.5010    1.1010', 28: '0.0000    0.0000    1.5010    1.5010'}, 
        'AngleAngle Coeffs': {1: '0.2738    -0.4825    0.2738    110.7700    110.7700    107.6600', 2: '-0.3157    -0.3157    -0.3157    107.6600    107.6600    107.6600', 3: '-1.3199    -1.3199    0.1184    112.6700    110.7700    110.7700', 4: '-1.3199    -1.3199    0.1184    112.6700    110.7700    110.7700', 5: '0.2738    -0.4825    0.2738    110.7700    110.7700    107.6600', 6: '3.9177    2.5926    0.1689    110.7700    111.2700    108.7280', 7: '0.2738    -0.4825    0.2738    110.7700    110.7700    107.6600', 8: '2.4259    2.4259    2.1283    107.6600    108.7280    108.7280', 9: '3.9177    2.5926    0.1689    110.7700    111.2700    108.7280', 10: '2.4259    2.4259    2.1283    107.6600    108.7280    108.7280', 11: '-0.0000    -0.0000    -0.0000    125.5320    108.4400    118.9855', 12: '0.0000    5.9863    0.0000    116.0640    116.0640    118.9000', 13: '0.0000    0.0000    0.0000    118.9000    117.9400    117.9400', 14: '0.0000    0.0000    0.0000    118.9000    120.7640    120.7640', 15: '0.0000    0.0000    0.0000    120.0700    122.9480    116.3230', 16: '-0.0001    -0.0001    -0.0001    120.5292    125.5320    125.5320', 17: '-0.0000    -0.0000    0.0000    118.9000    120.0500    120.0500', 18: '2.3794    3.0118    2.3794    111.0000    111.0000    107.6600', 19: '0.0000    0.0000    0.0000    111.0000    111.0000    111.0000', 20: '-0.3157    -0.3157    -0.3157    107.6600    107.6600    107.6600'}, 
        'AngleAngleTorsion Coeffs': {1: '-16.1640    112.6700    110.7700', 2: '-12.5639    110.7700    110.7700', 3: '-22.0450    112.6700    112.6700', 4: '-22.0450    112.6700    112.6700', 5: '-16.1640    112.6700    110.7700', 6: '-29.0420    112.6700    111.2700', 7: '-16.1640    112.6700    110.7700', 8: '-20.2006    110.7700    111.2700', 9: '-12.5639    110.7700    110.7700', 10: '-19.0059    111.2700    104.5000', 11: '-16.4441    104.5000    108.7280', 12: '-22.0450    112.6700    112.6700', 13: '-29.0420    112.6700    111.2700', 14: '-20.2006    110.7700    111.2700', 15: '-15.7082    109.0000    111.2700', 16: '-13.1502    109.0000    108.7280', 17: '-32.9362    109.0000    118.9855', 18: '-0.0004    108.4400    109.0000', 19: '-0.0000    116.0640    125.5320', 20: '0.0000    116.0640    108.4400', 21: '-4.8141    118.9000    117.9400', 22: '-0.0000    118.9000    118.9000', 23: '0.0016    116.0640    117.9400', 24: '-0.0000    116.0640    118.9000', 25: '0.3603    117.9400    117.9400', 26: '0.0006    117.9400    120.7640', 27: '-0.0000    118.9000    120.7640', 28: '-0.0000    120.7640    116.3230', 29: '0.0000    120.0700    120.7640', 30: '0.0006    122.9480    120.5292', 31: '-0.0003    122.9480    125.5320', 32: '-0.0000    120.0700    120.5292', 33: '0.0008    120.0700    125.5320', 34: '4.4444    120.0500    117.9400', 35: '-14.4097    118.9000    120.0500', 36: '-5.8888    120.0500    111.0000', 37: '0.0000    120.0500    111.0000'}, 
        'EndBondTorsion Coeffs': {1: '0.2486    0.2422    -0.0925    0.0814    0.0591    0.2219    1.5300    1.1010', 2: '0.2130    0.3107    0.0759    0.2130    0.3107    0.0759    1.1010    1.1010', 3: '-0.0732    -0.0000    -0.0000    -0.0732    -0.0000    -0.0000    1.5300    1.5300', 4: '-0.0732    -0.0000    -0.0000    -0.0732    -0.0000    -0.0000    1.5300    1.5300', 5: '0.2486    0.2422    -0.0925    0.0814    0.0591    0.2219    1.5300    1.1010', 6: '-0.3190    0.4411    -0.7174    1.1538    0.8409    -0.9138    1.5300    1.4200', 7: '0.2486    0.2422    -0.0925    0.0814    0.0591    0.2219    1.5300    1.1010', 8: '0.9681    0.9562    0.0444    0.5922    0.6708    0.8616    1.1010    1.4200', 9: '0.2130    0.3107    0.0759    0.2137    0.3126    0.0778    1.1010    1.1010', 10: '-0.2456    1.0517    -0.7795    0.4741    1.2635    0.5576    1.5300    1.4200', 11: '-0.1632    0.1545    -1.1430    -0.6062    1.3314    0.9629    1.4200    1.1010', 12: '-0.0732    -0.0000    -0.0000    -0.0732    -0.0000    -0.0000    1.5300    1.5300', 13: '-0.3190    0.4411    -0.7174    1.1538    0.8409    -0.9138    1.5300    1.4200', 14: '0.9681    0.9562    0.0444    0.5922    0.6708    0.8616    1.1010    1.4200', 15: '0.2560    0.8133    -0.0728    -1.2164    -0.1715    -0.0964    1.3750    1.5300', 16: '0.2271    2.2990    -0.4481    0.9593    0.9182    -0.6016    1.3750    1.1010', 17: '0.0881    -2.4294    -0.7416    -4.2416    10.1104    1.6831    1.4200    1.2160', 18: '-0.0004    -0.0009    -0.0007    -0.0010    -0.0015    -0.0022    1.4890    1.4200', 19: '0.0000    0.0000    0.0000    0.0000    -0.0000    0.0000    1.4170    1.2160', 20: '0.0000    0.0000    0.0000    0.0000    0.0000    -0.0000    1.4170    1.3750', 21: '0.0000    -6.8958    -0.0000    -0.0000    -0.4669    -0.0000    1.4170    1.0982', 22: '-0.1185    6.3204    0.0000    -0.1185    6.3204    0.0000    1.4170    1.4170', 23: '-0.0001    -0.0002    0.0001    0.0000    -0.0010    -0.0003    1.4890    1.0982', 24: '-0.0000    0.0000    -0.0000    0.0000    0.0000    -0.0000    1.4890    1.4170', 25: '-0.0002    -0.6895    -0.0005    -0.0002    -0.6895    -0.0005    1.0982    1.0982', 26: '-0.0003    -0.0012    -0.0015    0.0006    0.0004    0.0012    1.0982    1.3950', 27: '0.0000    -0.0000    -0.0000    -0.0000    0.0000    0.0000    1.4170    1.3950', 28: '0.0000    -0.0000    -0.0000    -0.0000    -0.0000    -0.0000    1.4170    1.0100', 29: '-0.0000    -0.0000    -0.0000    0.0000    0.0000    0.0000    1.3880    1.4170', 30: '-0.0006    -0.0015    -0.0019    -0.0002    0.0008    0.0004    1.0100    1.3880', 31: '-0.0006    -0.0016    -0.0004    0.0013    0.0005    0.0008    1.0100    1.2160', 32: '0.0003    -0.0005    0.0004    0.0001    -0.0001    0.0001    1.3950    1.3880', 33: '-0.0002    -0.0010    -0.0005    -0.0006    -0.0011    -0.0020    1.3950    1.2160', 34: '-0.0000    -1.7970    -0.0000    -0.0000    -0.4879    -0.0000    1.5010    1.0982', 35: '0.0000    -0.6918    0.0000    -0.0000    0.2421    -0.0000    1.4170    1.5010', 36: '-0.5835    1.1220    0.3978    1.3997    0.7756    0.0000    1.4170    1.1010', 37: '0.0000    0.0000    0.0000    -0.0000    -0.0000    -0.0000    1.4170    1.5010'}, 
        'MiddleBondTorsion Coeffs': {1: '-14.8790    -3.6581    -0.3138    1.5300', 2: '-14.2610    -0.5327    -0.4864    1.5300', 3: '-17.7870    -7.1877    -0.0000    1.5300', 4: '-17.7870    -7.1877    -0.0000    1.5300', 5: '-14.8790    -3.6581    -0.3138    1.5300', 6: '-21.8842    -7.6764    -0.6868    1.5300', 7: '-14.8790    -3.6581    -0.3138    1.5300', 8: '-16.7985    -1.2328    -0.2777    1.5300', 9: '-14.2610    -0.5327    -0.4864    1.5300', 10: '-5.9288    -2.7007    -0.3175    1.4200', 11: '-6.7996    -4.6544    -1.4093    1.4200', 12: '-17.7870    -7.1877    -0.0000    1.5300', 13: '-21.8842    -7.6764    -0.6868    1.5300', 14: '-16.7985    -1.2328    -0.2778    1.5300', 15: '9.9416    2.6421    2.2333    1.4200', 16: '7.7139    4.2542    -1.0128    1.4200', 17: '0.4554    7.3093    0.2840    1.3750', 18: '0.0001    -0.0012    -0.0013    1.3750', 19: '0.0000    2.4002    -0.0000    1.4890', 20: '-0.0000    2.4002    -0.0000    1.4890', 21: '0.0000    -1.1521    -0.0000    1.4170', 22: '27.5989    -2.3120    0.0000    1.4170', 23: '0.0006    0.0004    0.0009    1.4170', 24: '0.0000    3.8762    0.0000    1.4170', 25: '-0.0017    4.8198    -0.0019    1.4170', 26: '0.0006    5.2023    -0.0005    1.4170', 27: '0.0000    5.2012    0.0000    1.4170', 28: '-0.0000    2.4730    0.0000    1.3950', 29: '-0.0000    4.9027    -0.0000    1.3950', 30: '0.0001    6.3289    -0.0007    1.3880', 31: '0.0009    4.4689    -0.0008    1.3880', 32: '0.0009    6.3274    -0.0001    1.3880', 33: '-0.0024    4.4645    -0.0045    1.3880', 34: '0.0000    3.9421    -0.0000    1.4170', 35: '-0.0000    9.1792    -0.0000    1.4170', 36: '-5.5679    1.4083    0.3010    1.5010', 37: '-0.0000    -0.0000    0.0000    1.5010'}, 
        'BondBond13 Coeffs': {1: '-0.0000    1.5300    1.1010', 2: '-0.0037    1.1010    1.1010', 3: '-0.0000    1.5300    1.5300', 4: '-0.0000    1.5300    1.5300', 5: '-0.0000    1.5300    1.1010', 6: '0.0000    1.5300    1.4200', 7: '-0.0000    1.5300    1.1010', 8: '-0.0012    1.1010    1.4200', 9: '-0.0037    1.1010    1.1010', 10: '-0.0000    1.5300    1.4200', 11: '0.0006    1.4200    1.1010', 12: '-0.0000    1.5300    1.5300', 13: '-0.0000    1.5300    1.4200', 14: '-0.0012    1.1010    1.4200', 15: '-0.0000    1.3750    1.5300', 16: '-0.0008    1.3750    1.1010', 17: '0.0008    1.4200    1.2160', 18: '0.0011    1.4890    1.4200', 19: '0.0000    1.4170    1.2160', 20: '-0.0000    1.4170    1.3750', 21: '-6.2741    1.4170    1.0982', 22: '53.0000    1.4170    1.4170', 23: '-0.0002    1.4890    1.0982', 24: '0.0000    1.4890    1.4170', 25: '-1.7042    1.0982    1.0982', 26: '-0.0014    1.0982    1.3950', 27: '-0.0000    1.4170    1.3950', 28: '0.0000    1.4170    1.0100', 29: '-0.0000    1.3880    1.4170', 30: '0.0021    1.0100    1.3880', 31: '-0.0000    1.0100    1.2160', 32: '-0.0018    1.3950    1.3880', 33: '0.0009    1.3950    1.2160', 34: '0.8743    1.5010    1.0982', 35: '2.5085    1.4170    1.5010', 36: '-3.4826    1.4170    1.1010', 37: '0.0000    1.4170    1.5010'}, 
        'AngleTorsion Coeffs': {1: '-0.2454    -0.0000    -0.1136    0.3110    0.4515    -0.1988    112.6700    110.7700', 2: '-0.8084    0.5581    -0.2454    -0.8084    0.5581    -0.2454    110.7700    110.7700', 3: '0.3886    -0.3139    0.1389    0.3886    -0.3139    0.1389    112.6700    112.6700', 4: '0.3886    -0.3139    0.1389    0.3886    -0.3139    0.1389    112.6700    112.6700', 5: '-0.2454    -0.0000    -0.1136    0.3110    0.4515    -0.1988    112.6700    110.7700', 6: '0.5623    -0.3041    -0.4015    0.9667    -0.7569    -1.2332    112.6700    111.2700', 7: '-0.2454    -0.0000    -0.1136    0.3110    0.4515    -0.1988    112.6700    110.7700', 8: '2.3676    2.4930    -1.0112    -0.1897    0.4917    0.7277    110.7700    111.2700', 9: '-0.8084    0.5581    -0.2454    -0.8086    0.5568    -0.2465    110.7700    110.7700', 10: '-2.7466    1.4877    -0.8955    0.5676    0.9450    0.0703    111.2700    104.5000', 11: '-0.7782    0.4336    -0.6655    -1.8242    1.6385    0.5135    104.5000    108.7280', 12: '0.3886    -0.3139    0.1389    0.3886    -0.3139    0.1389    112.6700    112.6700', 13: '0.5623    -0.3041    -0.4015    0.9667    -0.7569    -1.2332    112.6700    111.2700', 14: '2.3676    2.4930    -1.0112    -0.1897    0.4917    0.7277    110.7700    111.2700', 15: '-0.0890    -0.9159    0.7229    -0.4625    1.4489    -0.6766    109.0000    111.2700', 16: '-0.3142    -0.8698    0.0974    -0.4994    2.8057    -0.0403    109.0000    108.7280', 17: '2.3584    1.0065    -0.0324    5.9737    2.7270    1.9053    109.0000    118.9855', 18: '0.0003    0.0005    0.0005    0.0002    0.0004    -0.0000    108.4400    109.0000', 19: '0.0000    0.0000    -0.0000    0.0006    0.0003    0.0001    116.0640    125.5320', 20: '0.0000    -0.0000    -0.0000    0.0002    0.0001    0.0000    116.0640    108.4400', 21: '-0.0000    2.5014    -0.0000    -0.0002    2.7146    -0.0000    118.9000    117.9400', 22: '1.9767    1.0239    -0.0000    1.9767    1.0239    -0.0000    118.9000    118.9000', 23: '0.0014    -0.0004    -0.0006    -0.0003    -0.0008    -0.0004    116.0640    117.9400', 24: '-0.0000    -0.0000    -0.0000    -0.0002    -0.0001    -0.0000    116.0640    118.9000', 25: '0.0002    2.4505    0.0002    0.0002    2.4505    0.0002    117.9400    117.9400', 26: '0.0008    0.0001    0.0003    -0.0006    -0.0003    0.0000    117.9400    120.7640', 27: '-0.0000    -0.0000    -0.0000    -0.0005    -0.0002    -0.0001    118.9000    120.7640', 28: '0.0000    0.0000    0.0000    -0.0003    -0.0001    -0.0000    120.7640    116.3230', 29: '-0.0000    -0.0000    -0.0000    -0.0005    -0.0002    -0.0001    120.0700    120.7640', 30: '0.0005    -0.0002    -0.0002    -0.0008    -0.0004    -0.0002    122.9480    120.5292', 31: '-0.0010    -0.0010    -0.0008    0.0015    0.0015    0.0015    122.9480    125.5320', 32: '0.0005    0.0005    0.0005    0.0000    0.0001    0.0002    120.0700    120.5292', 33: '0.0011    0.0008    0.0008    0.0004    0.0002    -0.0002    120.0700    125.5320', 34: '0.0000    -0.1242    -0.0000    -0.0002    3.4600    -0.0000    120.0500    117.9400', 35: '-0.0000    3.8987    -0.0000    0.0000    -4.4683    0.0000    118.9000    120.0500', 36: '0.2251    0.6548    0.1237    4.6266    0.1632    0.0461    120.0500    111.0000', 37: '0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    120.0500    111.0000'}}
    
    @classmethod
    def create_potential_lines(cls, atom_type_dict, joint_type_dict):
        lines = "\n"
        for p_name, j_type in cls.param_types.items():
            lines += "\n%s\n\n" % p_name
            if j_type == "Atoms":
                for k, v in atom_type_dict.items():
                    lines += "%s    %s\n" % (str(v).rjust(8), cls.params[p_name][cls.type_map[j_type][k]])
            else:
                for k, v in joint_type_dict[j_type[:-1]].items():
                    p = cls.params[p_name][cls.type_map[j_type][k]] if k in cls.type_map[j_type] else cls.params[p_name][cls.type_map[j_type][k[::-1]]]
                    lines += "%s    %s\n" % (str(v).rjust(8), p)
        return lines

