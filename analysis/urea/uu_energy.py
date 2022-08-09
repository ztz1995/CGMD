from pyMD import functions as f
from pyMD import cartesian_operation as cso
import pickle as pkl
import copy as cp
import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl
from pyMD import ibi


if __name__ == '__main__':

    with open("../../template/U1.pkl", "rb") as file:
        Uh = pkl.load(file)

    Uh2 = cp.deepcopy(Uh)

    Uh.set_centroid([0,0,0])
    Uh2.set_centroid([0,0,0])

    ds = np.arange(1., 14., 0.01)
    engs = list()
    for d in ds:
        Uh2.set_centroid([0, -d, 0.])
        # Uh2.set_centroid([0, 0., -d])
        eng = 0.
        for atom2 in Uh2.Atoms.values():
            for atom1 in Uh.Atoms.values():
                dist = cso.calculate_distance(atom1.coordinate, atom2.coordinate)
                eng += f.cal_aa_force(dist, (atom1.type, atom2.type), cut_off=14., q=True, dsf=True)
        engs.append(eng)
    engs = np.asarray(engs)
    plt.plot(ds, engs)
    plt.ylim(-15, 20)
    plt.axhline(y=-11.667, ls="--", c="red")
    plt.xlabel("Distance Angstrom", fontsize=14)
    # plt.ylabel("Energy kcal/mol", fontsize=14)
    plt.title("Urea-Urea interaction", fontsize=16)
    plt.show()
    print(min(engs))
    print(ds[np.argmin(engs)])

    dy =- f.cal_derivative(ds, engs)

    plt.plot(ds[200:], dy[200:])
    plt.ylim(-15, 20)
    print(ds[np.argmin(dy[200:])+200])
    print(min(dy[200:]))

    plt.show()

    # with open("data/uu_cgu_head-tail.pkl", "wb") as file:
    # with open("data/uu_cgu_parallel.pkl", "wb") as file:
    #     pkl.dump((ds,engs), file)
    #     pass
    # -11.86491883250703
    # 4.569999999999977

    # -15.978048287045823
    # 4.489999999999979

    # fp1 = "../../data/20200720_lr=0.05_1blk_50_cg_ex_hard/record/"
    # r1 = ibi.Record(fp1, restart=True)
    # # param1 = r1.get_param(22, 1)
    # param1 = r1.get_param(52, 1)
    #
    # x = ibi.IBI.cal_points["non_bond"]
    # y = param1["non_bond"][("U", "U")]
    # # plt.plot(x, y)
    # # plt.ylim(-15, 20)
    # # plt.show()
    # with open("data/uu_cg.pkl", "wb") as file:
    #     pkl.dump((x, y), file)