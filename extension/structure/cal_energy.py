from pyMD import functions as f
from pyMD import cartesian_operation as cso
import pickle as pkl
import copy as cp
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':

    with open("../../template/U1.pkl", "rb") as file:
        Uh = pkl.load(file)

    Uh2 = cp.deepcopy(Uh)

    Uh.set_centroid([0, 0, 0])
    Uh2.set_centroid([0, 0, 0])

    ds = np.arange(4., 5.5, 0.01)
    engs = list()
    for d in ds:
        Uh2.set_centroid([0, -d, 0.])
        eng = 0.
        for atom2 in Uh2.Atoms.values():
            for atom1 in Uh.Atoms.values():
                dist = cso.calculate_distance(atom1.coordinate, atom2.coordinate)
                eng += f.cal_aa_force(dist, (atom1.type, atom2.type), cut_off=14., q=True, dsf=True, a=0.01,
                                      dielectric=1)
        engs.append(eng)
    engs = np.asarray(engs)
    plt.plot(ds, engs)
    plt.xlabel("Distance ($\AA$)", fontdict={'family': 'Times New Roman', 'size': 14})
    plt.show()
    print(min(engs))
    print(ds[np.argmin(engs)])

# -11.86491883250703
# 4.569999999999977

# -15.978048287045823
# 4.489999999999979
