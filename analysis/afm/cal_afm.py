from pyMD import analysis as ans
import matplotlib.pyplot as plt
import numpy as np


def plot_2d(fig_name, _axis, _v, _save_path=None):
    # plt.figure(figsize=(6, 5))
    plt.figure(figsize=(3, 2))
    plt.axis("equal")
    # plt.imshow(_v, interpolation='bilinear', v max=1.6, vmin=0.6)  # cmap="gray",
    # plt.imshow(_v, interpolation='bilinear', vmax=1000, vmin=0)  # cmap="gray",
    # plt.imshow(_v, interpolation='bilinear', vmax=_v.max(), vmin=0)  # cmap="gray",
    # plt.imshow(_v, interpolation='bilinear', vmax=3E8, vmin=0)  # cmap="gray",
    plt.pcolormesh(_v, vmax=_v.max(), vmin=0)  # cmap="gray",
    plt.colorbar()
    plt.title(fig_name)
    plt.xlabel(_axis[1])
    plt.ylabel(_axis[0])
    if _save_path:
        plt.savefig("%s/%s.png" % (_save_path, fig_name), dpi=300)
    else:
        plt.show()
    plt.close()


def gradient_density(atom3d, _grid, _begin, _end, p_axis):
    shape = atom3d.lattice_parameter / grid
    shape = np.delete(shape, p_axis).astype(int)+1
    density = np.zeros(shape)
    d = _end - _begin
    for atom in atom3d.Atoms.values():
        at = atom.type
        if at == "TO(1)" or at == "TO(2)":
            continue
        co = atom.coordinate - atom3d.box_l
        if (co[p_axis] - _begin) * (co[p_axis] - _end) < 0:
            k = (_end - co[p_axis]) / d
            co = np.delete(co, p_axis)
            p = (co / grid).astype(int)
            density[p[0], p[1]] += atom3d.type_mass_dict[at] * k
    return density[:-1, :-1]


if __name__ == '__main__':
    fp = "/home/centos/work/SH6S/cgu_7blk_300chain/anneal_400/"
    start = 17200000

    fp_data = fp + "strain.data"
    fp_trj = fp + "record.lammpstrj."
    atom3ds = ans.atom3d_generator(fp_data, fp_trj, start=start, step=100000, frames=1)
    grid = 5
    begin = 80
    end = 30 + begin
    projection_axis = 1
    data = gradient_density(next(atom3ds), grid, begin, end, projection_axis)
    plot_2d("AFM", "xy", data)
