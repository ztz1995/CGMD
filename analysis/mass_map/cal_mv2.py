from pyMD import analysis as ans

if __name__ == '__main__':

    fp = "/home/centos/model/aa/1blk_50/0/"
    fp_data = fp + "1blk_50.data"
    fp_trj = fp + "1blk_50.lammpstrj"
    fp_info = fp + "cg_struct_info.pkl"

    bead_type = ["TO(1)", "TO(2)", "Es", "Ph", "U", "Me"]
    mv2_aa_dict = dict()
    mv2_cg_dict = dict()
    num_dict = dict()
    for _t in bead_type:
        mv2_aa_dict[_t] = 0.
        mv2_cg_dict[_t] = 0.
        num_dict[_t] = 0.

    atom3ds = ans.atom3d_generator_from_aa(fp_info, fp_data, fp_trj, start=15000000, step=1000, frames=100,
                                           one_file=True, fake=False, renew_v=True)
    for atom3d in atom3ds:
        for atom in atom3d.Atoms.values():
            _t = atom.type
            mv2_aa_dict[_t] += atom.mv2_aa
            mv2_cg_dict[_t] += atom.mv2_cg
            num_dict[_t] += 1
    for _t in bead_type:
        mv2_aa_dict[_t] /= num_dict[_t]
        mv2_cg_dict[_t] /= num_dict[_t]

    print(mv2_aa_dict)
    print(mv2_cg_dict)
