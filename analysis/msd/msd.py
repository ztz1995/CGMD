from pyMD import analysis as ans


if __name__ == '__main__':
    fp = "/home/centos/model/aa/1blk_40/long/1/"
    fp_data = fp + "1blk_40.data"
    fp_trj = fp + "trj/1blk_40.lammpstrj."
    fp_info = fp + "cg_struct_info.pkl"

    # atom3ds = ans.atom3d_generator_from_aa(fp_info, fp_data, fp_trj, start=0, step=100000, frames=11, one_file=False)
    starts = range(120)
    time = 2
    atom3d_gt = ans.Atom3dGenerator(fp_data, fp_trj, one_file=False, info_path=fp_info)
    msd = ans.DynamicAnalyzer.mean_squared_displacement(atom3d_gt, starts, time, max_step=100)
    print(msd)
