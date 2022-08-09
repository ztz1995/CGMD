from pyMD.ibi import IBI
from pyMD import analysis as ans


if __name__ == '__main__':
    label = "cgu"
    print(label)
    for i in range(2):
        print(i)
        fp = "/home/centos/work/1blk_50/cgu_1blk_50chain/%d/" % i
        fp_data = fp + "1blk_50.data"
        fp_trj = fp + "trj/1blk_50.lammpstrj."
        atom3ds = ans.Atom3dGenerator(fp_data, fp_trj, one_file=False, pad=9)

        ibi = IBI("Cal %s, 1blk_50" % label, atom3ds())
        out_path = "data/10ns/%s/%d/" % (label, i)
        ibi.cal_gr_new(atom3ds, out_path, range(10000000, 11000000, 1000), num_parallel=24, name=label)
        del ibi
        del atom3ds

    label = "cg"
    print(label)
    for i in range(8,10):
        print(i)
        fp = "/home/centos/work/1blk_50/cg_1blk_50chain/%d/" % i
        fp_data = fp + "1blk_50.data"
        fp_trj = fp + "trj/1blk_50.lammpstrj."
        atom3ds = ans.Atom3dGenerator(fp_data, fp_trj, one_file=False, pad=9)

        ibi = IBI("Cal %s, 1blk_50" % label, atom3ds())
        out_path = "data/10ns/%s/%d/" % (label, i)
        ibi.cal_gr_new(atom3ds, out_path, range(10000000, 11000000, 1000), num_parallel=24, name=label)

    label = "cgu2cg"
    print(label)
    for i in range(2):
        print(i)
        fp = "/home/centos/work/1blk_50/cgu_1blk_50chain/%d/" % i
        fp_data = fp + "1blk_50.data"
        fp_trj = fp + "trj/1blk_50.lammpstrj."
        atom3ds = ans.Atom3dGenCGU2CG(fp_data, fp_trj, one_file=False, pad=9)

        ibi = IBI("Cal %s, 1blk_50" % label, atom3ds())
        out_path = "data/10ns/%s/%d/" % (label, i)
        ibi.cal_gr_new(atom3ds, out_path, range(10000000, 11000000, 1000), num_parallel=24, name=label)
        del ibi
        del atom3ds

    # label = "cgu2cg"
    # fp = "/home/centos/Projects/CGMD/data/20200722_lr=0.05_1blk_50_cgu_ex_hard/iter_353/cal/"
    # fp_data = fp + "1blk_40.data"
    # fp_trj = fp + "300K.lammpstrj."
    # atom3ds = ans.Atom3dGenCGU2CG(fp_data, fp_trj, one_file=False, pad=0)
    # ibi = IBI("Cal %s, 1blk_50" % label, atom3ds())
    # out_path = "data/%s/" % label
    # # ibi.cal_gr_new(atom3ds, out_path, range(0, 500, 500), num_parallel=1, name=label)
    # ibi.cal_gr_new(atom3ds, out_path, range(0, 500000, 500), num_parallel=24, name=label)
    # del ibi
    # del atom3ds

    # label = "cgu"
    # fp = "/home/centos/Projects/CGMD/data/20200722_lr=0.05_1blk_50_cgu_ex_hard/iter_353/cal/"
    # fp_data = fp + "1blk_40.data"
    # fp_trj = fp + "300K.lammpstrj."
    # atom3ds = ans.Atom3dGenerator(fp_data, fp_trj, one_file=False, pad=0)
    # # atom3d = atom3ds(0, fake=False)
    # # from pyMD.file_parser import LmpParser
    # # LmpParser.create_data_file(atom3d, "/home/centos/test_cgu.data", improper=True)
    # ibi = IBI("Cal %s, 1blk_50" % label, atom3ds())
    # out_path = "data/%s/" % label
    # # ibi.cal_gr_new(atom3ds, out_path, range(0, 500000, 500), num_parallel=24, name=label)
    # ibi.cal_gr_new(atom3ds, out_path, [0], num_parallel=1, name=label)
    # del ibi
    # del atom3ds
    #
    # label = "cg"
    # fp = "/home/centos/Projects/CGMD/data/20200720_lr=0.05_1blk_50_cg_ex_hard/iter_52/cal/"
    # fp_data = fp + "1blk_40.data"
    # fp_trj = fp + "300K.lammpstrj."
    # atom3ds = ans.Atom3dGenerator(fp_data, fp_trj, one_file=False, pad=0)
    # ibi = IBI("Cal %s, 1blk_50" % label, atom3ds())
    # out_path = "data/%s/" % label
    # ibi.cal_gr_new(atom3ds, out_path, range(0, 500000, 500), num_parallel=24, name=label)

    # fp = "/home/centos/Projects/CGMD/data/20200722_lr=0.05_1blk_50_cgu_ex_hard/iter_353/cal/"
    # fp_data = fp + "1blk_40.data"
    # fp_trj = fp + "300K.lammpstrj."
    # atom3ds = ans.Atom3dGenCGU2CG(fp_data, fp_trj, one_file=False, pad=0)
    # cgu = atom3ds(0, fake=False)
    # print(cgu)
    #
    # fp = "/home/centos/Projects/CGMD/data/20200722_lr=0.05_1blk_50_cgu_ex_hard/iter_353/cal/"
    # fp_data = fp + "1blk_40.data"
    # fp_trj = fp + "300K.lammpstrj."
    # atom3ds = ans.Atom3dGenerator(fp_data, fp_trj, one_file=False, pad=0)
    # atom3ds(0)
    # atom3ds.atom3d.convert_to_default()
    # print(atom3ds.atom3d)
