from pyMD import collective_structure_class as csc
from pyMD.file_parser import LmpParser


if __name__ == '__main__':
    data_fp = "/home/centos/Projects/MD_lmp/data/20_03_05_lr_0.05_old_fene/iter_01/cal/iter_01.data"
    trj_fp = "/home/centos/Projects/MD_lmp/data/20_03_05_lr_0.05_old_fene/iter_01/cal/iter_01_300K.lammpstrj."
    out_fp = "cg_fgh.data"

    # fp = "/home/centos/work/7blk_500chain/1E-6_fene/"
    # data_fp = fp + "strain.data"
    # trj_fp = fp + "strain.lammpstrj.s1."
    # out_fp = "/home/centos/work/fg_7blk_500chain/1E-6/strain_2.data"
    parser = LmpParser()
    atom3d = parser.load_data_file(data_fp, "cg_fgh")
    parser.renew_coordinate_file(atom3d, trj_fp, 0)
    atom3d.convert_to_fg(h=True)
    atom3d.cal_charge()
    # atom3d.convert_to_in_cell()
    parser.create_data_file(atom3d, out_fp, q=True)

