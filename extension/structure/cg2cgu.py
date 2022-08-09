from pyMD.file_parser import LmpParser
import os


if __name__ == '__main__':
    parser = LmpParser()
    # cg_atom3d = parser.load_data_file("/home/centos/work/change_density/cg_1blk_200chain/init/strain.data")
    # fp1 = "/home/centos/work/fix/cg_1blk_204chain/init2/"
    # fp1 = "/home/centos/work/1blk_50/cg_7blk_300chain/init/"
    # fp1 = "/home/centos/work/1blk_50/cg_7blk_300chain_wallx/init/"
    fp1 = "/home/centos/ztz/stress_new/8blk_2000/dpd2/"
    cg_atom3d = parser.load_data_file(fp1 + "data.8blk_2000")
    cg_atom3d.renew_coordinate_file(fp1 + "trj/dpd.8blk_2000.lammpstrj.7300000")

    # cg_atom3d.cg_to_cgu(u_orien=[0., 0., 1])
    cg_atom3d.delete_bond_type(("U", "U"))
    cg_atom3d.cg_to_cgu()
    # fp2 = "/home/centos/work/fix/cgu_1blk_204chain/init2_vertical_SH6S/"
    # fp2 = "/home/centos/work/1blk_50/cgu_7blk_300chain/init/"
    # fp2 = "/home/centos/work/1blk_50/cgu_7blk_300chain_wallx/init/"
    fp2 = "/home/centos/ztz/stress_new/8blk_2000/cgu/init/"
    os.makedirs(fp2, exist_ok=True)
    # cg_atom3d.box_h[0] = 175.
    # cg_atom3d.box_l[0] = 0.
    # cg_atom3d.lattice_parameter[0] = 175.
    LmpParser.create_data_file(cg_atom3d, fp2 + "data.8blk_2000", q=False, improper=True)
