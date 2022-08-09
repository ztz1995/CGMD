from pyMD import collective_structure_class as csc
from pyMD import basic_structure_class as bsc
from pyMD.file_parser import LmpParser


if __name__ == '__main__':
    # data_fp = "/home/centos/Projects/MD_lmp/data/20_03_05_lr_0.05_old_fene/iter_01/cal/iter_01.data"
    # trj_fp = "/home/centos/Projects/MD_lmp/data/20_03_05_lr_0.05_old_fene/iter_01/cal/iter_01_300K.lammpstrj."

    # # fp = "/home/centos/work/7blk_500chain/1E-6_fene/"
    # fp = "/home/centos/work/200-125_1E-7_new/"
    # data_fp = fp + "strain.data"
    # trj_fp = fp + "strain.lammpstrj.s1."
    # out_fp = "/home/centos/work/fgh_7blk_200chain/strain.data"
    # # out_fp2 = "/home/centos/work/fgh_7blk_500chain/1E-6/strain_1.data"
    # parser = LmpParser()
    # atom3d = parser.load_data_file(data_fp, "cg_fgh")
    # parser.renew_coordinate_file(atom3d, trj_fp, 0)
    
    # # # atom3d.convert_to_in_cell()
    # # # parser.create_data_file(atom3d, out_fp, q=True)
    # # atom3d = parser.load_data_file(out_fp, q=False)
    # parser.create_data_file(atom3d, out_fp)

    # fp = "/home/centos/ztz/dpd/cg/dpd/1blk_50/init/"
    # data_fp = fp + "data.1blk_50"
    # trj_fp = fp + "trj/300K.lammpstrj."
    # out_fp = "/home/centos/ztz/dpd/cg/dpd/1blk_50/fgh/data.1blk_50"

    fp = "/home/centos/ztz/dpd/cgu/8blk_250/1E-6_1/"
    data_fp = fp + "data.8blk_250"
    trj_fp = fp + "init.8blk_250.lammpstrj."
    out_fp = "/home/centos/ztz/dpd/fgh/8blk_250/1E-6/data.8blk_250"

    parser = LmpParser()
    atom3d = parser.load_data_file(data_fp, "fg", map_dict={"Me":"c1"})
    parser.renew_coordinate_file(atom3d, trj_fp, 180000)

    atom3d.convert_to_fg(h=True)
    # atom3d.sort_atom_bond_id()
    # ma = atom3d.get_max_atom_id()+1
    # mb = atom3d.get_max_bond_id()+1
    # atom3d.append_element(bsc.Atom(ma, _type="TO(3)"))
    # ma+=1
    # atom3d.append_element(bsc.Atom(ma, _type="Es"))
    # ma+=1
    # atom3d.append_element(bsc.Bond(mb, ma-1, ma-2, _type=("Es", "TO(3)")))
    # mb+=1
    # atom3d.append_element(bsc.Atom(ma, _type="c2"))
    # ma+=1
    # atom3d.append_element(bsc.Bond(mb, ma-2, ma-1, _type=("Es", "c2")))
    # mb += 1
    # atom3d.append_element(bsc.Atom(ma, _type="c2"))
    # atom3d.append_element(bsc.Bond(mb, ma, ma-1, _type=("c2", "c2")))
    atom3d.cg_to_fgh()

    parser.create_data_file(atom3d, out_fp)