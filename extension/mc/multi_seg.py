# -*- coding: UTF-8 -*-
import sys
import os

project_path = "/home/centos/Projects/CGMD/"
sys.path.append(project_path)
import pyMD.monte_carlo as mc
from pyMD.file_parser import LmpParser
import pickle as pkl


if __name__ == '__main__':
    # atom3d = LmpParser.load_data_file("/home/centos/work/SH6S/cgu_1blk_2000chain/init/strain.data")
    # atom3d.renew_coordinate_file("/home/centos/work/SH6S/cgu_1blk_2000chain/1E-6/strain.lammpstrj.s1.0")

    # atom3d = LmpParser.load_data_file("/home/centos/work/1blk_50/cg_1blk_2000chain/init/strain.data")
    # atom3d.renew_coordinate_file("/home/centos/work/1blk_50/cg_1blk_2000chain/init/record.lammpstrj.10300000")

    fp = "/home/centos/ztz/stress/cgu/1blk_2000chain/SH6S_init/"
    atom3d = LmpParser.load_data_file(fp + "1blk_2000.data")
    atom3d.renew_coordinate_file(fp + "trj/1blk_2000.lammpstrj.011000000")

    # fp = "/home/centos/ztz/stress/cgu/1blk_16000chain/init_SH6S/"
    # atom3d = LmpParser.load_data_file(fp + "1blk_16000.data")
    # atom3d.renew_coordinate_file(fp + "trj/1blk_16000.lammpstrj.012000000")

    print(atom3d.lattice_parameter)
    chains = 2000
    # chains = 16000
    out_fp = "/home/centos/ztz/stress/cgu/1blk_2000chain/"
    # out_fp = "/home/centos/ztz/stress/cgu/1blk_16000chain/"
    # atom3d = LmpParser.load_data_file("/home/centos/work/1blk_50/cgu_1blk_160chain/init/1blk_160.data")
    # atom3d.renew_coordinate_file(
    #     "/home/centos/work/1blk_50/cgu_1blk_160chain/init/trj/1blk_160.lammpstrj.010000000")
    # chains = 160
    chain_model = mc.TOChain()
    chain_creator = mc.ChainCreator(chain_model, 5.13)
    # for i in [4, 2]:
    for i in [4]:
        test_box = mc.HardBox(chain_model, chain_creator)
        test_box.initialize_atom3d(atom3d, repeat_num=i)

        mc_s = mc.MonteCarlo(test_box)
        mc_s.simulation(298, 2000000)

        with open("%dblk_%d.pkl" % (i, chains // i), "wb") as file:
            pkl.dump(mc_s, file)

        # with open("%dblk_%d.pkl" % (i, chains // i), "rb") as file:
        #     mc_s = pkl.load(file)
        # mc_s.simulation(298, 1000000)
        # with open("%dblk_%d.pkl" % (i, chains // i), "wb") as file:
        #     pkl.dump(mc_s, file)
        ofp = out_fp + "8blk_250_new/%dblk_%d_init/" % (i, chains // i)
        os.makedirs(ofp, exist_ok=True)
        mc_s.state.create_data_file(ofp + "data.%dblk_%d" % (i, chains // i))



