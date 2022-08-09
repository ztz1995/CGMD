# -*- coding: UTF-8 -*-
import sys
import os

project_path = "/home/centos/Projects/CGMD/"
sys.path.append(project_path)
import pyMD.monte_carlo as mc
from pyMD.file_parser import LmpParser
from pyMD import analysis as ans
import pickle as pkl
import numpy as np


if __name__ == '__main__':

    # name = "para18"
    # name = "tilt68"
    name = "tilt8"
    # fp = "/home/centos/ztz/stress/line_test/para18/init/"
    fp = "/home/centos/ztz/stress/line_test/tilt8/init/"
    atom3d = LmpParser.load_data_file(fp + "data.%s" % name)
    # atom3d.renew_coordinate_file("/home/centos/ztz/stress/line_test/para18_intra/para18_8/1E-6/trj/para18_8.lammpstrj.s1.0")
    # atom3d.renew_coordinate_file(fp + "trj/init.tilt8.lammpstrj.420000")

    print(atom3d.lattice_parameter)
    chains = 160
    # chains = 16000
    # out_fp = "/home/centos/ztz/stress/line_test/para18_intra/"
    # out_fp = "/home/centos/ztz/stress/line_test/para18_inter/"
    # out_fp = "/home/centos/ztz/stress/line_test/para18_random3/"
    out_fp = "/home/centos/ztz/stress/line_test/tilt8_random/"

    # chain_model = mc.FreelyJointChain(5.13)
    chain_model = mc.TOChain()

    def exclude_func_para(lattice_parameter, linesy, linesz):
        def exclude(co):
            dz = lattice_parameter[2]/linesz
            dy = lattice_parameter[1]/linesy
            co = co % lattice_parameter
            z = (co[2] - dz/2) % dz
            y = co[1] % dy
            judge_z = z < 3 or z > dz - 3
            judge_y = (dy-22)/2 < y < (dy - (dy-22)/2)
            return judge_y and judge_z
        return exclude

    def exclude_func_tilt(lattice_parameter, linesy, linesz):
        def exclude(co):
            co = co % lattice_parameter
            xz = (lattice_parameter[0] - co[0] + co[2])*0.7071
            dz = np.sqrt(np.square(lattice_parameter[0])+np.square(lattice_parameter[2]))/linesz
            xz = xz % dz
            # print("xz", xz)
            # print("dz", dz)
            judge_z = (dz-21.01941)/2 < xz < (dz - (dz-21.01941)/2)
            y = co[1]
            dy = lattice_parameter[1]/linesy
            y = y % dy
            judge_y = -3 < y-dy/2 < 3
            return judge_y and judge_z
        return exclude


    chain_creator = mc.ChainCreator(chain_model, 5.13, exclude_func=exclude_func_tilt(atom3d.lattice_parameter, 4, 4))
    # chain_creator = mc.ChainCreator(chain_model, 5.13)
    # judger = ans.SJPara(fp + "segment.%s" % name, "x")
    # for i in [4, 2]:
    for i in [8]:
        # test_box = mc.HardBox(chain_model, chain_creator, judger, intra_p=1.)
        # test_box = mc.HardBox(chain_model, chain_creator, judger, intra_p=0.)
        # test_box = mc.HardBox(chain_model, chain_creator, judger, intra_p=0.5)
        test_box = mc.HardBox(chain_model, chain_creator)
        test_box.initialize_atom3d(atom3d, repeat_num=i)

        mc_s = mc.MonteCarlo(test_box)
        mc_s.simulation(298, 1000000)

        # with open("%s_%d.pkl" % (name, i), "wb") as file:
        #     pkl.dump(mc_s, file)

        # with open("%dblk_%d.pkl" % (i, chains // i), "rb") as file:
        #     mc_s = pkl.load(file)
        # mc_s.simulation(298, 1000000)
        # with open("%dblk_%d.pkl" % (i, chains // i), "wb") as file:
        #     pkl.dump(mc_s, file)
        ofp = out_fp + "%s_%d/init/" % (name, i)
        os.makedirs(ofp, exist_ok=True)
        mc_s.state.create_data_file(ofp + "data.%s_%d" % (name, i))
