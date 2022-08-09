import numpy as np
import pandas as pd
from pyMD import analysis as ans
from collections import OrderedDict
import os
import matplotlib.pyplot as plt
import pickle as pkl
from multiprocessing import Pool

if __name__ == '__main__':
    fp = "/home/centos/ztz/1blk_50/%s_8blk_250/equi/" % "cg"
    fp_data = fp + "strain.data"
    fp_trj = fp + "trj/strain.lammpstrj."
    atom3d_gt = ans.Atom3dGenSegCMV2(fp_data, fp_trj, one_file=False, pad=9)


