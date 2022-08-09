from pyMD import analysis as ans
from pyMD import file_parser as fps
import  pickle as pkl

if __name__ == '__main__':

    fp = "/home/centos/work/7blk_500chain/5E-8_fene/"
    data_fp = fp + "strain.data"
    trj_fp = fp + "strain.lammpstrj.s1."
    step = 400000
    atom3ds = ans.atom3d_generator(data_fp, trj_fp, start=10 * step, step=step, frames=20, _type="cg")

    with open("5E-8_fene/result.pkl", "rb") as f:
        result_dict = pkl.load(f)

    for i, atom3d in enumerate(atom3ds):
        strain = int((i+10) * step * 5 * 5E-8 * 100 + 0.1)
        print(strain)
        v = ans.StaticAnalyzer.grid_hard_domain(atom3d, grid=0.5, r=3., out_path="5E-8_fene/strain_%d.data" % int(strain))
        print(v)
        result_dict[strain] = v

    with open("5E-8_fene/result.pkl", "wb") as f:
        pkl.dump(result_dict, f)
