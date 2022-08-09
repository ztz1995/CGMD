from pyMD import analysis as ans
from pyMD.file_parser import LmpParser

if __name__ == '__main__':
    parser = LmpParser()
    name = "8blk_2000"
    fp1 = "/home/centos/ztz/true_strain_cg/x/1E-6/"
    atom3d = parser.load_data_file(fp1 + "data.8blk_2000")
    visited = atom3d.set_block_id(begin_block=1)
    lines = ["# atom chunk id for different block\n", "%d\n" % len(atom3d.Atoms)]
    for a_id in visited:
        lines.append("%d %.1f\n" % (a_id, atom3d.Atoms[a_id].block))
    with open(fp1 + "block.%s" % name, "w") as f:
        f.writelines(lines)
