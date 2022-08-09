from pyMD.file_parser import LmpParser
import pickle as pkl


if __name__ == '__main__':
    for i in range(10):
        fp = "/home/centos/model/aa/1blk_50/%d/" % i
        atom3d_aa = LmpParser.load_data_file(fp + "1blk_50.data")
        with open(fp + "cgu_struct_info.pkl", "rb") as file:
            info = pkl.load(file)
        # id_type
        for group_id, group in info["id_ids"].items():
            for atom_id in group:
                atom3d_aa.Atoms[atom_id].type = info["id_type"][group_id]

        for atom in atom3d_aa.Atoms.values():
            # if atom.type[0] not in ["T", "E", "U", "P", "c", "o", "n", "h", "M"]:
            #     print(atom)
            if atom.type[0].isdigit():
                atom.type = "Me"
                # connect = atom3d_aa.get_connected_atom_id(atom.id).pop()
                # print(atom3d_aa.Atoms[connect])
        atom3d_aa.recreate_atom_type_dict()
        LmpParser.create_data_file(atom3d_aa, fp + "1blk_50_cgu.data")