from pyMD import collective_structure_class as csc
from pyMD.file_parser import LmpParser


if __name__ == '__main__':
    atom3d = csc.Atom3D()
    pattern = ["TO(2)"] * 7
    param = {"Bond":{("TO(2)", "TO(2)"):[1.20593084e-04, 2.86654111e+01, 1.82028604e+04]}, "Angle":{}}
    atom3d.create_from_pattern(pattern, 200, param, 0.01, anglecheck=False)
    print(atom3d.lattice_parameter)
    LmpParser.create_data_file(atom3d, "test.data")
