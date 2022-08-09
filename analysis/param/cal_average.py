import pandas as pd
import os


if __name__ == '__main__':

    # label = "cgu"
    label = "cgu2cg"
    # label = "cg"
    path = "data/10ns/%s/" % label
    files = os.listdir(path + "0/")
    data_dict = dict()
    for fn in files:
        if fn[-3:] == "csv":
            data_dict[fn] = pd.read_csv(path + "0/" + fn)

    for i in range(1, 8):
        for fn in files:
            if fn[-3:] == "csv":
                new_data = pd.read_csv(path + "%d/" % i + fn)
                data_dict[fn]["y"] += new_data["y"]
                data_dict[fn]["gr"] += new_data["gr"]

    for fn in data_dict:
        data_dict[fn]["y"] /= 8
        data_dict[fn]["gr"] /= 8
        data_dict[fn].to_csv(path + fn, index=False)

    # label = "cg"
    # path = "data/%s/" % label
    # files = os.listdir(path + "0/")
    # data_dict = dict()
    # for fn in files:
    #     if fn[-3:] == "csv":
    #         data_dict[fn] = pd.read_csv(path + "0/" + fn)
    #
    # for i in range(1, 8):
    #     for fn in files:
    #         if fn[-3:] == "csv":
    #             new_data = pd.read_csv(path + "%d/" % i + fn)
    #             data_dict[fn]["y"] += new_data["y"]
    #             data_dict[fn]["gr"] += new_data["gr"]
    #
    # for fn in data_dict:
    #     data_dict[fn]["y"] /= 8
    #     data_dict[fn]["gr"] /= 8
    #     data_dict[fn].to_csv(path + fn, index=False)