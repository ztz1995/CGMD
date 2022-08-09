import pickle as pkl
import numpy as np

if __name__ == '__main__':

    # title = "aa"
    h_cut = 2.5
    # for i in [0, 1, 2, 3, 5, 7, 8, 9]:
    #     with open("data/%s_cut_%.1f_%d.pkl" % (title, h_cut, i), "rb") as f:
    #         h_dict = pkl.load(f)
    #     h_dict["h1"] = np.asarray(h_dict["h1"]) / 50 / 2
    #     h_dict["h2"] = np.asarray(h_dict["h2"]) / 50 / 2
    #     h_dict["hs"] = np.asarray(h_dict["hs"]) / 50 / 2
    #     with open("data/%s_cut_%.1f_%d.pkl" % (title, h_cut, i), "wb") as f:
    #         pkl.dump(h_dict, f)

    with open("data/%s_cut_%.1f_%d.pkl" % ("cgu", h_cut, 369), "rb") as f:
        h_dict = pkl.load(f)
    h_dict["h1"] = np.asarray(h_dict["h1"]) / 50 / 2
    h_dict["h2"] = np.asarray(h_dict["h2"]) / 50 / 2
    h_dict["hs"] = np.asarray(h_dict["hs"]) / 50 / 2
    with open("data/%s_cut_%.1f_%d.pkl" % ("cgu", h_cut, 369), "wb") as f:
        pkl.dump(h_dict, f)