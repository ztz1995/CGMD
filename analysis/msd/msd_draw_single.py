import pandas as pd
import matplotlib.pyplot as plt

if __name__ == '__main__':
    plt.figure(figsize=(5, 4), dpi=300)
    plt.axes(xscale="log", yscale="log")

    df1 = pd.read_csv("aa_cg_msd.csv")
    df2 = pd.read_csv("cgu_msd.csv")
    df3 = pd.read_csv("cg_msd.csv")
    df4 = pd.read_csv("aa_cg_cm_msd.csv")
    df5 = pd.read_csv("cg_cm_msd.csv")
    df6 = pd.read_csv("cgu_cm_msd.csv")
    df7 = pd.read_csv("aa_seg_msd.csv")
    df8 = pd.read_csv("cg_seg_msd.csv")
    df9 = pd.read_csv("cgu_seg_msd.csv")

    bg = 10
    ed = 67

    t = "cm"
    x = df4.time[bg:ed] / 1E6
    # cgu = df6[t][bg:ed]/df4[t][bg:ed]
    # cg = df5[t][bg:ed]/df4[t][bg:ed]
    cgu = df6[t][bg:ed]
    aa = df4[t][bg:ed]
    cg = df5[t][bg:ed]

    plt.plot(x, cgu, label="CGU")
    plt.plot(x, aa, label="AA")
    plt.plot(x, cg, label="CG")
    plt.title("Mean-Squared Displacement",fontsize=14)
    plt.xlabel("t/ns",fontsize=12)
    plt.ylabel("msd",fontsize=12)

    plt.show()
