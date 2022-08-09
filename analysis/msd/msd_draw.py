import pandas as pd
import matplotlib.pyplot as plt

if __name__ == '__main__':
    # plt.figure(figsize=(5, 4), dpi=300)
    # plt.axes(xscale="log", yscale="log")
    # plt.axes(xscale="log")

    df1 = pd.read_csv("aa_cg_msd.csv")
    df2 = pd.read_csv("cgu_msd.csv")
    df3 = pd.read_csv("cg_msd.csv")
    df4 = pd.read_csv("aa_cg_cm_msd.csv")
    df5 = pd.read_csv("cg_cm_msd.csv")
    df6 = pd.read_csv("cgu_cm_msd.csv")
    df7 = pd.read_csv("aa_seg_msd.csv")
    df8 = pd.read_csv("cg_seg_msd.csv")
    df9 = pd.read_csv("cgu_seg_msd.csv")

    # bg = 10
    # ed = 65
    # x = df1.time[bg:ed] / 1E6
    # cgu_h = df2["Ph"][bg:ed]/df1["Ph"][bg:ed]
    # cgu_s = df2["TO(2)"][bg:ed]/df1["TO(2)"][bg:ed]
    # # cg_h = df3["Ph"][bg:ed]/df1["Ph"][bg:ed]
    # # cg_s = df3["TO(2)"][bg:ed]/df1["TO(2)"][bg:ed]
    # plt.plot(x, cgu_h, label="cgu_h")
    # plt.plot(x, cgu_s, label="cgu_s")
    # # plt.plot(x, cg_h, label="cg_h")
    # # plt.plot(x, cg_s, label="cg_s")
    # plt.legend()

    # plt.plot(df2.time/1E6, df2["Ph"])
    # plt.plot(df2.time/1E6, df2["TO(2)"])
    # plt.plot(df2.time/1E6, df2["c1"])
    # plt.plot(df1.time/1E6, df1["Ph"])
    # plt.plot(df1.time/1E6, df1["TO(2)"])
    # plt.plot(df1.time/1E6, df1["U"])
    bg = 10
    ed = 68
    # ed = 70
    # plt.plot(df4.time[bg:ed]/1E6, df4["cm"][bg:ed], label="aa")
    # plt.plot(df5.time[bg:ed]/1E6, df5["cm"][bg:ed], label="cg")
    # plt.plot(df6.time[bg:ed]/1E6, df6["cm"][bg:ed], label="cgu")


    # bg = 10
    # ed = 65
    # t = "cm"
    fig, (ax0, ax1) = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(8, 4), dpi=300)

    t = "hard"
    x = df9.time[bg:ed] / 1E6
    # cgu = df6[t][bg:ed]/df4[t][bg:ed]
    # cg = df5[t][bg:ed]/df4[t][bg:ed]
    cgu = df9[t][bg:ed] / df7[t][bg:ed]
    cg = df8[t][bg:ed] / df7[t][bg:ed]
    ax0.plot(x, cgu, label="cgu/aa")
    ax0.plot(x, cg, label="cg/aa")
    ax0.set_title("Hard Segment",fontsize=14)
    ax0.set_xscale("log")
    ax0.set_yscale("log")
    ax0.set_xlabel("t/ns",fontsize=12)
    ax0.set_ylabel("msd",fontsize=12)

    t = "soft"
    x = df9.time[bg:ed] / 1E6
    # cgu = df6[t][bg:ed]/df4[t][bg:ed]
    # cg = df5[t][bg:ed]/df4[t][bg:ed]
    cgu = df9[t][bg:ed] / df7[t][bg:ed]
    cg = df8[t][bg:ed] / df7[t][bg:ed]
    ax1.plot(x, cgu, label="cgu/aa")
    ax1.plot(x, cg, label="cg/aa")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_title("Soft Segment",fontsize=14)
    ax1.set_xlabel("t/ns",fontsize=12)

    # plt.tick_params(labelsize=18)
    plt.legend()
    # fig.suptitle("Mean-Squared Displacement",fontsize=18)
    plt.tight_layout()
    # plt.show()

    # df = pd.read_csv("aa_cgu_msd.csv")
    # plt.axes(xscale="log", yscale="log")
    # plt.plot(df.time/1E6, df["c1"])
    # plt.plot(df.time/1E6, df["TO(2)"])
    # plt.show()

    # # plt.axes(xscale="log", yscale="log")
    # plt.plot(df.time/1E6, df["Ph"])
    # plt.plot(df.time/1E6, df["TO(2)"])
    # # plt.show()
    #
    # # plt.axes(xscale="log", yscale="log")
    # # plt.plot(df.time/1E6, df["Ph"])
    # plt.plot(df.time/1E6, df["TO(2)"])
    plt.show()
