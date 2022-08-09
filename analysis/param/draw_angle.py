import matplotlib.pyplot as plt
import pandas as pd
import matplotlib

if __name__ == '__main__':
    plt.style.use(["science", "ieee", "no-latex"])
    plt.rc("font", family='Times New Roman')

    labels = ["AA", "CGU", "CG"]
    names = ["aa2cg", "cgu2cg", "cg"]
    plt.figure(dpi=1000)

    for i in range(3):
        l = labels[i]
        name = names[i]
        # data = pd.read_csv("data/Angle_Ph_U_Ph_%s.csv" % l.lower())
        data = pd.read_csv("data/10ns/%s/Angle_Ph_U_Ph_%s.csv" % (name, name))
        x = data.x
        y = data.gr
        # y /= y.sum() * (x[1] - x[0])
        # x = x[::5]
        # y = y[::5]
        plt.plot(x, y, label=l)
    plt.ylabel("Probability", fontsize=14)
    plt.xlabel("Angle (degrees)", fontsize=14)
    plt.legend(loc="upper left")
    plt.show()
