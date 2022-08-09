import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt


def survey(results, category_names, ax):
    """
    Parameters
    ----------
    results : dict
        A mapping from question labels to a list of answers per category.
        It is assumed all lists contain the same number of entries and that
        it matches the length of *category_names*.
    category_names : list of str
        The category labels.
    """
    labels = list(results.keys())
    data = np.array(list(results.values()))
    data_cum = data.cumsum(axis=1)
    category_colors = plt.get_cmap('RdYlGn')(
        # category_colors = plt.get_cmap('RdYlBu')(
        # category_colors = plt.get_cmap('ocean')(
        # category_colors = plt.get_cmap('gray')(
        np.linspace(0.15, 0.75, data.shape[1] - 1)[::-1])
    # RdYlGn terrain

    ax.invert_yaxis()
    ax.xaxis.set_visible(False)
    ax.set_xlim(0, np.sum(data, axis=1).max())
    # labels = [0, 1, 2]
    for i, (colname, color) in enumerate(zip(category_names, category_colors)):
        widths = data[:, i]
        starts = data_cum[:, i] - widths
        ax.barh(labels, widths, left=starts, height=0.5,
                label=colname, color=color)
        xcenters = starts + widths / 2

        r, g, b, _ = color
        # text_color = 'white' if r * g * b < 0.3 else 'black'
        text_color = 'black'
        for y, (x, c) in enumerate(zip(xcenters, widths)):
            if c > 0.1:
                ax.text(x, y+0.03, "%.2f" % c, ha='center', va='center',
                        color=text_color, fontsize=10)
    ax.legend(ncol=len(category_names), bbox_to_anchor=(-0.02,-0.2, 1.04, 0.28),
              loc="lower left", fontsize=font_size - 5, mode="expand")
    plt.yticks(fontsize=font_size)
    return ax


if __name__ == '__main__':
    # plt.style.use(["science", "ieee", "no-latex"])
    plt.rc("font", family='Times New Roman')
    inch_cm = 1 / 2.54
    # inch_cm = 1
    fig = plt.figure(dpi=1000, figsize=(8 * inch_cm, 4 * inch_cm))
    ax1 = fig.add_subplot(111)
    font_size = 12

    plt.subplots_adjust(left=0.1, bottom=0.5, right=0.9, wspace=0.1)

    results = dict()
    sum_res = np.zeros(7)
    for i in [0, 1, 2, 3, 5, 7, 8, 9]:
        # for i in [0, 2, 3, 7, 8, 9]:
        with open("data/aa_u_5.0_%d.pkl" % i, "rb") as f:
            a = pkl.load(f)
            sum_res += a
        # results["AA_%d"%i] = a[:-2]
    # print(sum_res)
    results["AA"] = list(sum_res / 8)[:]

    title = "cgu"
    with open("data/%s_u_5.0_new.pkl" % title, "rb") as f:
        results["CGU"] = pkl.load(f)[:]

    title = "cg"
    with open("data/%s_u_5.0_new.pkl" % title, "rb") as f:
        results["CG"] = pkl.load(f)[:]

    category_names = [str(_) for _ in range(7)]

    survey(results, category_names, ax1)
    # plt.title("neighbor urea group number cut_off=5.0", fontdict={'fontsize': 18})
    plt.show()
