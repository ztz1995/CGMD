import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt


def survey(results, category_names):
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
    # category_colors = plt.get_cmap('RdYlGn')(
    category_colors = plt.get_cmap('gray')(
        np.linspace(0.15, 0.75, data.shape[1]-1)[::-1])
    # RdYlGn terrain
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.invert_yaxis()
    ax.xaxis.set_visible(False)
    ax.set_xlim(0, np.sum(data, axis=1).max())

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
            if c > 0.05:
                ax.text(x, y, "%.2f" % c, ha='center', va='center',
                        color=text_color, fontsize=14)
    ax.legend(ncol=len(category_names), bbox_to_anchor=(0., -0.18, 1., -0.18),
              loc='lower left', fontsize=14, mode="expand")
    plt.yticks(fontsize=14)
    return fig, ax


if __name__ == '__main__':
    results = dict()
    sum_res = np.zeros(7)
    # for i in [0, 1, 2, 3, 5, 7, 8, 9]:
    for i in [0, 2, 3, 7, 8, 9]:
        with open("data/aa_u_5.0_%d.pkl" % i, "rb") as f:
            a = pkl.load(f)
            sum_res += a
        # results["AA_%d"%i] = a[:-2]
    # print(sum_res)
    results["AA-1blk"] = list(sum_res / 6)[:]

    i = 369
    with open("data/cgu_u_5.0_%d.pkl" % i, "rb") as f:
        results["CGU-1blk"] = pkl.load(f)[:]

    # i = 0
    # with open("data/cgu_7blk_u_5.0_%d.pkl" % i, "rb") as f:
    #     results["CGU-7blk"] = pkl.load(f)[:-1]

    i = 0
    with open("data/cgu_7blk_u_5.0_%d.pkl" % i, "rb") as f:
        # results["CGU-7blk"] = pkl.load(f)[:]
        results["CGU-8blk"] = pkl.load(f)[:]

    i = 52
    with open("data/cg_u_5.0_%d.pkl" % i, "rb") as f:
        results["CG-1blk"] = pkl.load(f)[:]

    i = 0
    with open("data/cg_7blk_u_5.0_%d.pkl" % i, "rb") as f:
        # results["CG-7blk"] = pkl.load(f)[:]
        results["CG-8blk"] = pkl.load(f)[:]

    category_names = [str(_) for _ in range(7)]

    survey(results, category_names)
    plt.title("neighbor urea group number cut_off=5.0", fontdict={'fontsize': 18})
    plt.show()
