import matplotlib.pyplot as plt

plt.style.use(["science", "ieee", "no-latex"])
plt.rc("font", family='Times New Roman')


# pts[[3, 14]] += 1.8#将索引为3个和14的元素加1.8处理成两个离散点
#
#
# fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, dpi=100)
#
#
# ax1.plot(pts)
# ax2.plot(pts)
#
# ax1.set_ylim(.4, 2.)  # 子图1设置y轴范围，只显示部分图
# ax2.set_ylim(0, .28)  # 子图2设置y轴范围，只显示部分图
#
#
# ax1.spines['bottom'].set_visible(False)#关闭子图1中底部脊
# ax2.spines['top'].set_visible(False)##关闭子图2中顶部脊
# ax2.set_xticks(range(0,31,1))
#
#
# d = .85  #设置倾斜度
# #绘制断裂处的标记
# kwargs = dict(marker=[(-1, -d), (1, d)], markersize=15,
#               linestyle='none', color='r', mec='r', mew=1, clip_on=False)
# ax1.plot([0, 1], [0, 0],transform=ax1.transAxes, **kwargs)
# ax2.plot([0, 1], [1, 1], transform=ax2.transAxes, **kwargs)
#
# plt.tight_layout()
# plt.show()