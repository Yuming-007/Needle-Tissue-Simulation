"""单点针尖 vs 平底圆盘针尖的对比图（只读已有数据，不跑仿真）。

三幅小图（各自独立的纵轴）：
  (a) 小压深刚度随网格尺寸 h 的变化：单点 = 第 1 步 P9（results/step1/t1_mesh_dependence.json，压深 1 mm 割线刚度）；
      圆盘 = V3-3（results/step1_5/v3_3_dense.json，a_eff = 10 mm，点距 s = h，线性小压深），虚线 = Sneddon 解析解 k_S。
  (b) 落点依赖：6 个落点的刚度相对"对准顶点"时的变化（V2，results/step1_5/v2_placement.json），h = 10、5 mm。
  (c) 大压深下单元翻转的压深随 h 的变化（V4，results/step1_5/v4_large_depth_{A,B}.json；圆盘 a_eff = 10 mm）。
颜色：参考调色板的类别色第 1、2 位（蓝 = 圆盘，橙 = 单点），并用不同标记形状 + 直接标注区分，不只靠颜色。
运行：python3 plot_point_vs_disk.py → results/step1_5/point_vs_disk.png
"""
import sys
sys.dont_write_bytecode = True
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

for f in ["/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"]:
    if os.path.exists(f):
        font_manager.fontManager.addfont(f)
        plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name()
plt.rcParams["axes.unicode_minus"] = False

R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../results")
DISK, POINT = "#2a78d6", "#eb6834"
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
KS = 2 * 0.010 * 5000 / (1 - 0.45 ** 2)


def style(ax, title):
    ax.set_facecolor(SURF)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    for s in ["left", "bottom"]:
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    ax.set_title(title, loc="left", fontsize=11, color=INK, pad=10)


fig, axs = plt.subplots(1, 3, figsize=(15, 4.8), facecolor=SURF)

# (a) 刚度 vs h
p9 = json.load(open(os.path.join(R, "step1/t1_mesh_dependence.json")))
pt = sorted({(r["h_mm"], r["k_secant"]) for r in p9 if abs(r["depth_mm"] - 1.0) < 1e-9})
d3 = json.load(open(os.path.join(R, "step1_5/v3_3_dense.json")))
dk = sorted([(r["h_mm"], r["k_over_kS"] * KS) for r in d3 if abs(r["a_mm"] - 10) < 1e-9 and abs(r["s_over_h"] - 1) < 1e-9])
ax = axs[0]
style(ax, "(a) 网格加细时刚度是否收敛")
ax.plot([x for x, _ in pt], [y for _, y in pt], "-o", color=POINT, lw=2, ms=8, mec=SURF, mew=1.5)
ax.plot([x for x, _ in dk], [y for _, y in dk], "-s", color=DISK, lw=2, ms=8, mec=SURF, mew=1.5)
ax.axhline(KS, color=INK2, lw=1.2, ls="--")
ax.text(0.4, KS + 5, "解析解（Sneddon）", va="bottom", fontsize=9, color=INK2)
ax.text(6.3, 12, "单点：越细越软，趋向 0（不收敛）", fontsize=9, color=INK)
ax.text(0.4, 208, "圆盘（半径 10 mm）：越细越接近解析解", fontsize=9, color=INK)
ax.set_xlabel("单元尺寸 h（mm，越往左越细）", color=INK2, fontsize=9)
ax.set_ylabel("小压深刚度（N/m）", color=INK2, fontsize=9)
ax.set_xlim(0, 13.5); ax.set_ylim(0, 230)

# (b) 落点依赖
v2 = json.load(open(os.path.join(R, "step1_5/v2_placement.json")))
ax = axs[1]
style(ax, "(b) 针尖落在不同位置时刚度变化多少")
groups = [("单点\nh = 10", "point", 10.0, POINT, "o"), ("圆盘\nh = 10", "disk", 10.0, DISK, "s"),
          ("单点\nh = 5", "point", 5.0, POINT, "o"), ("圆盘\nh = 5", "disk", 5.0, DISK, "s")]
for i, (lab, tip, h, c, m) in enumerate(groups):
    ks = [r["k"] for r in v2 if r["tip"] == tip and abs(r["h_mm"] - h) < 1e-9]
    rel = 100 * (np.array(ks) / ks[0] - 1)
    ax.plot([i, i], [rel.min(), rel.max()], color=c, lw=2, alpha=0.5)
    ax.scatter(np.full(len(rel), i), rel, color=c, marker=m, s=60, edgecolor=SURF, linewidth=1.5, zorder=3)
    ax.text(i + 0.12, rel.max(), f"最大差 {ks and (max(ks)/min(ks)-1)*100:.0f}%", fontsize=9, color=INK, va="bottom")
ax.set_xticks(range(4)); ax.set_xticklabels([g[0] for g in groups], fontsize=9, color=INK2)
ax.set_ylabel("相对\u201c对准网格顶点\u201d时的刚度变化（%）", color=INK2, fontsize=9)
ax.set_xlim(-0.5, 3.9); ax.set_ylim(-5, 80)

# (c) 翻转压深 vs h
v4 = json.load(open(os.path.join(R, "step1_5/v4_large_depth_A.json"))) + json.load(open(os.path.join(R, "step1_5/v4_large_depth_B.json")))
ip = sorted([(r["h_mm"], r["dinv_mm"]) for r in v4 if r["tip"] == "point" and r["v_mm_s"] == 5])
idk = sorted([(r["h_mm"], r["dinv_mm"]) for r in v4 if r["tip"] == "disk" and abs(r["a_mm"] - 10) < 1e-9 and r["v_mm_s"] == 5])
ax = axs[2]
style(ax, "(c) 压多深时网格单元被压坏（翻转）")
ax.plot([x for x, _ in ip], [y for _, y in ip], "-o", color=POINT, lw=2, ms=8, mec=SURF, mew=1.5)
ax.plot([x for x, _ in idk], [y for _, y in idk], "-s", color=DISK, lw=2, ms=8, mec=SURF, mew=1.5)
ax.text(5.4, 4.6, "单点：约 1.5 倍单元尺寸，越细越早坏", fontsize=9, color=INK)
ax.text(2.3, 16.3, "圆盘：约 1.1–1.2 倍圆盘半径，基本不随网格变化", fontsize=9, color=INK)
ax.set_xlabel("单元尺寸 h（mm，越往左越细）", color=INK2, fontsize=9)
ax.set_ylabel("单元翻转时的压深（mm）", color=INK2, fontsize=9)
ax.set_xlim(2, 12); ax.set_ylim(0, 20)

fig.suptitle("单点针尖 vs 平底圆盘针尖（规则网格测试台，E = 5 kPa，ν = 0.45）", fontsize=13, color=INK, x=0.01, ha="left")
plt.tight_layout(rect=[0, 0, 1, 0.94])
out = os.path.join(R, "step1_5/point_vs_disk.png")
plt.savefig(out, dpi=130, facecolor=SURF)
print(out)
