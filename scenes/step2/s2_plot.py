"""第 2 步的图：力–深度曲线（和 Okamura 2004 定性对比）、力–时间、顶面中心节点回弹、耦合点数。

读 results/step2/ 下 s2_puncture.py 的结果：s2_puncture_dt0.01.npz（主结果，SparseLDL + β = 0.1）和
s2_puncture_dt0.01_beta0.npz（不加阻尼的对照）。刺穿前的力–压深用 F = c0 + a1·d + a2·d² 拟合：a1·d + a2·d² 是 Okamura 2004 式 3 的形式；
加了阻尼后，针一碰到组织就有一个黏性力台阶（约 β·k·v），用常数 c0 表示（不加阻尼时 c0 ≈ 0）。
运行：python3 s2_plot.py → results/step2/s2_force_depth.png + 终端打印拟合系数
"""
import sys
sys.dont_write_bytecode = True
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../results/step2")


def load(name):
    return dict(np.load(os.path.join(OUT, name)))


def fit_quadratic(r):
    kp = int(np.argmax(r["n_cp"] > 0))
    d, F = r["depth"][:kp], r["F"][:kp, 2]
    m = d > 0
    A = np.c_[np.ones(m.sum()), d[m], d[m] ** 2]
    (c0, a1, a2), *_ = np.linalg.lstsq(A, F[m], rcond=None)
    res = F[m] - A @ [c0, a1, a2]
    r2 = 1 - (res ** 2).sum() / ((F[m] - F[m].mean()) ** 2).sum()
    return kp, c0, a1, a2, r2


if __name__ == "__main__":
    runs = [("SparseLDL, β = 0.1 (main)", load("s2_puncture_dt0.01.npz"), "C0", "-"),
            ("SparseLDL, β = 0 (no damping)", load("s2_puncture_dt0.01_beta0.npz"), "C1", "--")]
    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    for name, r, c, ls in runs:
        kp, c0, a1, a2, r2 = fit_quadratic(r)
        imax = int(np.argmax(r["depth"]))
        print(f"{name}: 刺穿压深 {r['depth'][kp]:.2f} mm，刺穿前峰值 {r['F'][kp-1, 2]:.3f} N；"
              f"刺穿前拟合 F = {c0:.4f} + {a1:.4f}·d + {a2:.5f}·d²（N, mm），r² = {r2:.4f}", flush=True)
        ins, ret = slice(0, imax + 1), slice(imax, None)
        ax[0, 0].plot(r["depth"][ins], r["F"][ins, 2], c=c, ls=ls, label=name + " — insert")
        ax[0, 0].plot(r["depth"][ret], r["F"][ret, 2], c=c, ls=ls, alpha=0.4, lw=3, label=name + " — retract")
        dd = np.linspace(0.01, r["depth"][kp - 1], 50)
        ax[0, 0].plot(dd, c0 + a1 * dd + a2 * dd ** 2, c="k", ls=":", lw=1)
        ax[0, 1].plot(r["t"], r["F"][:, 2], c=c, ls=ls, label=name)
        ax[1, 0].plot(r["t"] - r["t"][kp], r["probe_z"], c=c, ls=ls, marker=".", ms=3, label=name)
        ax[1, 1].plot(r["t"], r["n_cp"], c=c, ls=ls, label=name)
    ax[0, 0].set(xlabel="tip depth below original surface [mm]", ylabel="axial force on needle Fz [N]",
                 title="Force–depth (dotted: fit c0 + a1·d + a2·d²)")
    ax[0, 1].set(xlabel="time [s]", ylabel="Fz [N]", title="Force–time (insert 5 mm/s → hold 1 s → retract)")
    ax[1, 0].set(xlabel="time after puncture [s]", ylabel="z of top-centre node [mm]", xlim=(-0.3, 0.8),
                 title="Surface rebound after puncture")
    ax[1, 1].set(xlabel="time [s]", ylabel="number of coupling points", title="Coupling points along the shaft")
    for a in ax.ravel():
        a.grid(alpha=0.3); a.legend(fontsize=8)
    fig.suptitle("Step 2 puncture: W50, rigid needle, E = 5 kPa, ν = 0.45, threshold 0.22 N, no friction / cutting yet")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "s2_force_depth.png"), dpi=120)
    print("→ results/step2/s2_force_depth.png", flush=True)
