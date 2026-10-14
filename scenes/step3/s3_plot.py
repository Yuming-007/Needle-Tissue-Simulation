"""第 3 步的图：3A 侧向力 / 耦合点数随深度；3B 侧向力–时间、侧向力–平移距离、保持结束时各耦合点的力；R3 和 R1 对照。
读 results/step3/s3_{deep,lateral}_{ldl,async}.npz（s3_insertion.py 的结果）。
第二张图：耦合点间距 4 / 5 / 6 mm 时各点的力和合力（s3_lateral_ldl[_tip4|_tip6].npz），以及底面反力和针上的力的对比。
运行：python3 s3_plot.py → results/step3/s3_insertion.png、s3_spacing.png
"""
import sys
sys.dont_write_bytecode = True
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../results/step3")
load = lambda n: dict(np.load(os.path.join(OUT, n)))

if __name__ == "__main__":
    dl, da = load("s3_deep_ldl.npz"), load("s3_deep_async.npz")
    ll, la = load("s3_lateral_ldl.npz"), load("s3_lateral_async.npz")
    fig, ax = plt.subplots(2, 3, figsize=(16, 9))

    a = ax[0, 0]
    imax = int(np.argmax(dl["depth"]))
    ins = slice(0, imax + 1)
    a.plot(dl["depth"][ins], np.hypot(*dl["F"][ins, :2].T) * 1e3, label="lateral |F_xy| (R1)")
    a.plot(dl["depth"][ins], np.abs(dl["F"][ins, 2]) * 1e3, "--", label="|F_z| (R1)")
    a.set(xlabel="tip depth [mm]", ylabel="force on needle [mN]", ylim=(0, 30),
          title="3A insertion to 50 mm: forces vs depth (insert phase)")
    b = a.twinx(); b.plot(dl["depth"][ins], dl["n_cp"][ins], "k:", lw=1); b.set_ylabel("coupling points (dotted)")
    a.legend(loc="upper left", fontsize=8)

    a = ax[0, 1]
    a.plot(dl["t"], dl["step_ms"], lw=0.5, label="R1 SparseLDL")
    a.plot(da["t"], da["step_ms"], lw=0.5, label="R3 Async")
    a.axhline(10, c="k", ls=":", lw=1); a.axhline(33.3, c="r", ls=":", lw=1)
    a.set(xlabel="time [s]", ylabel="animate() time per step [ms]", ylim=(0, 200),
          title="3A step time, no rendering (R1 run in parallel with 3 others, R3 alone)")
    a.legend(fontsize=8)

    a = ax[0, 2]
    for r, name, ls in [(ll, "R1", "-"), (la, "R3", "--")]:
        a.plot(r["t"], r["F"][:, 0] * 1e3, ls=ls, label=f"F_x ({name})")
        a.plot(r["t"], r["F"][:, 1] * 1e3, ls=ls, label=f"F_y ({name})")
    b = a.twinx(); b.plot(ll["t"], ll["dx"], "k:", lw=1); b.set_ylabel("needle offset dx [mm] (dotted)")
    a.set(xlabel="time [s]", ylabel="force on needle [mN]", title="3B lateral translation: force vs time")
    a.legend(fontsize=8)

    a = ax[1, 0]
    m = ll["n_cp"] > 0
    a.plot(ll["dx"][m], ll["F"][m, 0] * 1e3, label="R1")
    a.plot(la["dx"][la["n_cp"] > 0], la["F"][la["n_cp"] > 0, 0] * 1e3, "--", label="R3")
    a.set(xlabel="needle offset dx [mm]", ylabel="F_x on needle [mN]",
          title="3B F_x vs offset (loop width = viscous force, beta = 0.1)")
    a.legend(fontsize=8)

    a = ax[1, 1]
    i_end = int(np.where(np.isclose(ll["dx"], ll["dx"].max()))[0][-1])
    cd, cf = ll["cp_depth"][i_end], ll["cp_f"][i_end]
    k = ~np.isnan(cd)
    a.barh(cd[k], cf[k, 0] * 1e3, height=2.0)
    a.invert_yaxis()
    a.set(xlabel="F_x on needle at each coupling point [mN]", ylabel="coupling point depth [mm]",
          title=f"3B per-point force at end of hold (sum {np.nansum(cf[k, 0])*1e3:.0f} mN)")
    a.axvline(0, c="k", lw=0.5)

    a = ax[1, 2]
    a.plot(ll["t"], ll["probe_dx"], lw=1)
    a.plot(ll["t"], ll["dx"], "k:", lw=1)
    a.legend([f"node {r:.0f} mm from needle, {d:.0f} mm deep" for r, d in zip(ll["probe_r"], ll["probe_depth"])]
             + ["needle offset"], fontsize=7)
    a.set(xlabel="time [s]", ylabel="x displacement [mm]", title="3B tissue displacement near the needle")

    for a in ax.ravel():
        a.grid(alpha=0.3)
    fig.suptitle("Step 3 insertion guidance: W50, rigid needle, E = 5 kPa, nu = 0.45, beta = 0.1, no friction")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "s3_insertion.png"), dpi=110)
    print("→ results/step3/s3_insertion.png", flush=True)

    runs = [("4 mm", load("s3_lateral_ldl_tip4.npz")), ("5 mm", ll), ("6 mm", load("s3_lateral_ldl_tip6.npz"))]
    fig, ax = plt.subplots(1, 3, figsize=(16, 5.5))
    a = ax[0]
    for j, (name, r) in enumerate(runs):
        i_end = int(np.where(np.isclose(r["dx"], r["dx"].max()))[0][-1])
        cd, cf = r["cp_depth"][i_end], r["cp_f"][i_end]
        k = ~np.isnan(cd)
        a.barh(cd[k] + (j - 1) * 0.9, cf[k, 0] * 1e3, height=0.9,
               label=f"spacing {name}: {k.sum()} points, sum {np.nansum(cf[k, 0])*1e3:.0f} mN")
    a.invert_yaxis(); a.axvline(0, c="k", lw=0.5)
    a.set(xlabel="F_x on needle at each coupling point [mN]", ylabel="coupling point depth [mm]",
          title="3B per-point force at end of hold vs coupling spacing")
    a.legend(fontsize=8)
    a = ax[1]
    for name, r in runs:
        m = r["n_cp"] > 0
        a.plot(r["dx"][m], r["F"][m, 0] * 1e3, label=f"spacing {name}")
    a.set(xlabel="needle offset dx [mm]", ylabel="total F_x on needle [mN]", title="3B total force vs offset, by spacing")
    a.legend(fontsize=8)
    a = ax[2]
    ok = np.isfinite(ll["R_bottom"][:, 0])
    a.plot(ll["t"], ll["F"][:, 0] * 1e3, label="F_x on needle (from constraint forces)")
    a.plot(ll["t"][ok], ll["R_bottom"][ok, 0] * 1e3, "--", label="R_x: support force on tissue (bottom nodes)")
    a.plot(ll["t"], ll["F"][:, 2] * 1e3, label="F_z on needle")
    a.plot(ll["t"][ok], ll["R_bottom"][ok, 2] * 1e3, "--", label="R_z (during translation: 2nd-order error of K(x)u)")
    a.set(xlabel="time [s]", ylabel="force [mN]", title="3B equilibrium: needle force vs bottom reaction (5 mm)")
    a.legend(fontsize=8)
    for a in ax.ravel():
        a.grid(alpha=0.3)
    fig.suptitle("Step 3 checks: per-point force depends on spacing, total force and bottom reaction agree")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "s3_spacing.png"), dpi=110)
    print("→ results/step3/s3_spacing.png", flush=True)
