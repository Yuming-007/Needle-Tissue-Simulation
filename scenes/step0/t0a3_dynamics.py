"""0a-3 动力学基线。见 docs/plan/step0_foundation.md §4。

实时网格：8 cm 立方体，h = 10 mm（n = 9），ν = 0.45（验证用参数），E = 5 kPa，ρ = 1000 kg/m³，重力关闭。
底面全固定，侧面和顶面自由。求解器：EulerImplicitSolver + SparseLDLSolver（common.add_dynamic_solver）。

A 自由释放：初始状态 = 弱体力（竖直向下 0.3 m/s² 的"重力"）下的静力解，速度为 0；然后撤掉体力，自由振动。
   体力下的静变形形状接近竖向最低阶模态（Rayleigh–Ritz），比局部压头变形激发的高阶模态少得多
   （第一版用压头变形作初始状态，多模态叠加使"一个周期后振幅比"无法和单模态理论比较）。
   观测顶面中心节点的 z(t)。
   A1：dt = 1e-4、无阻尼 → 最低阶振动的周期 T1（FFT 峰值 + 过零点）。
   A2：dt = 1e-3 … 5e-2、无阻尼 → 数值阻尼：第一个周期后的振幅比、衰减到 1% 的时间；
       和后向 Euler 的理论值比较：对无阻尼振子每步放大因子 |λ| = 1/sqrt(1+(ω dt)²)。
   A3：dt = 0.01，加 Rayleigh 阻尼对照。
B 准静态加载：压头节点（r ≤ a，只约束 z）以速度 v 下压到 δ_max，然后保持。
   驱动方式：PartialFixedProjectiveConstraint [0,0,1] 只清零 z 方向的增量；把压头节点的 z 速度设为 -v，
   隐式方程右端的 h·K·v 项使压头运动一致地进入每一步的线性化；到达 δ_max 的那一步起速度设为 0。
   反力：第 k+1 步之后读 force，得到 x_{k+1} 上的内力（见 common.add_dynamic_solver 的说明）。
   和同一网格、同一压深的静力解 F_static(δ) 比较；分两段统计：δ < 1 mm（含起步瞬态）和 δ = 1–2 mm。
   B2：v = 5 mm/s、dt = 0.01，加 Rayleigh 阻尼对照。force 里只有弹性内力；加阻尼时压头受到的总反力还要加上
       刚度阻尼力 -r_s (K v)（用 common.ForceProbe 有限差分计算）。质量阻尼 -r_m M v 的压头分量很小（约 1e-5 N），不计入。
运行：python3 t0a3_dynamics.py → results/step0/t0a3_dynamics.json
"""
import sys
sys.dont_write_bytecode = True
import os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
import Sofa, Sofa.Core, Sofa.Simulation
import common

L, N, E, NU, RHO, A_PUNCH = 0.08, 9, 5000.0, 0.45, 1000.0, 0.010
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../results/step0/t0a3_dynamics.json")

X = common.grid_nodes(L, N)
BOT = np.where(np.isclose(X[:, 2], 0))[0]
TOP = np.isclose(X[:, 2], L)
R = np.hypot(X[:, 0] - L / 2, X[:, 1] - L / 2)
PUNCH = np.where(TOP & (R <= A_PUNCH + 1e-9))[0]
CENTER = int(np.where(TOP & (R < 1e-9))[0][0])


def static_punch(delta):
    """静力解：压头压深 delta。返回 (节点位置, 压头反力)。"""
    root = common.make_root()
    common.add_static_solver(root, linear="ldl", abs_tol=1e-9)
    t, mo, _ = common.add_tissue(root, L=L, n=N, E=E, nu=NU, rho=RHO)
    cons = [t.addObject("FixedProjectiveConstraint", indices=BOT.tolist()),
            t.addObject("PartialFixedProjectiveConstraint", indices=PUNCH.tolist(), fixedDirections=[0, 0, 1])]
    Sofa.Simulation.init(root)
    p = X.copy(); p[PUNCH, 2] -= delta
    mo.position.value = p
    Sofa.Simulation.animate(root, root.dt.value)
    x = mo.position.array().copy()
    f = common.internal_force(root, t, cons)
    return x, float(f[PUNCH, 2].sum())


def static_bodyforce(g):
    """静力解：竖直体力 g（m/s²，负值向下），底面固定。返回节点位置。"""
    root = common.make_root(gravity=(0.0, 0.0, g))
    common.add_static_solver(root, linear="ldl", abs_tol=1e-10)
    t, mo, _ = common.add_tissue(root, L=L, n=N, E=E, nu=NU, rho=RHO)
    t.addObject("FixedProjectiveConstraint", indices=BOT.tolist())
    Sofa.Simulation.init(root)
    Sofa.Simulation.animate(root, root.dt.value)
    return mo.position.array().copy()


def dyn_scene(dt, rm=0.0, rs=0.0, punch=False):
    root = common.make_root(dt=dt)
    common.add_dynamic_solver(root, rayleigh_mass=rm, rayleigh_stiffness=rs)
    t, mo, _ = common.add_tissue(root, L=L, n=N, E=E, nu=NU, rho=RHO)
    t.addObject("FixedProjectiveConstraint", indices=BOT.tolist())
    if punch:
        t.addObject("PartialFixedProjectiveConstraint", indices=PUNCH.tolist(), fixedDirections=[0, 0, 1])
    Sofa.Simulation.init(root)
    return root, mo


def free_release(x0, z_eq, dt, t_end, rm=0.0, rs=0.0):
    root, mo = dyn_scene(dt, rm, rs)
    mo.position.value = x0
    nsteps = int(round(t_end / dt))
    z = np.empty(nsteps + 1); z[0] = x0[CENTER, 2]
    tt = []
    for k in range(nsteps):
        t0 = time.perf_counter()
        Sofa.Simulation.animate(root, dt)
        tt.append(time.perf_counter() - t0)
        z[k + 1] = mo.position.array()[CENTER, 2]
    return np.arange(nsteps + 1) * dt, z - z_eq, float(np.median(tt))


def period_fft(t, u):
    u = u - u.mean()
    n = 1 << (len(u) * 4 - 1).bit_length()
    spec = np.abs(np.fft.rfft(u * np.hanning(len(u)), n))
    freq = np.fft.rfftfreq(n, t[1] - t[0])
    i = np.argmax(spec[1:]) + 1
    return 1.0 / freq[i]


def period_zero_cross(t, u):
    s = np.where(np.diff(np.sign(u)) != 0)[0]
    tc = t[s] - u[s] * (t[s + 1] - t[s]) / (u[s + 1] - u[s])   # 线性插值的过零时刻
    return 2 * np.mean(np.diff(tc[:6])) if len(tc) >= 3 else float("nan"), tc


def quasi_static(v, dt, delta_max, t_hold, delta_tab, F_tab, rm=0.0, rs=0.0):
    root, mo = dyn_scene(dt, rm, rs, punch=True)
    n_load = int(round(delta_max / v / dt))
    n_hold = int(round(t_hold / dt))
    vel = np.zeros_like(X); vel[PUNCH, 2] = -v
    mo.velocity.value = vel
    probe = common.ForceProbe(L=L, n=N, E=E, nu=NU, rho=RHO) if rs > 0 else None
    xs, Fs, Fd, ts, tt = [], [], [], [], []
    for k in range(n_load + n_hold):
        if k == n_load:                        # 到达 δ_max：停止
            vv = mo.velocity.array().copy(); vv[PUNCH, 2] = 0.0
            mo.velocity.value = vv
        t0 = time.perf_counter()
        Sofa.Simulation.animate(root, dt)
        tt.append(time.perf_counter() - t0)
        if k > 0:                              # 这一步读到的 force 对应上一步结束时的位置
            Fs.append(float(mo.force.array()[PUNCH, 2].sum()))
        xs.append(mo.position.array().copy())
        if probe is not None:                  # 刚度阻尼力（向上为正）：-r_s (K v) 的 z 分量之和
            Fd.append(float(-rs * probe.K_times(xs[-1], mo.velocity.array().copy())[PUNCH, 2].sum()))
        ts.append((k + 1) * dt)
    xs, ts = xs[:-1], ts[:-1]                  # 最后一个位置没有对应的力
    F_damp = np.array(Fd[:-1]) if probe is not None else np.zeros(len(Fs))
    delta = np.array([X[PUNCH[0], 2] - x[PUNCH[0], 2] for x in xs])
    Fs = np.array(Fs) + F_damp                 # 总反力 = 弹性 + 刚度阻尼
    Fst = np.interp(delta, delta_tab, F_tab)
    return dict(t=np.array(ts), delta=delta, F=Fs, F_damp=F_damp, F_static=Fst, n_load=n_load, t_step=float(np.median(tt)))


if __name__ == "__main__":
    out = {"params": dict(L=L, n=N, h_mm=1000 * L / (N - 1), E=E, nu=NU, rho=RHO, a=A_PUNCH)}

    # ---------- A 自由释放 ----------
    g0 = -0.3
    x0 = static_bodyforce(g0)
    z0 = x0[CENTER, 2] - X[CENTER, 2]
    print(f"初始状态：体力 {g0} m/s² 下的静变形，顶面中心初始位移 {z0*1e3:.4f} mm", flush=True)

    t, u, tstep = free_release(x0, X[CENTER, 2], 1e-4, 0.6)
    T_fft = period_fft(t, u)
    T_zc, _ = period_zero_cross(t, u)
    T1 = T_fft
    w1 = 2 * np.pi / T1
    c = np.sqrt(E / RHO)
    print(f"A1 dt=1e-4 无阻尼：T1(FFT) = {T_fft:.4f} s，T1(过零点) = {T_zc:.4f} s；估计值 4L/c = {4*L/c:.3f} s", flush=True)
    out["A1"] = dict(T1_fft=T_fft, T1_zero_cross=T_zc, estimate_4L_over_c=4 * L / c,
                     t=t[::10].tolist(), u=u[::10].tolist())

    out["A2"] = []
    for dt in [1e-3, 5e-3, 1e-2, 2e-2, 5e-2]:
        t, u, tstep = free_release(x0, X[CENTER, 2], dt, 2.0)
        env = np.abs(u)
        # 第一个周期后的振幅比：取 [T1/2, 3T1/2] 内的最大 |u| / 初始 |u|
        m = (t >= 0.5 * T1) & (t <= 1.5 * T1)
        ratio = env[m].max() / abs(u[0]) if m.any() else float("nan")
        below = np.where(env < 0.01 * abs(u[0]))[0]
        # 1% 以下并且之后不再超过
        t_settle = float("nan")
        for i in below:
            if env[i:].max() < 0.01 * abs(u[0]):
                t_settle = t[i]; break
        theory = (1 + (w1 * dt) ** 2) ** (-(T1 / dt) / 2)
        T_zc, _ = period_zero_cross(t, u)
        r = dict(dt=dt, amp_ratio_after_T1=ratio, theory_single_mode=theory, t_settle_1pct=t_settle,
                 period_measured=T_zc, t_step=tstep, t=t.tolist() if dt >= 5e-3 else t[::5].tolist(),
                 u=u.tolist() if dt >= 5e-3 else u[::5].tolist())
        out["A2"].append(r)
        print(f"A2 dt={dt:.0e}: 一个周期后振幅比 = {ratio:.3f}（单模态理论 {theory:.3f}），衰减到 1% 用时 {t_settle:.3f} s，"
              f"测得周期 {T_zc:.4f} s，每步 {tstep*1e3:.2f} ms", flush=True)

    out["A3"] = []
    for rm, rs in [(0.1, 0.0), (0.0, 0.01), (0.0, 0.1), (0.1, 0.1)]:
        t, u, tstep = free_release(x0, X[CENTER, 2], 1e-2, 2.0, rm, rs)
        env = np.abs(u)
        m = (t >= 0.5 * T1) & (t <= 1.5 * T1)
        ratio = env[m].max() / abs(u[0])
        below = [i for i in np.where(env < 0.01 * abs(u[0]))[0] if env[i:].max() < 0.01 * abs(u[0])]
        t_settle = t[below[0]] if below else float("nan")
        out["A3"].append(dict(rayleighMass=rm, rayleighStiffness=rs, amp_ratio_after_T1=ratio, t_settle_1pct=t_settle,
                              t=t.tolist(), u=u.tolist()))
        print(f"A3 dt=1e-2 Rayleigh m={rm} s={rs}: 一个周期后振幅比 = {ratio:.3f}，衰减到 1% 用时 {t_settle:.3f} s", flush=True)

    # ---------- B 准静态加载 ----------
    dmax = 2e-3
    delta_tab = np.linspace(0, dmax, 9)
    F_tab = np.array([0.0] + [static_punch(d)[1] for d in delta_tab[1:]])
    print("静力曲线 F(δ)：" + "  ".join(f"{d*1e3:.2f}mm:{F:.4f}N" for d, F in zip(delta_tab, F_tab)), flush=True)
    out["static_curve"] = dict(delta=delta_tab.tolist(), F=F_tab.tolist())
    out["B"] = []
    cases = [(v, dt, 0.0, 0.0) for v in [1e-3, 5e-3, 10e-3] for dt in [1e-2, 5e-3]]
    cases += [(5e-3, 1e-2, rm, rs) for rm, rs in [(0.1, 0.0), (0.0, 0.01), (0.0, 0.1)]]
    for v, dt, rm, rs in cases:
            r = quasi_static(v, dt, dmax, 1.0, delta_tab, F_tab, rm, rs)
            nl = r["n_load"] - 1
            rel = r["F"] / np.maximum(r["F_static"], 1e-12) - 1
            # 加载阶段：δ ≥ 0.25 mm 之后的相对偏差（开头的力很小，相对值没有意义）
            mask = (np.arange(len(rel)) < nl) & (r["delta"] >= 0.25e-3)
            abs_dev = r["F"] - r["F_static"]
            F_end_load = r["F"][nl - 1]
            hold = r["F"][nl:]
            Fs_final = F_tab[-1]
            settle = [i for i in range(len(hold)) if np.all(np.abs(hold[i:] / Fs_final - 1) < 0.01)]
            seg = {}
            for name, (lo, hi) in {"0.25-1mm": (0.25e-3, 1e-3), "1-2mm": (1e-3, dmax + 1e-12)}.items():
                mm = (np.arange(len(rel)) < nl) & (r["delta"] >= lo) & (r["delta"] < hi)
                seg[name] = dict(max_abs_rel=float(np.abs(rel[mm]).max()), mean_rel=float(rel[mm].mean()))
            d = dict(v=v, dt=dt, rayleighMass=rm, rayleighStiffness=rs, segments=seg,
                     max_rel_dev_load=float(np.abs(rel[mask]).max()),
                     mean_rel_dev_load=float(rel[mask].mean()),
                     max_abs_dev_load=float(np.abs(abs_dev[:nl]).max()),
                     F_end_load=float(F_end_load), F_static_end=float(Fs_final),
                     hold_settle_1pct=float(settle[0] * dt) if settle else float("nan"),
                     F_hold_final=float(hold[-1]), t_step=r["t_step"],
                     t=r["t"].tolist(), delta=r["delta"].tolist(), F=r["F"].tolist(), F_damp=r["F_damp"].tolist(),
                     F_static=r["F_static"].tolist())
            out["B"].append(d)
            print(f"B v={v*1e3:.0f}mm/s dt={dt:.0e} Rayleigh m={rm} s={rs}: δ<1mm 最大 {seg['0.25-1mm']['max_abs_rel']:.3f}；"
                  f"δ=1-2mm 最大 {seg['1-2mm']['max_abs_rel']:.3f}、平均 {seg['1-2mm']['mean_rel']:+.4f}；"
                  f"停止时 F/F_static = {F_end_load/Fs_final:.4f}；保持后 {d['hold_settle_1pct']:.2f} s 内进入 ±1%，"
                  f"最终 {d['F_hold_final']/Fs_final:.5f}；每步 {r['t_step']*1e3:.2f} ms", flush=True)
    json.dump(out, open(RES, "w"))
