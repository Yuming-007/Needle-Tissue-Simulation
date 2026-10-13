"""第 2 步：组织阻尼的两项对照检查（为什么加 β = 0.1、加了以后实时求解器是否算得对、α 要不要加）。见 docs/plan/step2_results.md §3。

背景：不加阻尼时，刺穿后被压缩约 6 mm 的组织在 2–3 步内弹回；实时求解器 AsyncSparseLDLSolver 用的是后台线程算好的
  旧分解（W50 上分解一次约 85 ms ≈ 8–9 步之前的矩阵），刚度矩阵每步变化太大，几步内发散。加 Rayleigh 刚度阻尼
  β = 0.1 s（Bui 2019 的文献值）后回弹变慢（约 0.1 s 回弹一半），旧分解够用。
算例（轨迹、场景和 s2_puncture.py 完全相同，只换求解器 / 阻尼）：
  ① async_beta0   ：实时求解器，无阻尼 —— 预期发散（记录发散的步数）；
  ② async_beta0.1 ：实时求解器，β = 0.1 —— 和离线主结果 s2_puncture_dt0.01.npz（SparseLDL，β = 0.1）逐步比较；
  ③ ldl_alpha0.01 ：离线求解器，β = 0.1 + 质量阻尼 α = 0.01 1/s（Bui 2019 的 α）—— 和主结果比较，看 α 的影响。
异步求解器的结果和线程时序有关，每次运行略有不同（第 0 步的结论）。
运行：python3 s2_damping_checks.py → results/step2/s2_puncture_check_*.npz + 终端摘要（约 3 分钟，三个算例依次运行）
"""
import sys
sys.dont_write_bytecode = True
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import s2_puncture as S


def puncture_index(rec):
    return int(np.argmax(rec["n_cp"] > 0))


def describe(name, rec):
    kp = puncture_index(rec)
    pz = rec["probe_z"]
    after = lambda n: f"{pz[kp + n]:+.2f}" if kp + n < len(pz) else "—"
    print(f"   {name}：{len(rec['t'])} 步，刺穿压深 {rec['depth'][kp]:.2f} mm；顶面节点 z 刺穿前 {pz[kp-1]:+.2f}，"
          f"刺穿后 3 / 10 / 50 步 {after(3)} / {after(10)} / {after(50)} mm", flush=True)


def compare(name, rec, ref):
    n = min(len(rec["t"]), len(ref["t"]))
    kp = puncture_index(ref)
    dF = np.abs(rec["F"][:n, 2] - ref["F"][:n, 2])
    dz = np.abs(rec["probe_z"][:n] - ref["probe_z"][:n])
    rel = dF[:kp - 2] / np.maximum(np.abs(ref["F"][:kp - 2, 2]), 0.02)       # 刺穿前（力 < 0.02 N 时按 0.02 N 算相对值）
    same_cp = np.array_equal(rec["n_cp"][:n], ref["n_cp"][:n])
    print(f"   {name} 和主结果比：轴向力最大差 {dF.max():.1e} N（刺穿前最大相对差 {100*rel.max():.1f}%），"
          f"顶面节点 z 最大差 {dz.max():.3f} mm，耦合点数逐步相同：{same_cp}", flush=True)


if __name__ == "__main__":
    ref = dict(np.load(os.path.join(S.OUT, "s2_puncture_dt0.01.npz")))
    assert float(ref["beta"]) == S.BETA and str(ref["linear"]) == "ldl", "先运行 python3 s2_puncture.py 0.01"
    print(f"== 阻尼对照检查（dt = 0.01，主结果 = SparseLDL + β = {S.BETA}）", flush=True)
    describe("主结果 ldl β=0.1", ref)
    r1 = S.run(0.01, beta=0.0, linear="async", tag="check_async_beta0")
    describe("① async β=0", r1)
    r2 = S.run(0.01, beta=S.BETA, linear="async", tag="check_async_beta0.1")
    describe("② async β=0.1", r2)
    compare("② async β=0.1", r2, ref)
    r3 = S.run(0.01, beta=S.BETA, linear="ldl", tag="check_ldl_alpha0.01", alpha=0.01)
    describe("③ ldl β=0.1 α=0.01", r3)
    compare("③ ldl β=0.1 α=0.01", r3, ref)
