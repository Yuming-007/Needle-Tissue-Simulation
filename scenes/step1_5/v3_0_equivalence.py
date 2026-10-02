"""第 1.5 步 V3-0：新 LCP 工具（lcp_tool.py）的等价检查。见 step1_5_results.md V3-0。

(a) 同一参照问题、只换计算方式：重跑 V1 的动态接触（v1_multipoint_lcp.dynamic）拿到压深 0.05 mm 时实际使用的
    投影节点、权重、方向 n 和间隙 q，用 SOFA 的切线刚度 K 一次分解求 W = Bᵀ K⁻¹ B、解 LCP，和 V1 存下的
    参照合力（节点柔度 ±1 mN 对称差分）比较。候选检查量级 0.1%（未冻结）。
(b) 完整新流程（未变形顶面上自己求投影和权重，n = +z，q = −δ）vs V1 参照：差别包含"变形后 vs 未变形的雅可比"。
(c) 四分之一模型 vs 全模型（新流程，点距 s = h）：检验对称边界条件和对称面上点的载荷系数 φ。
h = 10、5 mm，a_eff = 10 mm 圆盘。
运行：python3 v3_0_equivalence.py → results/step1_5/v3_0_equivalence.json
"""
import sys
sys.dont_write_bytecode = True
import os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import lcp_tool as T
import v1_multipoint_lcp as V1

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../results/step1_5")
A = 0.010


def check_a(mesh, n, v1_ref):
    V1.PROBE_MM = [0.05]
    p = V1.dynamic(n)[0]
    rows = p["rows"]
    X = mesh.X
    B = np.zeros((3 * len(X), len(rows))); q = np.zeros(len(rows))
    for r, x in enumerate(rows):
        nv = np.array(x["n"])
        for j, w in zip(x["nodes"], x["w"]):
            B[3 * j:3 * j + 3, r] += w * nv
        q[r] = nv @ (np.array(x["xp"]) - sum(w * X[j] for j, w in zip(x["nodes"], x["w"])))
    Bf = B[mesh.free]
    W = Bf.T @ mesh.lu.solve(Bf); W = 0.5 * (W + W.T)
    mu, g = T.solve_lcp(W, q)
    nv = np.array([x["n"] for x in rows])
    F = (mu[:, None] * nv).sum(0)[2]
    act_new = mu > 1e-9 * mu.max()
    act_v1 = np.array(v1_ref["active_ref"])
    return dict(depth_mm=p["depth_mm"], F_new=float(F), F_v1=v1_ref["F_ref"][2], rel=float(F / v1_ref["F_ref"][2] - 1),
                lam_rel_max=float(np.abs(mu - np.array(v1_ref["lam_ref"])).max() / max(v1_ref["lam_ref"])),
                same_active=bool((act_new == act_v1).all()))


if __name__ == "__main__":
    v1 = json.load(open(os.path.join(OUT, "v1_multipoint_lcp.json")))
    out = []
    for h, n in [(0.010, 9), (0.005, 17)]:
        a = time.perf_counter(); full = T.Mesh(h, quarter=False); t_full = time.perf_counter() - a
        a = time.perf_counter(); qm = T.Mesh(h, quarter=True); t_q = time.perf_counter() - a
        ref = [d for r in v1 if r["n"] == n for d in r["probes"] if abs(d["depth_mm"] - 0.05) < 1e-6][0]
        ra = check_a(full, n, ref)
        k_v1 = ref["F_ref"][2] / (ref["depth_mm"] * 1e-3)
        rf = T.disk_stiffness(full, A, h)
        rq = T.disk_stiffness(qm, A, h)
        d = dict(h_mm=h * 1e3, a=ra, k_v1=k_v1, full=rf, quarter=rq, rel_b=rf["k"] / k_v1 - 1, rel_c=rq["k"] / rf["k"] - 1,
                 t_mesh_full=t_full, t_mesh_quarter=t_q, ndof_full=full.ndof, ndof_quarter=qm.ndof)
        out.append(d)
        print(f"== h = {h*1e3:.0f} mm（自由度 全 {full.ndof} / 四分之一 {qm.ndof}；建网格 + 分解 {t_full:.1f} / {t_q:.1f} s）", flush=True)
        print(f"  (a) 同一雅可比：F 新 {ra['F_new']:.6e} N，V1 {ra['F_v1']:.6e} N，相对差 {ra['rel']:+.2e}；逐点 λ 最大差 / λmax {ra['lam_rel_max']:.1e}；"
              f"受力点集合一致 {ra['same_active']}", flush=True)
        print(f"  (b) 新流程（未变形雅可比）全模型 k = {rf['k']:.4f} N/m，V1 k = {k_v1:.4f} N/m，相对差 {d['rel_b']:+.2e}；"
              f"受力点 {rf['n_active']}/{rf['n_points']}，互补残差 {rf['complementarity']:.1e}", flush=True)
        print(f"  (c) 四分之一模型 k = {rq['k']:.4f} N/m，相对全模型 {d['rel_c']:+.2e}；受力点 {rq['n_active']}/{rq['n_points']}", flush=True)
        json.dump(out, open(os.path.join(OUT, "v3_0_equivalence.json"), "w"), indent=1)
