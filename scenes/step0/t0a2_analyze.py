"""0a-2 结果汇总：读 results/step0/t0a2_punch_{coarse,quarter}.json，输出比值、有限尺寸效应、收敛阶。"""
import sys
sys.dont_write_bytecode = True
import os, json
import numpy as np

R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../results/step0")
C = json.load(open(os.path.join(R, "t0a2_punch_coarse.json")))
Q = json.load(open(os.path.join(R, "t0a2_punch_quarter.json")))
E, a = 5000.0, 0.01
kS = lambda nu: 2 * a * E / (1 - nu ** 2)


def get(D, **kw):
    r = [d for d in D if all(np.isclose(d.get(k), v) if isinstance(v, float) else d.get(k) == v for k, v in kw.items())
         and d.get("delta_over_a", 0.01) == 0.01]
    assert len(r) == 1, (kw, len(r))
    return r[0]


print("== 全模型，8 cm，实时网格 h = 10 mm 和 h = 5 mm ==")
for nu in [0.3, 0.45, 0.49]:
    for h in [10.0, 5.0]:
        d = get(C, L=0.08, h_mm=h, nu=nu, bonded=False)
        print(f"nu={nu} h={h:4.1f}: k={d['k']:6.1f} kS={kS(nu):6.1f}  k/kS(a)={d['ratio_nominal']:.3f}  k/kS(a_eff)={d['ratio_eff']:.3f}")

out = {}
print("\n== 四分之一模型：有限尺寸效应（h = 10 mm），k(L)/k(32 cm) 和 1/L 外推的 k∞ ==")
for nu in [0.3, 0.45, 0.49]:
    ks = {L: get(Q, L=L, h_mm=10.0, nu=nu)["k"] for L in [0.08, 0.16, 0.24, 0.32]}
    c = (ks[0.16] - ks[0.32]) / (1 / 0.16 - 1 / 0.32); kinf = ks[0.32] - c / 0.32
    out[("size", nu)] = ks[0.08] / kinf
    print(f"nu={nu}: " + "  ".join(f"L={L*100:.0f}cm {ks[L]/ks[0.32]:.4f}" for L in ks) +
          f"  | k(8cm)/k∞ ≈ {ks[0.08]/kinf:.3f}")
    k5 = {L: get(Q, L=L, h_mm=5.0, nu=nu)["k"] for L in [0.08, 0.16]}
    print(f"        h=5mm: k(8cm)/k(16cm) = {k5[0.08]/k5[0.16]:.4f}   (h=10mm: {ks[0.08]/ks[0.16]:.4f})")

print("\n== 四分之一模型：网格收敛（8 cm，h = 10/5/2.5 mm），观测收敛阶 p 和 Richardson 外推 ==")
for nu in [0.3, 0.45, 0.49]:
    k = [get(Q, L=0.08, h_mm=h, nu=nu)["k"] for h in [10.0, 5.0, 2.5]]
    p = np.log2((k[0] - k[1]) / (k[1] - k[2]))
    kx = k[2] - (k[1] - k[2]) / (2 ** p - 1)
    out[("conv", nu)] = (p, kx)
    print(f"nu={nu}: k = {k[0]:.1f}, {k[1]:.1f}, {k[2]:.1f} N/m  p ≈ {p:.2f}  k(h→0) ≈ {kx:.1f}  "
          f"k(h→0)/kS = {kx/kS(nu):.3f}  | 再除以尺寸因子 → 半空间估计 {kx/out[('size', nu)]/kS(nu):.3f}")

print("\n== 四分之一 vs 全模型（同网格，四面体切分方向的影响）==")
for nu in [0.3, 0.45, 0.49]:
    for h in [10.0, 5.0]:
        print(f"nu={nu} h={h}: k_Q/k_F = {get(Q, L=0.08, h_mm=h, nu=nu)['k']/get(C, L=0.08, h_mm=h, nu=nu, bonded=False)['k']:.4f}")

print("\n== 其他（全模型 8 cm h = 5 mm）==")
for nu in [0.3, 0.49]:
    b, f = get(C, L=0.08, h_mm=5.0, nu=nu, bonded=True), get(C, L=0.08, h_mm=5.0, nu=nu, bonded=False)
    print(f"nu={nu}: bonded/frictionless = {b['k']/f['k']:.4f}")
lin = [d for d in C if "delta_over_a" in d] + [get(C, L=0.08, h_mm=5.0, nu=0.3, bonded=False)]
for d in sorted(lin, key=lambda d: d.get("delta_over_a", 0.01)):
    print(f"delta/a = {d.get('delta_over_a', 0.01):.3f}: k = {d['k']:.2f} N/m")
