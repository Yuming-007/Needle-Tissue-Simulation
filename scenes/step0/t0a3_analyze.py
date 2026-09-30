"""0a-3 结果汇总：读 results/step0/t0a3_dynamics.json。

自由释放：用 matrix pencil 法从 z(t) 中分离模态（频率、衰减率 σ、振幅），主模态和后向 Euler 的理论值比较。
  对无阻尼振子 x'' = -ω²x，SOFA 的格式（v_{n+1} = v_n + h a_{n+1}，x_{n+1} = x_n + h v_{n+1}）的特征值
  z = 1/(1 - iωh)，所以 σ_th = -ln|z|/h = ln(1+(ωh)²)/(2h)，f_th = arg(z)/(2πh)。
  （t0a3_dynamics.py 打印的"一个周期后振幅比"被第二模态干扰，不作为结论依据。）
加载：分段统计（δ < 1 mm 含起步瞬态；δ = 1–2 mm），以及 Rayleigh 刚度阻尼的预测 ΔF/F ≈ r_s·v/δ。
"""
import sys
sys.dont_write_bytecode = True
import os, json
import numpy as np

R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../results/step0")
d = json.load(open(os.path.join(R, "t0a3_dynamics.json")))


def pencil(t, u, M, tmax):
    m = t <= tmax
    t, u = t[m], u[m]
    h = t[1] - t[0]
    n = len(u); Lp = n // 2
    Y = np.array([u[i:i + Lp + 1] for i in range(n - Lp)])
    _, _, Vt = np.linalg.svd(Y, full_matrices=False)
    V = Vt[:M].T
    z = np.linalg.eigvals(np.linalg.pinv(V[:-1]) @ V[1:])
    Z = np.vander(z, n, increasing=True).T
    a = np.linalg.lstsq(Z, u.astype(complex), rcond=None)[0]
    lam = np.log(z) / h
    modes = [(l.imag / (2 * np.pi), -l.real, 2 * abs(ai)) for l, ai in zip(lam, a) if l.imag > 0.5]
    over = [(0.0, -l.real, abs(ai)) for l, ai in zip(lam, a) if abs(l.imag) <= 0.5 and l.real < 0]  # 过阻尼（实根）
    return sorted(modes, key=lambda x: -x[2]), sorted(over, key=lambda x: -x[2])


out = {}
t, u = np.array(d["A1"]["t"]), np.array(d["A1"]["u"])
modes, _ = pencil(t, u, 10, 0.6)
f1 = modes[0][0]; w1 = 2 * np.pi * f1
print(f"A1（dt = 1e-4）：主模态 f1 = {f1:.3f} Hz（T1 = {1/f1:.4f} s），振幅 {modes[0][2]*1e3:.4f} mm；"
      f"第二模态 {modes[1][0]:.2f} Hz，振幅 {modes[1][2]*1e3:.4f} mm；估计值 4L/c = {d['A1']['estimate_4L_over_c']:.3f} s")
out["f1"] = f1

print("\nA2 数值阻尼（无 Rayleigh）：")
out["A2"] = []
for r in d["A2"]:
    t, u, h = np.array(r["t"]), np.array(r["u"]), r["dt"]
    modes, _ = pencil(t, u, min(10, int(0.6 / h) // 3 * 2), 0.6 if h < 2e-2 else 1.0)
    f, s, a = modes[0]
    z = 1 / (1 - 1j * w1 * h)
    s_th, f_th = -np.log(abs(z)) / h, np.angle(z) / h / (2 * np.pi)
    zeta = s / (2 * np.pi * f)
    decay_T = np.exp(-s / f)
    out["A2"].append(dict(dt=h, f=f, sigma=s, f_th=f_th, sigma_th=s_th, zeta_eff=zeta, decay_per_period=decay_T))
    print(f"  dt = {h:.0e}：f = {f:5.2f} Hz（理论 {f_th:5.2f}），σ = {s:6.2f}/s（理论 {s_th:6.2f}），"
          f"等效阻尼比 ζ = {zeta:.3f}，每周期衰减到 {decay_T:.3f}")

print("\nA3 dt = 0.01 加 Rayleigh：")
for r in d["A3"]:
    t, u = np.array(r["t"]), np.array(r["u"])
    modes, over = pencil(t, u, 10, 0.6)
    osc = f"振荡模态 f = {modes[0][0]:.2f} Hz σ = {modes[0][1]:.2f}/s 振幅 {modes[0][2]*1e3:.4f} mm" if modes else "无振荡模态"
    ov = f"；过阻尼分量 σ = {over[0][1]:.2f}/s（时间常数 {1/over[0][1]:.3f} s）振幅 {over[0][2]*1e3:.4f} mm" if over else ""
    zeta_r = r["rayleighStiffness"] * w1 / 2 + r["rayleighMass"] / (2 * w1)
    print(f"  r_m = {r['rayleighMass']}, r_s = {r['rayleighStiffness']}（对 f1 的物理阻尼比 ζ = r_m/(2ω)+r_sω/2 = {zeta_r:.3f}）：{osc}{ov}")

print("\nB 准静态加载（反力 vs 静力解）：")
for r in d["B"]:
    s1, s2 = r["segments"]["0.25-1mm"], r["segments"]["1-2mm"]
    pred = r["rayleighStiffness"] * r["v"] / 1.5e-3
    extra = f"｜r_s·v/δ(1.5 mm) 预测 {pred:+.4f}" if r["rayleighStiffness"] > 0 else ""
    print(f"  v = {r['v']*1e3:4.0f} mm/s dt = {r['dt']:.0e} r_m = {r['rayleighMass']} r_s = {r['rayleighStiffness']}："
          f"δ<1 mm 最大 {s1['max_abs_rel']:.3f}；δ = 1–2 mm 最大 {s2['max_abs_rel']:.3f} 平均 {s2['mean_rel']:+.4f}{extra}；"
          f"停止后 {r['hold_settle_1pct']:.2f} s 进入 ±1%；每步 {r['t_step']*1e3:.1f} ms")
