"""第 1.5 步 V1 补充：h = 10 mm 圆盘，动态接触合力与 LCP 参照之差随压深的变化（刚接触后的惯性过渡）。
复用 v1_multipoint_lcp 的 dynamic / node_compliance / lcp_reference，只增加探测压深。
运行：python3 v1_transient.py → results/step1_5/v1_transient.json
"""
import sys
sys.dont_write_bytecode = True
import os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import v1_multipoint_lcp as M

M.PROBE_MM = [0.002, 0.004, 0.006, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05]
n = 9
X0 = M.common.grid_nodes(M.L, n)
pr = M.dynamic(n)
S = sorted({j for p in pr for x in p["rows"] for j in x["nodes"]})
C = M.node_compliance(n, np.array(S))
out = []
for p in pr:
    ref = M.lcp_reference(p["rows"], S, C, X0)
    nv = np.array([x["n"] for x in p["rows"]])
    Fr = float((ref["lam"][:, None] * nv).sum(0)[2])
    out.append(dict(depth_mm=p["depth_mm"], F_dyn=p["F"][2], F_lcp=Fr, rel=p["F"][2] / Fr - 1))
    print(f"压深 {p['depth_mm']:.3f} mm：F_dyn {p['F'][2]:.5e} N，F_LCP {Fr:.5e} N，相对差 {p['F'][2]/Fr-1:+.2e}", flush=True)
json.dump(out, open(os.path.join(M.OUT, "v1_transient.json"), "w"), indent=1)
