"""第 1.5 步 M0：参数化生成非均匀四面体工作网格（gmsh）。见 docs/plan/step1_5_mesh_design.md。

组织块：边长 L 的立方体 [0, L]³（米），底面 z = 0 固定（在 SOFA 侧设置），顶面 z = L 是接触面。
尺寸场（gmsh 官方教程 t10 的组合）：
  - MathEval：到竖直线段（从顶面中心 (L/2, L/2, L) 往下深 D）的距离
      d = sqrt((x − cx)² + (y − cy)² + dz²)，dz = max(0, (L − D) − z) = ((L−D−z) + |L−D−z|)/2；
    用公式写，不需要在体内嵌入几何实体（设计方案 §8）；
  - Threshold：d ≤ R_in 时尺寸 h_local，d ≥ R_in + W 时 h_far，中间线性；W = (h_far − h_local) / g；
  - 作为背景网格；关闭 MeshSizeExtendFromBoundary / MeshSizeFromPoints / MeshSizeFromCurvature（t10 的做法）。
3D 算法 Delaunay（Mesh.Algorithm3D = 1），生成后 Mesh.Optimize + Mesh.OptimizeNetgen（M0 比较：Netgen 优化把 B 的二面角从 13.5–154° 改善到 22.6–137.6°，HXT 算法上 Netgen 优化无效）。
输出（meshes/step1_5/<名字>.*）：
  - .msh：msh 2.2 ASCII，只含物理体（四面体），供 SOFA MeshGmshLoader 使用；
  - .npz：节点 X（米，0 起连续编号，只保留四面体用到的节点）、四面体 tets（有符号体积为正的顺序）、参数；
  - _quality.json：质量报告（见 quality()）。
运行：PYTHONPATH=~/sofa/tools/gmsh-4.15.2-Linux64-sdk/lib python3 make_graded_mesh.py [候选名 ...]
"""
import sys
sys.dont_write_bytecode = True
import os, json, time
import numpy as np

SDK = os.path.expanduser("~/sofa/tools/gmsh-4.15.2-Linux64-sdk/lib")
if SDK not in sys.path:
    sys.path.insert(0, SDK)
import gmsh

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../meshes/step1_5")
L = 0.08

# 第一轮候选（a_eff = 10 mm 固定，设计方案 §4、ChatGPT 复核 2026-10-03）。单位：米。
CANDIDATES = {
    "A-":      dict(a=0.010, h_local=0.005, h_far=0.012, R_in=0.010, D=0.010, g=0.35),
    "A":       dict(a=0.010, h_local=0.005, h_far=0.012, R_in=0.015, D=0.015, g=0.35),
    "A+":      dict(a=0.010, h_local=0.005, h_far=0.012, R_in=0.020, D=0.020, g=0.35),
    "A_far16": dict(a=0.010, h_local=0.005, h_far=0.016, R_in=0.015, D=0.015, g=0.35),
    "A_g05":   dict(a=0.010, h_local=0.005, h_far=0.012, R_in=0.015, D=0.015, g=0.50),
    "B":       dict(a=0.010, h_local=0.004, h_far=0.012, R_in=0.015, D=0.015, g=0.35),
    "E":       dict(a=0.010, h_local=0.0025, h_far=0.012, R_in=0.015, D=0.015, g=0.35),   # 离线参照
    # M2 之后补充（2026-10-03）：加密区更大，但远场更粗、过渡更快，寻找精度和实时的更好折中
    "A_f16g05":   dict(a=0.010, h_local=0.005, h_far=0.016, R_in=0.015, D=0.015, g=0.50),
    "A+_f16g05":  dict(a=0.010, h_local=0.005, h_far=0.016, R_in=0.020, D=0.020, g=0.50),
    "A+_f20g05":  dict(a=0.010, h_local=0.005, h_far=0.020, R_in=0.020, D=0.020, g=0.50),
    "A+_f20g07":  dict(a=0.010, h_local=0.005, h_far=0.020, R_in=0.020, D=0.020, g=0.70),
    # 离线参照：处处 5 mm 的均匀 gmsh 网格（h_far = h_local），用来衡量"非均匀"本身带来的偏差
    "U5":         dict(a=0.010, h_local=0.005, h_far=0.005, R_in=0.015, D=0.015, g=0.35),
    # M3 筛选（2026-10-03，用户 + ChatGPT 讨论）：A+ 和 B 系列；B+ 的加密区和 A+ 相同，A++ / B++ 检查加密区是否够大
    "A++":        dict(a=0.010, h_local=0.005, h_far=0.012, R_in=0.025, D=0.025, g=0.35),
    "B+":         dict(a=0.010, h_local=0.004, h_far=0.012, R_in=0.020, D=0.020, g=0.35),
    "B++":        dict(a=0.010, h_local=0.004, h_far=0.012, R_in=0.025, D=0.025, g=0.35),
    "U4":         dict(a=0.010, h_local=0.004, h_far=0.004, R_in=0.015, D=0.015, g=0.35),   # 离线参照
    # 第 2–6 步的工作网格（用户 2026-10-09）：A+ 的加密区沿竖直针道延长到最大插入深度 50 mm，一次覆盖后续各步
    "W50":        dict(a=0.010, h_local=0.005, h_far=0.012, R_in=0.020, D=0.050, g=0.35),
    "W50c":       dict(a=0.010, h_local=0.005, h_far=0.012, R_in=0.020, D=0.050, g=0.35, center_node=True),   # W50 + 顶面中心节点
}


def target_size(X, p):
    cx = cy = L / 2
    dz = np.maximum(0.0, (L - p["D"]) - X[:, 2])
    d = np.sqrt((X[:, 0] - cx) ** 2 + (X[:, 1] - cy) ** 2 + dz ** 2)
    W = (p["h_far"] - p["h_local"]) / p["g"]
    return np.clip(p["h_local"] + (d - p["R_in"]) * p["g"], p["h_local"], p["h_far"]), d


def generate(name, p, optimize_netgen=True, algo3d=1):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(name)
    gmsh.model.occ.addBox(0, 0, 0, L, L, L)
    gmsh.model.occ.synchronize()
    if p.get("center_node"):
        # 在顶面中心嵌入一个节点（针的入口点），让单点针尖正好压在节点上（2026-10-09）
        top = [t for d, t in gmsh.model.getEntities(2)
               if abs(gmsh.model.occ.getCenterOfMass(2, t)[2] - L) < 1e-9]
        pc = gmsh.model.occ.addPoint(L / 2, L / 2, L, p["h_local"])   # 和立方体同一个内核（OCC），否则嵌入不生效
        gmsh.model.occ.synchronize()
        gmsh.model.mesh.embed(0, [pc], 2, top[0])
    vols = [t for _, t in gmsh.model.getEntities(3)]
    gmsh.model.addPhysicalGroup(3, vols, 1, name="tissue")
    W = max((p["h_far"] - p["h_local"]) / p["g"], 1e-6)
    cx = cy = L / 2
    zc = L - p["D"]
    expr = (f"Sqrt((x-{cx})^2+(y-{cy})^2+((({zc}-z)+Abs({zc}-z))/2)^2)")
    f1 = gmsh.model.mesh.field.add("MathEval")
    gmsh.model.mesh.field.setString(f1, "F", expr)
    f2 = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(f2, "InField", f1)
    gmsh.model.mesh.field.setNumber(f2, "SizeMin", p["h_local"])
    gmsh.model.mesh.field.setNumber(f2, "SizeMax", p["h_far"])
    gmsh.model.mesh.field.setNumber(f2, "DistMin", p["R_in"])
    gmsh.model.mesh.field.setNumber(f2, "DistMax", p["R_in"] + W)
    gmsh.model.mesh.field.setAsBackgroundMesh(f2)
    for opt in ["Mesh.MeshSizeExtendFromBoundary", "Mesh.MeshSizeFromPoints", "Mesh.MeshSizeFromCurvature"]:
        gmsh.option.setNumber(opt, 0)
    gmsh.option.setNumber("Mesh.Algorithm3D", algo3d)
    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 1 if optimize_netgen else 0)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.option.setNumber("Mesh.SaveAll", 0)
    t0 = time.perf_counter()
    gmsh.model.mesh.generate(3)
    t_gen = time.perf_counter() - t0
    # 数组
    ntags, coords, _ = gmsh.model.mesh.getNodes()
    coords = coords.reshape(-1, 3)
    etags, enodes = gmsh.model.mesh.getElementsByType(4)
    enodes = enodes.reshape(-1, 4)
    used = np.unique(enodes)
    pos = {int(t): i for i, t in enumerate(ntags)}
    X_all = coords[[pos[int(t)] for t in used]]
    remap = {int(t): i for i, t in enumerate(used)}
    tets = np.vectorize(lambda t: remap[int(t)])(enodes)
    gq = {q: np.asarray(gmsh.model.mesh.getElementQualities(etags.tolist(), q)) for q in ["minSICN", "gamma"]}
    os.makedirs(OUT, exist_ok=True)
    gmsh.write(os.path.join(OUT, f"{name}.msh"))
    gmsh.finalize()
    # 有符号体积：保证为正（记录翻转了多少个）
    vol = signed_volumes(X_all, tets)
    n_neg = int((vol <= 0).sum())
    tets[vol < 0] = tets[vol < 0][:, [0, 2, 1, 3]]
    return X_all, tets, gq, t_gen, n_neg


def signed_volumes(X, tets):
    a, b, c, d = (X[tets[:, i]] for i in range(4))
    return np.einsum("ij,ij->i", np.cross(b - a, c - a), d - a) / 6.0


def dihedral_angles(X, tets):
    """每个四面体 6 个二面角（度）。"""
    P = X[tets]
    faces = [(1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2)]       # 第 k 个面不含顶点 k
    nrm = []
    for k, (i, j, l) in enumerate(faces):
        n = np.cross(P[:, j] - P[:, i], P[:, l] - P[:, i])
        # 外法向：指向远离对顶点 k
        s = np.sign(np.einsum("ij,ij->i", n, P[:, i] - P[:, k]))[:, None]
        nrm.append(n * s / np.linalg.norm(n, axis=1)[:, None])
    ang = []
    for a in range(4):
        for b in range(a + 1, 4):
            c = -np.einsum("ij,ij->i", nrm[a], nrm[b])
            ang.append(np.degrees(np.arccos(np.clip(c, -1, 1))))
    return np.stack(ang, 1)


def quality(name, X, tets, p, gq, t_gen, n_neg):
    vol = signed_volumes(X, tets)
    dih = dihedral_angles(X, tets)
    E = np.concatenate([np.linalg.norm(X[tets[:, i]] - X[tets[:, j]], axis=1)[:, None]
                        for i in range(4) for j in range(i + 1, 4)], 1)
    hmean = E.mean(1)
    cen = X[tets].mean(1)
    h_t, d = target_size(cen, p)
    # 面相邻单元的尺寸比
    fmap = {}
    for e, t in enumerate(tets):
        for f in ((t[0], t[1], t[2]), (t[0], t[1], t[3]), (t[0], t[2], t[3]), (t[1], t[2], t[3])):
            fmap.setdefault(tuple(sorted(f)), []).append(e)
    pairs = np.array([v for v in fmap.values() if len(v) == 2])
    ratio = np.maximum(hmean[pairs[:, 0]], hmean[pairs[:, 1]]) / np.minimum(hmean[pairs[:, 0]], hmean[pairs[:, 1]])
    fine = d <= p["R_in"]
    top = np.isclose(X[:, 2], L)
    rtop = np.hypot(X[top, 0] - L / 2, X[top, 1] - L / 2)
    pct = lambda v, q: float(np.percentile(v, q))
    rep = dict(name=name, params=p, nodes=int(len(X)), tets=int(len(tets)), t_generate=t_gen,
               n_flipped_orientation=n_neg, min_volume=float(vol.min()),
               dihedral_min=float(dih.min()), dihedral_max=float(dih.max()),
               frac_dihedral_outside_20_160=float(((dih.min(1) < 20) | (dih.max(1) > 160)).mean()),
               minSICN_min=float(gq["minSICN"].min()), minSICN_p1=pct(gq["minSICN"], 1),
               gamma_min=float(gq["gamma"].min()), gamma_p1=pct(gq["gamma"], 1), gamma_median=pct(gq["gamma"], 50),
               neighbor_size_ratio_p99=pct(ratio, 99), neighbor_size_ratio_max=float(ratio.max()),
               fine_region_tets=int(fine.sum()),
               fine_edge_mean_over_target=float(hmean[fine].mean() / p["h_local"]) if fine.any() else None,
               fine_edge_p5_p95_mm=[pct(E[fine].ravel() * 1e3, 5), pct(E[fine].ravel() * 1e3, 95)] if fine.any() else None,
               top_nodes_within_a=int((rtop <= p["a"] + 1e-9).sum()), top_nodes_within_Rin=int((rtop <= p["R_in"] + 1e-9).sum()))
    return rep


if __name__ == "__main__":
    names = sys.argv[1:] or list(CANDIDATES)
    for name in names:
        p = CANDIDATES[name]
        X, tets, gq, t_gen, n_neg = generate(name, p)
        np.savez(os.path.join(OUT, f"{name}.npz"), X=X, tets=tets, **{k: v for k, v in p.items() if k != "center_node"})
        rep = quality(name, X, tets, p, gq, t_gen, n_neg)
        json.dump(rep, open(os.path.join(OUT, f"{name}_quality.json"), "w"), indent=1)
        print(f"{name:8s} 节点 {rep['nodes']:5d}，四面体 {rep['tets']:6d}（加密区 {rep['fine_region_tets']}），{t_gen:.1f} s；"
              f"二面角 {rep['dihedral_min']:.1f}–{rep['dihedral_max']:.1f}°（20–160° 之外的单元 {rep['frac_dihedral_outside_20_160']:.1%}），"
              f"minSICN 最小 {rep['minSICN_min']:.3f}（1% 分位 {rep['minSICN_p1']:.3f}），gamma 最小 {rep['gamma_min']:.3f}；"
              f"相邻尺寸比 p99 {rep['neighbor_size_ratio_p99']:.2f} / 最大 {rep['neighbor_size_ratio_max']:.2f}；"
              f"加密区平均边长 / 目标 {rep['fine_edge_mean_over_target']:.2f}；顶面 r ≤ a 的节点 {rep['top_nodes_within_a']}；"
              f"方向翻转 {n_neg}", flush=True)
