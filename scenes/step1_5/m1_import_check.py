"""第 1.5 步 M1：gmsh 非均匀网格导入 SOFA 后的基础检查。见 docs/plan/step1_5_mesh_design.md §5。

对每个候选网格（meshes/step1_5/*.npz，由 scenes/mesh/make_graded_mesh.py 生成）：
  1. 拓扑：SOFA 拓扑容器里的四面体数、节点坐标和 npz 一致；有符号体积全为正；底面（z = 0）节点数；
  2. 表面：Tetra2TriangleTopologicalMapping 提取的顶面三角形，按顶点顺序算出的法向是否朝组织内部（−z），
     和规则网格（第 1 步约定 flipNormals=False 时朝内）一致；
  3. MeshGmshLoader 读 .msh（msh 2.2 ASCII）得到的节点和四面体与 npz 相同（集合意义下：坐标逐点一致、四面体逐个一致）；
  4. 小压深接触：a_eff = 10 mm 圆盘（点距 s = 标称 h_local），R1，0.2 mm/s 准静态压到 0.05 mm，
     记录轴向力、侧向力 / 轴向力、合力矩、受力点数（只看能否正常出力，侧向漂移的正式检查在 M3 ①）。
规则网格 h = 5 mm（swapping=True）作对照。
运行：python3 m1_import_check.py [网格名 ...] → results/step1_5/m1_import_check.json
"""
import sys
sys.dont_write_bytecode = True
import os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../step1"))
import numpy as np
import Sofa, Sofa.Simulation
import common, contact
import t1_vertex as V

contact.load_plugins()
L, LEN, dt, A = V.L, V.LEN, 0.01, 0.010
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../results/step1_5")
MESHDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../meshes/step1_5")
NAMES = ["A-", "A", "A+", "A_far16", "A_g05", "B", "E"]


def tet_volumes(x, tets):
    a, b, c, d = (x[tets[:, i]] for i in range(4))
    return np.einsum("ij,ij->i", np.cross(b - a, c - a), d - a) / 6.0


def topology_check(mesh):
    root = common.make_root(dt=dt, loop="free")
    tis = contact.add_tissue_with_surface(root, L=L, E=V.E, nu=V.NU, mesh=mesh) if mesh is not None else \
        contact.add_tissue_with_surface(root, n=17, L=L, E=V.E, nu=V.NU)
    Sofa.Simulation.init(root)
    X = tis["dofs"].position.array().copy()
    tets = np.array(tis["node"].topo.tetrahedra.value)
    tris = np.array(tis["surface"].tris.triangles.value)
    Sofa.Simulation.unload(root)
    r = dict(nodes=int(len(X)), tets=int(len(tets)), bottom_nodes=int(len(tis["bottom"])),
             min_volume=float(tet_volumes(X, tets).min()))
    if mesh is not None:
        r["same_positions"] = bool(np.allclose(X, mesh["X"], atol=1e-12))
        r["same_tets"] = bool((tets == mesh["tets"]).all())
    top = tris[np.all(np.isclose(X[tris][:, :, 2], L, atol=1e-9), axis=1)]
    nz = np.cross(X[top[:, 1]] - X[top[:, 0]], X[top[:, 2]] - X[top[:, 0]])[:, 2]
    r.update(top_triangles=int(len(top)), top_normals_inward=int((nz < 0).sum()), top_normals_outward=int((nz > 0).sum()),
             surface_triangles=int(len(tris)))
    return r


def loader_check(name, mesh):
    root = common.make_root()
    root.addObject("RequiredPlugin", name="Sofa.Component.IO.Mesh")
    ld = root.addObject("MeshGmshLoader", name="loader", filename=os.path.join(MESHDIR, name + ".msh"))
    Sofa.Simulation.init(root)
    P = np.array(ld.position.value); T = np.array(ld.tetrahedra.value)
    ntri = len(ld.triangles.value)
    Sofa.Simulation.unload(root)
    # 和 npz 比较：npz 只保留四面体用到的节点、按 gmsh 节点编号升序重新编号；msh 2.2 中节点按编号顺序写出
    ok_n = len(P) == len(mesh["X"]) and np.allclose(P, mesh["X"], atol=1e-12)
    key = lambda tt: np.sort(np.sort(tt, axis=1), axis=0)
    ok_t = len(T) == len(mesh["tets"]) and bool((np.sort(T, axis=1) == np.sort(mesh["tets"], axis=1)).all())
    vol = tet_volumes(P, T) if len(T) else np.array([0.0])
    return dict(loader_nodes=int(len(P)), loader_tets=int(len(T)), loader_triangles=int(ntri),
                loader_same_nodes=bool(ok_n), loader_same_tets=bool(ok_t), loader_min_volume=float(vol.min()))


def contact_check(mesh, h, gap=0.0001, v=0.0002, depth=0.00005):
    pts = contact.tip_patch_points(LEN, A, h, "disk")
    zf, T = V.piecewise(L + gap, [(gap / v, -v), (depth / v, -v)])
    root = common.make_root(dt=dt, loop="free")
    traj = lambda t: (np.array([L / 2, L / 2, zf(t) + LEN]), V.Q)
    nd, base, body, tip, local = common.add_kinematic_needle(root, length=LEN, n_body=V.NB, pose0=traj(0.0))
    pmo, pg = contact.add_tip_patch(nd, pts)
    root.addObject(common.NeedleDriver(base, traj, dt, children=[(body, local), (tip, local[-1:]), (pmo, pts)], name="driver"))
    root.addObject("CollisionLoop")
    tis = contact.add_tissue_with_surface(root, L=L, E=V.E, nu=V.NU, linear="ldl", mesh=mesh) if mesh is not None else \
        contact.add_tissue_with_surface(root, n=17, L=L, E=V.E, nu=V.NU, linear="ldl")
    _, bg = contact.add_needle_geometry(nd, body, tip, V.NB)
    contact.add_tip_contact(root, pg, bg, tis, distance=0.02)
    Sofa.Simulation.init(root)
    for k in range(int(round(T / dt))):
        Sofa.Simulation.animate(root, dt)
    lam = np.array(root.csolver.constraintForces.value)
    z_tip = zf(T)
    F, M, _ = common.needle_wrench([pmo], lam, dt, np.array([L / 2, L / 2, z_tip]))
    gs = int(root.csolver.currentIterations.value)
    Sofa.Simulation.unload(root)
    d = L - z_tip
    return dict(n_points=int(len(pts)), k=float(F[2] / d), lateral_over_axial=float(np.hypot(F[0], F[1]) / abs(F[2])),
                M_tip=M.tolist(), cop_mm=(np.array([-M[1], M[0]]) / F[2] * 1e3).tolist(),
                n_active=int((lam > 1e-12 * max(lam.max(), 1e-30)).sum()), gs_iter=gs)


if __name__ == "__main__":
    out = []
    kS = 2 * A * V.E / (1 - V.NU ** 2)
    for name in (sys.argv[1:] or ["regular_h5"] + NAMES):
        mesh = None if name == "regular_h5" else common.load_mesh(name)
        h = 0.005 if mesh is None else float(mesh["h_local"])
        r = dict(name=name)
        r.update(topology_check(mesh))
        if mesh is not None:
            r.update(loader_check(name, mesh))
        r.update(contact_check(mesh, h))
        out.append(r)
        print(f"{name:10s} 节点 {r['nodes']:5d} 四面体 {r['tets']:6d} 底面节点 {r['bottom_nodes']:3d} 最小体积 {r['min_volume']:.2e}；"
              f"顶面三角形 {r['top_triangles']}（法向朝内 {r['top_normals_inward']} / 朝外 {r['top_normals_outward']}）"
              + (f"；与 npz 一致 {r['same_positions']}/{r['same_tets']}；MeshGmshLoader：节点 {r['loader_nodes']} 四面体 {r['loader_tets']} "
                 f"三角形 {r['loader_triangles']}，与 npz 一致 {r['loader_same_nodes']}/{r['loader_same_tets']}" if mesh is not None else "")
              + f"；圆盘 {r['n_points']} 点：k = {r['k']:.2f} N/m（k/k_S = {r['k']/kS:.3f}），侧向/轴向 {r['lateral_over_axial']:.1e}，"
                f"压力中心偏移 ({r['cop_mm'][0]:+.2f}, {r['cop_mm'][1]:+.2f}) mm，受力点 {r['n_active']}，GS {r['gs_iter']}", flush=True)
        json.dump(out, open(os.path.join(OUT, "m1_import_check" + ("_extra" if len(sys.argv) > 1 else "") + ".json"), "w"), indent=1)
