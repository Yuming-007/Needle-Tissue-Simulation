"""第 1 步：A2 轻量检查——SOFA 自带的拉格朗日接触 vs 插件 A1，同一个顶点接触工况。

A2 配置参考 probe–tissue 教程（docs/resources/inventory.md）：碰撞流水线 + BruteForceBroadPhase + BVHNarrowPhase
+ LocalMinDistance（alarmDistance、contactDistance）+ 接触管理器（FrictionContactConstraint，mu = 0）；
针尖 PointCollisionModel，组织表面 TriangleCollisionModel / LineCollisionModel / PointCollisionModel。
contactDistance 会让针尖停在离表面 contactDistance 的地方，所以按"中心节点的实际压深"比较力。
针以 1 mm/s 压到 5 mm（R1，dt = 0.01）。
运行：python3 t1_a2_check.py → results/step1/t1_a2_check.json
"""
import sys
sys.dont_write_bytecode = True
import os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
import Sofa, Sofa.Simulation
import common, contact
import t1_vertex as V

L, dt, v, gap, depth = V.L, 0.01, 0.001, 0.001, 0.005
CD = 1e-5


def run_a2():
    zf, T = V.piecewise(L + gap, [(gap / v, -v), (depth / v, -v)])
    root = common.make_root(dt=dt, loop="free")
    root.addObject("CollisionPipeline")
    root.addObject("BruteForceBroadPhase"); root.addObject("BVHNarrowPhase")
    root.addObject("LocalMinDistance", alarmDistance=0.002, contactDistance=CD)
    root.addObject("CollisionResponse", response="FrictionContactConstraint", responseParams="mu=0")
    traj = lambda t: (np.array([L / 2, L / 2, zf(t) + V.LEN]), V.Q)
    nd, base, body, tip, local = common.add_kinematic_needle(root, length=V.LEN, n_body=V.NB, pose0=traj(0.0))
    root.addObject(common.NeedleDriver(base, traj, dt, children=[(body, local), (tip, local[-1:])], name="driver"))
    tip.getContext().addObject("PointCollisionModel", name="tipColl", group=1)
    tis = contact.add_tissue_with_surface(root, n=V.N, L=L, E=V.E, nu=V.NU, flip_normals=True)
    s = tis["surface"]
    for m in ["TriangleCollisionModel", "LineCollisionModel", "PointCollisionModel"]:
        s.addObject(m, group=2)
    Sofa.Simulation.init(root)
    F, nz = [], []
    for k in range(int(round(T / dt))):
        Sofa.Simulation.animate(root, dt)
        lam = np.array(root.csolver.constraintForces.value)
        Fw, _, _ = common.needle_wrench([body, tip], lam, dt, np.asarray(base.free_position.value)[0][:3])
        F.append(Fw.tolist()); nz.append(float(tis["dofs"].position.array()[V.C][2]))
    Sofa.Simulation.unload(root)
    return np.array(F), np.array(nz)


if __name__ == "__main__":
    Fa2, nz2 = run_a2()
    zf, T = V.piecewise(L + gap, [(gap / v, -v), (depth / v, -v)])
    r1 = V.run(dt, zf, T)
    Fa1 = np.array(r1["F_needle"]); nz1 = np.array(r1["node"])[:, 2]
    d1, d2 = L - nz1, L - nz2
    m = d2 >= 0.001
    Fi = np.interp(d2[m], d1, Fa1[:, 2])
    rel = Fa2[m, 2] / Fi - 1
    lat = np.abs(Fa2[m][:, :2]).max() / np.abs(Fa2[m, 2]).max()
    print(f"A2（SOFA 自带接触）vs A1（插件）：中心节点压深 ≥ 1 mm 时，同一压深下针根力相对差 最大 {np.abs(rel).max():.2e}，"
          f"平均 {rel.mean():+.2e}；A2 侧向力 / 轴向力 {lat:.2e}；A2 最终压深 {d2[-1]*1e3:.4f} mm（针尖压深 {depth*1e3:.1f} mm，"
          f"contactDistance {CD*1e3:.3f} mm）", flush=True)
    json.dump(dict(rel_max=float(np.abs(rel).max()), rel_mean=float(rel.mean()), lateral=float(lat),
                   F_a2=Fa2.tolist(), node_z_a2=nz2.tolist(), F_a1=Fa1.tolist(), node_z_a1=nz1.tolist()),
              open(os.path.join(V.OUT, "t1_a2_check.json"), "w"))
