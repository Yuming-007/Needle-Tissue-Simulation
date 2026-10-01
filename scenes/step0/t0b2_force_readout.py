"""0b-2 纯运动学刚性针：通过约束读针根力和力矩（任意方向）。见 docs/plan/step0_foundation.md §5.4。

场景：FreeMotionAnimationLoop + BlockGaussSeidelConstraintSolver；针 = common.add_kinematic_needle（mode="end"）。
质点：EulerImplicitSolver（二阶）+ SparseLDLSolver + LinearSolverConstraintCorrection（和组织以后的配置相同），
      质量 m，RestShapeSpringsForceField（各向同性刚度 k）把它拉向锚点（MechanicalObject 的 rest_position）。
约束：SOFA 自带的 BilateralLagrangianConstraint，object1 = 针上的点（针尖或针身点），object2 = 质点。
      J_1 = -I，J_2 = +I，δ = x2 - x1（BilateralLagrangianConstraint.inl:160-186, 214-226）
      ⇒ 质点受到 +λ/dt，针上的点受到 -λ/dt。
针以 v = 5 mm/s 沿方向 d 运动 2 s（10 mm），再停 0.5 s。

五种工况：
  ① 针竖直（针轴 -z），针尖连质点，d = -z（轴向）
  ② 针竖直，针尖连质点，d = +x（横向）
  ③ 针竖直，针尖连质点，d = (x - z)/√2（45° 斜向）
  ④ 针竖直，针尖和针身中点各连一个质点（k 不同），锚点预先偏移（针尖的向 +x 偏 5 mm，中点的向 +y 偏 5 mm），
     d = +x，使两点受力方向不同
  ⑤ 针绕 y 轴倾斜 30° 后重复 ② 和 ④
另外工况 ② 和 ④ 在 dt = 0.005 / 0.01 / 0.02 下重复。

参照（独立于 λ）：由质点自身的状态算出它受到的约束力 f_c = m·a + k·(x - 锚点)，
a = (v_{n+1} - v_n)/dt；针上对应的点受 -f_c；针根解析值 F = Σ(-f_c)，τ = Σ (x_i - p) × (-f_c)。
读数：common.needle_wrench（针身点 / 针尖的 J 和 λ）。
运行：python3 t0b2_force_readout.py → results/step0/t0b2_force_readout.json
"""
import sys
sys.dont_write_bytecode = True
import os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
import Sofa, Sofa.Core, Sofa.Simulation
import common

LENGTH, V, T_MOVE, T_HOLD, MASS = 0.10, 5e-3, 2.0, 0.5, 0.01
P0 = np.array([0.04, 0.04, 0.10])
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../results/step0/t0b2_force_readout.json")


def make_traj(q, d):
    d = np.asarray(d, float); d /= np.linalg.norm(d)

    def f(t):
        return P0 + V * min(t, T_MOVE) * d, q
    return f


def run(tilt_deg, d, attach, dt):
    """attach：[(子节点名 'tip' / 'body', 点号, 刚度 k, 锚点偏移[3]), ...]"""
    q = common.quat_from_axis_angle([0, 1, 0], np.radians(tilt_deg))
    traj = make_traj(q, d)
    root = common.make_root(dt=dt, loop="free")
    nd, base, body, tip, local = common.add_kinematic_needle(root, length=LENGTH, pose0=traj(0.0))
    root.addObject(common.NeedleDriver(base, traj, dt, children=[(body, local), (tip, local[-1:])], name="driver"))
    R0 = common.quat_to_mat(q)
    parts = []
    for j, (which, idx, k, off) in enumerate(attach):
        mo_n = tip if which == "tip" else body
        loc = local[-1] if which == "tip" else local[idx]
        x0 = P0 + R0 @ loc
        pn = root.addChild(f"Particle{j}")
        pn.addObject("EulerImplicitSolver")
        pn.addObject("SparseLDLSolver", template="CompressedRowSparseMatrixMat3x3d")
        pm = pn.addObject("MechanicalObject", template="Vec3d", position=[x0.tolist()],
                          rest_position=[(x0 + np.asarray(off)).tolist()])
        pn.addObject("UniformMass", totalMass=MASS)
        pn.addObject("RestShapeSpringsForceField", stiffness=k, points=[0])
        pn.addObject("LinearSolverConstraintCorrection")
        root.addObject("BilateralLagrangianConstraint", name=f"bilat{j}", object1=mo_n.getLinkPath(),
                       object2=pm.getLinkPath(), first_point=[0 if which == "tip" else idx], second_point=[0])
        parts.append(dict(mo=pm, k=k, anchor=x0 + np.asarray(off), needle_mo=mo_n, needle_idx=0 if which == "tip" else idx))
    Sofa.Simulation.init(root)
    nsteps = int(round((T_MOVE + T_HOLD) / dt))
    v_prev = [p["mo"].velocity.array()[0].copy() for p in parts]
    rec = dict(F=[], tau=[], F_ref=[], tau_ref=[], gap=[], t=[])
    for n in range(nsteps):
        Sofa.Simulation.animate(root, dt)
        p_base = np.asarray(base.free_position.value)[0][:3]
        F, tau, per = common.needle_wrench([body, tip], root.csolver.constraintForces.value, dt, p_base)
        F_ref, tau_ref, gap = np.zeros(3), np.zeros(3), 0.0
        for j, p in enumerate(parts):
            x = p["mo"].position.array()[0]; v = p["mo"].velocity.array()[0]
            f_c = MASS * (v - v_prev[j]) / dt + p["k"] * (x - p["anchor"])   # 质点受到的约束力
            v_prev[j] = v.copy()
            xn = p["needle_mo"].free_position.array()[p["needle_idx"]]
            F_ref += -f_c
            tau_ref += np.cross(xn - p_base, -f_c)
            gap = max(gap, np.linalg.norm(x - xn))                           # 约束误差（质点和针上的点的距离）
        for key, val in zip(["F", "tau", "F_ref", "tau_ref", "gap", "t"], [F, tau, F_ref, tau_ref, gap, (n + 1) * dt]):
            rec[key].append(val.tolist() if isinstance(val, np.ndarray) else val)
    F, Fr, T, Tr = (np.array(rec[k]) for k in ["F", "F_ref", "tau", "tau_ref"])
    scaleF = np.abs(Fr).max(); scaleT = np.abs(Tr).max()
    i_end = int(round(T_MOVE / dt)) - 1                                       # 运动结束时
    return dict(rec=rec, err_F=float(np.abs(F - Fr).max() / scaleF), err_tau=float(np.abs(T - Tr).max() / max(scaleT, 1e-30)),
                gap=float(max(rec["gap"])), F_end=F[i_end].tolist(), tau_end=T[i_end].tolist(),
                F_ref_end=Fr[i_end].tolist(), tau_ref_end=Tr[i_end].tolist())


if __name__ == "__main__":
    tipK = [("tip", 0, 10.0, [0, 0, 0])]
    two = [("tip", 0, 10.0, [0.005, 0, 0]), ("body", 10, 20.0, [0, 0.005, 0])]
    s2 = 1 / np.sqrt(2)
    cases = [("1 axial", 0, [0, 0, -1], tipK), ("2 lateral", 0, [1, 0, 0], tipK), ("3 oblique45", 0, [s2, 0, -s2], tipK),
             ("4 two-point", 0, [1, 0, 0], two), ("5a tilted lateral", 30, [1, 0, 0], tipK), ("5b tilted two-point", 30, [1, 0, 0], two)]
    out = []
    for name, tilt, d, att in cases:
        for dt in ([0.005, 0.01, 0.02] if name.startswith(("2", "4")) else [0.01]):
            r = run(tilt, d, att, dt)
            r.update(case=name, dt=dt)
            out.append(r)
            print(f"{name:20s} dt={dt:.3f}: 力误差 {r['err_F']:.1e}  力矩误差 {r['err_tau']:.1e}  约束误差 {r['gap']:.1e} m | "
                  f"运动结束时 F = {np.round(r['F_end'], 5)} N（参照 {np.round(r['F_ref_end'], 5)}），"
                  f"τ = {np.round(r['tau_end'], 6)} N·m（参照 {np.round(r['tau_ref_end'], 6)}）", flush=True)
    for r in out:
        r["rec"] = {k: v[::5] for k, v in r["rec"].items()}                 # 保存时降采样
    json.dump(out, open(RES, "w"))
