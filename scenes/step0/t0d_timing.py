"""0d 第一阶段：计时筛选。见 docs/plan/step0_foundation.md §7.1。

组织：8 cm 立方体，规则网格 n = 7 / 9 / 11，FTCfix polar，ν = 0.45，E = 5 kPa，ρ = 1000，底面固定，重力关闭。
动画循环：FreeMotionAnimationLoop + BlockGaussSeidelConstraintSolver（common.make_root(loop="free")）。
负载：
  W0：只有组织（为了让 FreeMotion 路径一致，针仍然在场，但不加约束）。
  W2：组织里沿竖直线、每个单元一个材料点（BarycentricMapping 嵌入；线的位置偏离网格节点和面，
      使点落在四面体内部），针在组织内约 7 cm；纯运动学的针身上对应的点用 BilateralLagrangianConstraint
      连到这些材料点（每点 3 行）。针以 5 mm/s 沿轴向向下运动（相当于完全粘住的针道）。
求解路线：
  R1   SparseLDLSolver + LinearSolverConstraintCorrection（基线：每步重新分解）
  R3   AsyncSparseLDLSolver + LinearSolverConstraintCorrection（另一个线程做分解，用上一次的结果）
  R2   PCGLinearSolver + WarpPreconditioner（RotationMatrixSystem，每 assemblingRate 步重新组装 / 分解）
       + LinearSolverConstraintCorrection。注意：FTC 力场没有实现 BaseRotationFinder，所以没有转动修正，
       预条件矩阵只是"每隔几步刷新一次的旧矩阵"。
  R4   SparseLDLSolver（自由运动）+ PrecomputedConstraintCorrection(rotations=True)（约束柔度预先计算）
计时：Sofa.Simulation.animate 的墙钟时间，跑 60 步，去掉前 10 步，取中位数和 90 分位数。
     计时只依赖每步的计算量，和 dt 无关；RTF = dt / t_step 对 dt = 0.01 和 0.02 分别给出。
运行：python3 t0d_timing.py → results/step0/t0d_timing.json
"""
import sys
sys.dont_write_bytecode = True
import os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
import Sofa, Sofa.Core, Sofa.Simulation
import common

L, E, NU, RHO, DT = 0.08, 5000.0, 0.45, 1000.0, 0.01
NEEDLE_IN = 0.07           # 针在组织内的长度
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../results/step0/t0d_timing.json")
SCRATCH = "/tmp/claude-1000/-home-yuming-sofa-needle-project/982afcea-fafe-44ba-8bab-2e129c167ce6/scratchpad/precomputed"


def add_solver(t, route, assembling_rate=10, pcg_iters=50, pcg_tol=1e-10):
    t.addObject("EulerImplicitSolver", name="odesolver", rayleighMass=0.0, rayleighStiffness=0.0)
    if route in ("R1", "R4"):
        t.addObject("SparseLDLSolver", name="linsolver", template="CompressedRowSparseMatrixMat3x3d")
    elif route == "R3":
        t.addObject("AsyncSparseLDLSolver", name="linsolver", template="CompressedRowSparseMatrixMat3x3d")
    elif route == "R2":
        t.addObject("PreconditionedMatrixFreeSystem", name="matrixFreeSystem", assemblingRate=1)
        t.addObject("PCGLinearSolver", name="linsolver", iterations=pcg_iters, tolerance=pcg_tol,
                    linearSystem="@matrixFreeSystem", preconditioner="@preconditioner")
        t.addObject("RotationMatrixSystem", name="rotationMatrix", assemblingRate=assembling_rate)
        t.addObject("WarpPreconditioner", name="preconditioner", linearSystem="@rotationMatrix", linearSolver="@initSolver")
        t.addObject("MatrixLinearSystem", template="CompressedRowSparseMatrixMat3x3d", name="system")
        t.addObject("SparseLDLSolver", name="initSolver", template="CompressedRowSparseMatrixMat3x3d", linearSystem="@system")


def build(n, route, load, root=None, **kw):
    h = L / (n - 1)
    root = common.make_root(dt=DT, loop="free", root=root)
    tissue_parent = root.addChild("TissueRoot")
    add_solver(tissue_parent, route, **kw)
    t, mo, _ = common.add_tissue(tissue_parent, L=L, n=n, E=E, nu=NU, rho=RHO)
    X = common.grid_nodes(L, n)
    t.addObject("FixedProjectiveConstraint", indices=np.where(np.isclose(X[:, 2], 0))[0].tolist())
    if route == "R4":
        os.makedirs(SCRATCH, exist_ok=True)
        t.addObject("PrecomputedConstraintCorrection", rotations=True, recompute=True, fileDir=SCRATCH)
    else:
        t.addObject("LinearSolverConstraintCorrection", linearSolver="@../linsolver")
    # 约束点：竖直线 x = L/2 + 0.31h，y = L/2 + 0.17h，z_i = L - (i + 0.5) h
    n_pts = int(NEEDLE_IN / h)
    xl, yl = L / 2 + 0.31 * h, L / 2 + 0.17 * h
    pts = np.array([[xl, yl, L - (i + 0.5) * h] for i in range(n_pts)])
    base0 = np.array([xl, yl, L + 0.02])
    local = pts - base0                                   # 针沿世界 -z，姿态为单位四元数，局部 = 世界偏移
    traj = lambda tt: (base0 + np.array([0, 0, -0.005 * tt]), np.array([0, 0, 0, 1.0]))
    nd, base, body, tip, _ = common.add_kinematic_needle(root, pose0=traj(0.0), local_points=local)
    root.addObject(common.NeedleDriver(base, traj, DT, children=[(body, local), (tip, local[-1:])], name="driver"))
    if load == "W2":
        emb = t.addChild("Embedded")
        emo = emb.addObject("MechanicalObject", name="dofs", template="Vec3d", position=pts.tolist())
        emb.addObject("BarycentricMapping", input="@../dofs", output="@dofs")
        root.addObject("BilateralLagrangianConstraint", object1=body.getLinkPath(), object2=emo.getLinkPath(),
                       first_point=list(range(n_pts)), second_point=list(range(n_pts)))
    return root, n_pts


def timeit(n, route, load, nsteps=60, skip=10, **kw):
    t0 = time.perf_counter()
    root, n_pts = build(n, route, load, **kw)
    Sofa.Simulation.init(root)
    t_init = time.perf_counter() - t0
    ts = []
    for k in range(nsteps):
        a = time.perf_counter()
        Sofa.Simulation.animate(root, DT)
        ts.append(time.perf_counter() - a)
    ts = np.array(ts[skip:])
    lam = np.array(root.csolver.constraintForces.value)
    Sofa.Simulation.unload(root)
    return dict(n=n, nodes=n ** 3, route=route, load=load, n_pts=n_pts if load == "W2" else 0,
                rows=3 * n_pts if load == "W2" else 0, t_init=t_init,
                t_med_ms=1e3 * float(np.median(ts)), t_p90_ms=1e3 * float(np.percentile(ts, 90)),
                lam_last=lam.tolist(), **kw)


if __name__ == "__main__":
    out = []
    for n in [7, 9, 11]:
        for route in ["R1", "R3", "R4"]:   # R2 不可行（见文档）
            for load in ["W0", "W2"]:
                if route == "R4" and load == "W0":
                    continue
                try:
                    r = timeit(n, route, load)
                except Exception as e:                    # 记录失败的组合，继续
                    r = dict(n=n, nodes=n ** 3, route=route, load=load, error=str(e))
                    print(f"n={n} {route} {load}: 失败 {e}", flush=True)
                    out.append(r); continue
                out.append(r)
                print(f"n={n:2d} ({n**3:4d} 节点) {route} {load} 约束行 {r['rows']:2d}: "
                      f"t_step 中位数 {r['t_med_ms']:7.2f} ms，90% {r['t_p90_ms']:7.2f} ms；"
                      f"RTF(dt=0.01) = {10/r['t_med_ms']:.2f}，RTF(dt=0.02) = {20/r['t_med_ms']:.2f}；初始化 {r['t_init']:.1f} s",
                      flush=True)
                json.dump(out, open(RES, "w"), indent=1)
