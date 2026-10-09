"""工作网格 W50 的一次性流畅度检查：针已插入 50 mm 时的计算量（为第 3 步插入预估）。

W50 = A+ 的加密区（h_local 5 mm、半径 20 mm）沿竖直针道延长到 50 mm 深（meshes/step1_5/W50.npz）。
负载：针身上每 5 mm 一个材料点（BarycentricMapping 嵌入组织），共 10 个点，用 BilateralLagrangianConstraint 和纯运动学
针身上的对应点相连（30 行，相当于刺穿后沿针道的插入约束；间距 ≥ 单元尺寸，见 Martin 2023 的经验），针以 5 mm/s 沿轴向运动。
求解路线 R3（AsyncSparseLDLSolver），dt = 0.01；统计每步计算时间（中位数、p95、最大），换算成大致帧率（不含绘制）。
只是流畅度的粗查（2026-10-09 起的轻量验证方式），不是插入机制本身。
运行：python3 w50_insertion_timing.py → 终端输出
"""
import sys
sys.dont_write_bytecode = True
import os, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../step1"))
import numpy as np
import Sofa, Sofa.Simulation
import common, contact
import t1_vertex as V

contact.load_plugins()
L, dt = V.L, 0.01


def run(name="W50", depth=0.050, spacing=0.005, nsteps=300, skip=30):
    root = common.make_root(dt=dt, loop="free")
    tis = contact.add_tissue_with_surface(root, L=L, E=V.E, nu=V.NU, linear="async", mesh=common.load_mesh(name))
    n = int(round(depth / spacing))
    x0, y0 = L / 2 + 0.0007, L / 2 + 0.0004          # 略偏离中心线，避免点正好落在节点 / 面上
    pts = np.array([[x0, y0, L - (i + 0.5) * spacing] for i in range(n)])
    base0 = np.array([x0, y0, L + 0.02])
    local = pts - base0
    traj = lambda tt: (base0 + np.array([0, 0, -0.005 * tt]), np.array([0, 0, 0, 1.0]))
    nd, base, body, tip, _ = common.add_kinematic_needle(root, pose0=traj(0.0), local_points=local)
    root.addObject(common.NeedleDriver(base, traj, dt, children=[(body, local), (tip, local[-1:])], name="driver"))
    emb = tis["node"].addChild("Embedded")
    emo = emb.addObject("MechanicalObject", name="dofs", template="Vec3d", position=pts.tolist())
    emb.addObject("BarycentricMapping", input="@../dofs", output="@dofs")
    root.addObject("BilateralLagrangianConstraint", object1=body.getLinkPath(), object2=emo.getLinkPath(),
                   first_point=list(range(n)), second_point=list(range(n)))
    Sofa.Simulation.init(root)
    ts = []
    for k in range(nsteps):
        a = time.perf_counter(); Sofa.Simulation.animate(root, dt); ts.append(time.perf_counter() - a)
    Sofa.Simulation.unload(root)
    ts = 1e3 * np.array(ts[skip:])
    print(f"{name}：针插入 {depth*1e3:.0f} mm，{n} 个约束点（{3*n} 行）：每步计算 中位数 {np.median(ts):.2f} / p95 {np.percentile(ts, 95):.2f} / "
          f"最大 {ts.max():.2f} ms，约 {1e3/np.mean(ts):.0f} 帧/秒（不含绘制）", flush=True)


if __name__ == "__main__":
    run()
