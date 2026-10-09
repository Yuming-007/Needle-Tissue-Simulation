"""第 1.5 步 M3：非均匀工作网格上的动画流畅度（GUI 实测）。见 docs/plan/step1_5_results.md "工作网格 M3 筛选"。

流畅度标准（用户 2026-10-03）：只要求 SOFA 动画播放流畅不卡，**不要求按真实速度播放**。
  判据（用户认可）：帧率约 ≥ 30 帧/秒，即每帧（计算 + 绘制）平均 ≤ 约 33 ms，且没有明显的长时间停顿。
场景：gmsh 非均匀网格（参数 1：网格名，默认 A+），FTCfix polar，ν = 0.45，E = 5 kPa，底面固定；
  AsyncSparseLDLSolver，dt = 0.01；针尖（参数 3）：point = 单点（默认，2026-10-09 起），disk = a_eff = 10 mm 圆盘
  （点距 = 标称 h_local）；插件单边接触（不刺穿）。
  针竖直对准顶面中心，在表面上方 5 mm 和最大压深（参数 2，默认 8 mm）之间以 5 mm/s（仿真时间）往复。
  默认 8 mm：在 A+ / B / B+ 的大压深包络之内（B+ 的 d_0.1 = 9.5 mm；第一版用 10 mm 时 B+ 每个周期都被压到翻转，
  求解器迭代暴增，每步 > 0.6 s）。
不做节拍等待：能算多快就播放多快（慢于真实时间也可以）。
终端每 500 帧打印：计算时间、帧时间（相邻两帧开始时刻之差 = 计算 + 绘制 + 其他）的平均 / 中位数 / p95 / 最大，
  帧率 = 1 / 平均帧时间，超过 33 ms 和 100 ms 的帧所占比例；每 100 步打印针尖合力和压深。

运行（由用户在自己的屏幕上运行）：
  RUNSOFA_OPTS="-g glfw" tools/run_gui.sh scenes/step1_5/m3_gui_scene.py W50 [最大压深 mm] [point|disk]   （按空格开始）
无界面初筛（Claude 用）：RUNSOFA_OPTS="-g batch -n 1000" tools/run_gui.sh scenes/step1_5/m3_gui_scene.py A+
"""
import sys
sys.dont_write_bytecode = True
import os, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../step1"))
import numpy as np
import Sofa, Sofa.Core
import common, contact
import t1_vertex as V


class SmoothMonitor(Sofa.Core.Controller):
    def __init__(self, **kw):
        super().__init__(**kw)
        self.k, self.t_begin, self.prev_begin = 0, None, None
        self.compute, self.frame = [], []

    def onAnimateBeginEvent(self, e):
        now = time.perf_counter()
        if self.prev_begin is not None:
            self.frame.append(now - self.prev_begin)
        self.prev_begin = self.t_begin = now

    def onAnimateEndEvent(self, e):
        self.compute.append(time.perf_counter() - self.t_begin)
        self.k += 1
        if self.k % 500 == 0 and len(self.frame) > 20:
            c = 1e3 * np.array(self.compute[20:]); f = 1e3 * np.array(self.frame[20:])
            q = lambda a: f"平均 {a.mean():.2f} / 中位数 {np.median(a):.2f} / p95 {np.percentile(a, 95):.2f} / 最大 {a.max():.2f} ms"
            print(f"[流畅度] 前 {self.k} 帧（去掉前 20 帧）：计算 {q(c)}；帧时间 {q(f)}；帧率 {1e3/f.mean():.1f} 帧/秒；"
                  f"超过 33 ms 的帧 {100*(f > 33.3).mean():.1f}%，超过 100 ms 的帧 {100*(f > 100).mean():.2f}%", flush=True)


class ForcePrinter(Sofa.Core.Controller):
    def __init__(self, root, pmo, base, dt, **kw):
        super().__init__(**kw); self.root, self.pmo, self.base, self.dt, self.k = root, pmo, base, dt, 0

    def onAnimateEndEvent(self, e):
        self.k += 1
        if self.k % 100 == 0:
            lam = np.array(self.root.csolver.constraintForces.value)
            F, _, _ = common.needle_wrench([self.pmo], lam, self.dt, np.asarray(self.base.free_position.value)[0][:3])
            d = V.L - (np.asarray(self.base.position.value)[0][2] - V.LEN)
            print(f"[针尖] 第 {self.k} 步：压深 {d*1e3:+6.2f} mm（负值 = 在表面之上），针尖合力 F_z = {F[2]:.3f} N", flush=True)


def createScene(root):
    args = [a for a in sys.argv[1:] if not a.endswith(".py")]
    name = args[0] if args else "A+"
    depth = float(args[1]) * 1e-3 if len(args) > 1 else 0.008     # 参数 2：最大压深（mm），默认 8 mm
    tip_kind = args[2] if len(args) > 2 else "point"                # 参数 3：point（默认，2026-10-09 起）或 disk
    mesh = common.load_mesh(name)
    h = float(mesh["h_local"])
    dt, v, top, a_eff = 0.01, 0.005, 0.005, 0.010
    period = 2 * (top + depth) / v

    def zf(t):
        ph = t % period
        s = v * ph if ph < period / 2 else (top + depth) - v * (ph - period / 2)
        return V.L + top - s
    common.make_root(dt=dt, loop="free", root=root)
    traj = lambda t: (np.array([V.L / 2, V.L / 2, zf(t) + V.LEN]), V.Q)
    nd, base, body, tip, local = common.add_kinematic_needle(root, length=V.LEN, n_body=V.NB, pose0=traj(0.0))
    pts = contact.tip_patch_points(V.LEN, a_eff, h, "disk") if tip_kind == "disk" else np.array([[0.0, 0.0, -V.LEN]])
    pmo, pg = contact.add_tip_patch(nd, pts)
    root.addObject(common.NeedleDriver(base, traj, dt, children=[(body, local), (tip, local[-1:]), (pmo, pts)], name="driver"))
    root.addObject("CollisionLoop")
    tis = contact.add_tissue_with_surface(root, L=V.L, E=V.E, nu=V.NU, linear="async", mesh=mesh)
    _, bg = contact.add_needle_geometry(nd, body, tip, V.NB)
    contact.add_tip_contact(root, pg, bg, tis, distance=0.02)
    root.addObject("VisualStyle", displayFlags="showVisualModels showBehaviorModels")
    tis["fem"].drawing = False
    # 显示模型放在单独的子节点（第 1 步的教训：放在表面节点里会替换接触用的 MO 和映射）
    vis = tis["surface"].addChild("Visual")
    vis.addObject("OglModel", name="visual", src="@../tris", color=[0.9, 0.6, 0.6, 0.6])
    vis.addObject("IdentityMapping", input="@../dofs", output="@visual")
    common.add_needle_visual(nd, V.LEN)
    root.addObject(SmoothMonitor(name="monitor"))
    root.addObject(ForcePrinter(root, pmo, base, dt, name="forces"))
    print(f"[场景] 网格 {name}：{len(mesh['X'])} 个节点，{len(mesh['tets'])} 个四面体，h_local = {h*1e3:.0f} mm，针尖 {len(pts)} 个点，"
          f"最大压深 {depth*1e3:.1f} mm，针尖 {tip_kind}", flush=True)
    return root
