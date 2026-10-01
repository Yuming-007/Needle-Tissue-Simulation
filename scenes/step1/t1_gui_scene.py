"""第 1 步 GUI 目测场景：针尖压组织表面（tenting），后退时组织回弹。

测试台：8 cm 立方体，n = 9，swapping=True，FTCfix polar，ν = 0.45，E = 5 kPa，底面固定；实时配置
（AsyncSparseLDLSolver，dt = 0.01）；插件 A1 的针尖–表面单边接触（不刺穿）。
针竖直对准顶面中心，在表面上方 5 mm 和压深 10 mm 之间以 5 mm/s 往复。按真实时间播放（全局时钟对齐）。
终端每 100 步打印针根力（向上为正）、压深和计时。

运行：RUNSOFA_OPTS="-g glfw" tools/run_gui.sh scenes/step1/t1_gui_scene.py   （按空格开始）
"""
import sys
sys.dont_write_bytecode = True
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../step0"))
import numpy as np
import Sofa, Sofa.Core
import common, contact
import t1_vertex as V
from t0d_gui_scene import FrameMonitor


class ForcePrinter(Sofa.Core.Controller):
    def __init__(self, root, h, dt, **kw):
        super().__init__(**kw); self.root, self.h, self.dt, self.k = root, h, dt, 0

    def onAnimateEndEvent(self, e):
        self.k += 1
        if self.k % 100 == 0:
            lam = np.array(self.root.csolver.constraintForces.value)
            F, _, _ = common.needle_wrench([self.h["body"], self.h["tip"]], lam, self.dt,
                                           np.asarray(self.h["base"].free_position.value)[0][:3])
            d = V.L - self.h["tip"].position.array()[0][2]
            print(f"[1] 第 {self.k} 步：针尖压深 {d*1e3:+6.2f} mm（负值 = 在表面之上），针根力 {F[2]:.3f} N", flush=True)


def createScene(root):
    dt, v, top, depth = 0.01, 0.005, 0.005, 0.010
    period = 2 * (top + depth) / v

    def zf(t):
        ph = t % period
        s = v * ph if ph < period / 2 else (top + depth) - v * (ph - period / 2)
        return V.L + top - s
    common.make_root(dt=dt, loop="free", root=root)
    traj = lambda t: (np.array([V.L / 2, V.L / 2, zf(t) + V.LEN]), V.Q)
    nd, base, body, tip, local = common.add_kinematic_needle(root, length=V.LEN, n_body=V.NB, pose0=traj(0.0))
    root.addObject(common.NeedleDriver(base, traj, dt, children=[(body, local), (tip, local[-1:])], name="driver"))
    root.addObject("CollisionLoop")
    tis = contact.add_tissue_with_surface(root, n=V.N, L=V.L, E=V.E, nu=V.NU, linear="async")
    tg, bg = contact.add_needle_geometry(nd, body, tip, V.NB)
    contact.add_tip_contact(root, tg, bg, tis)
    root.addObject("VisualStyle", displayFlags="showVisualModels showBehaviorModels")
    tis["fem"].drawing = False
    # 显示模型必须放在单独的子节点里：一个节点只允许一个状态对象和一个映射，
    # 加在表面节点里会替换掉接触用的表面 MO 和映射（第一版的错误：接触柔度变成 0，λ = 1.8e308）。
    vis = tis["surface"].addChild("Visual")
    vis.addObject("OglModel", name="visual", src="@../tris", color=[0.9, 0.6, 0.6, 0.6])
    vis.addObject("IdentityMapping", input="@../dofs", output="@visual")
    common.add_needle_visual(nd, V.LEN)       # 针的显示模型（细圆柱，显示半径 0.8 mm）
    root.addObject(FrameMonitor(dt, name="monitor"))
    root.addObject(ForcePrinter(root, dict(base=base, body=body, tip=tip), dt, name="forces"))
    return root
