"""第 3 步：插入导向的 GUI 演示场景。场景搭建和 s3_insertion.py 相同（实时求解器 R3 + β = 0.1）。

方案（参数 1）：
  lateral（默认，3B）：刺穿 → 插到 30 mm → 整根针沿 +x 平移 3 mm → 停 2 s → 平移回原位 → 退出（约 22 s 仿真时间）。
           看点：组织跟着针横向移动（离针越深、越靠近固定底面，移动越小），耦合点始终在针上；终端打印针上的 F_x。
  deep（3A）：插到 50 mm 再退回，耦合点增加到 10 个再逐个删除。
视角：正前方略高（视线沿 +y，x 方向的平移看得最清楚），窗口 1280×960；可以用鼠标旋转 / 缩放。
终端：每 20 步打印压深、横向平移、针上的力、耦合点数；每 500 帧打印帧率。
运行（项目目录下）：RUNSOFA_OPTS="-g glfw" tools/run_gui.sh scenes/step3/s3_gui_scene.py [lateral|deep] [平移幅度 mm，默认 3]   （按空格开始）
录屏用的是 lateral 10（10 mm 更容易看清；离线检查 6 / 10 mm 时仍然线性、没有单元变坏，见 step3_results.md）。
"""
import sys
sys.dont_write_bytecode = True
import os
HERE = os.path.dirname(os.path.abspath(__file__))
for d in ["..", "../step1", "../step1_5", "../step2"]:
    sys.path.insert(0, os.path.join(HERE, d))
import numpy as np
import Sofa, Sofa.Core
import common
import s2_puncture as S
import s2_gui_scene as G2
import s3_insertion as I
from m3_gui_scene import SmoothMonitor


class Printer(Sofa.Core.Controller):
    def __init__(self, s, zf, xf, **kw):
        super().__init__(**kw)
        self.s, self.zf, self.xf, self.k, self.punctured = s, zf, xf, 0, False

    def onAnimateEndEvent(self, e):
        self.k += 1
        s, root = self.s, self.s["root"]
        t = self.k * s["dt"]
        lam = np.array(root.csolver.constraintForces.value)
        F, _, _ = common.needle_wrench([s["pmo"], s["body"]], lam, s["dt"], np.asarray(s["base"].position.value)[0][:3])
        n = len(S.coupling_points(s["body"]))
        if (n > 0) != self.punctured:
            self.punctured = n > 0
            print(f"[状态] 第 {self.k} 步：{'刺穿 → 插入状态' if n else '耦合点已全部删除 → 回到刺穿前状态'}", flush=True)
        if self.k % 20 == 0:
            print(f"[针] 第 {self.k} 步：压深 {(S.L - self.zf(t))*1e3:+6.2f} mm，横向平移 {self.xf(t)*1e3:.2f} mm，"
                  f"F = ({F[0]*1e3:+.1f}, {F[1]*1e3:+.1f}, {F[2]*1e3:+.1f}) mN，耦合点 {n} 个", flush=True)


def createScene(root):
    args = [a for a in sys.argv[1:] if not a.endswith(".py")]
    mode = args[0] if args else "lateral"
    amp = float(args[1]) * 1e-3 if len(args) > 1 else I.LAT_AMP          # 平移幅度（mm），默认 3 mm
    zf, xf, T = I.trajectories(mode, amp)
    s = S.build(root, dt=0.01, linear="async", zf=zf, xf=xf)
    tis = s["tis"]
    root.addObject("VisualStyle", displayFlags="showVisualModels")   # 耦合点由插件绘制（drawcollision）
    root.algo.drawPointsScale = 0.001
    tis["fem"].drawing = False
    vis = tis["surface"].addChild("Visual")
    vis.addObject("OglModel", name="visual", src="@../tris", color=[0.9, 0.6, 0.6, 0.45])
    vis.addObject("IdentityMapping", input="@../dofs", output="@visual")
    common.add_needle_visual(root.Needle, S.LEN, radius=0.0004)
    G2.add_camera(root, look_at=(S.L / 2, S.L / 2, 0.045), azimuth=0.0, elevation=12.0, distance=0.17)
    root.addObject(SmoothMonitor(name="monitor"))
    root.addObject(Printer(s, zf, xf, name="printer"))
    print(f"[场景] 第 3 步插入导向，方案 {mode}：R3 + β = {S.BETA}，W50，单点针尖，阈值 {S.THRESHOLD} N，"
          f"耦合点间距 {S.TIP_DIST*1e3:.0f} mm，平移幅度 {amp*1e3:g} mm，总时长 {T:.1f} s", flush=True)
    return root
