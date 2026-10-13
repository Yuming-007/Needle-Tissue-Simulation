"""第 2 步：刺穿机制的 GUI 演示场景（由用户在自己的屏幕上运行）。场景搭建和 s2_puncture.py 相同。

看点：
  ① 针尖压下表面形成漏斗形凹坑（刺穿前，绿色线 = 针尖–表面接触对）；
  ② 力达到阈值（0.22 N，约 5.1 mm 压深）时刺穿：针尖接触消失，表面在约 0.5 s 内平滑回弹（组织阻尼 β = 0.1）；
  ③ 刺穿后沿针道每 5 mm 出现一个耦合点（蓝色小球 = 组织上的耦合点，紫色 = 它在针身上的投影，两者几乎重合）；
  ④ 插到 30 mm 停 1 s，再退回：耦合点在针尖经过时逐个消失，全部消失后回到刺穿前状态。整个过程循环一次约 13 s（仿真时间）。
视角：斜前上方的相机（add_camera，可以用鼠标旋转 / 缩放），窗口 1280×960（ViewerSetting），用于录屏。
终端：每 20 步打印压深、针上的轴向力、耦合点数、状态；每 500 帧打印帧率。
方案（参数 1，见 docs/plan/step2_results.md §3）：
  A    ：实时求解器（AsyncSparseLDLSolver）+ 组织 Rayleigh 刚度阻尼 β = 0.1（默认，第 2 步采用的模型）；
  C    ：离线求解器（SparseLDLSolver，每步重新分解）+ β = 0.1 —— 和离线验证完全相同，结果和 A 几乎一样，但慢（约 10 帧/秒）；
  orig ：实时求解器，不加阻尼 —— 刺穿后几步内发散（只用来看问题本身）。
运行（项目目录下）：RUNSOFA_OPTS="-g glfw" tools/run_gui.sh scenes/step2/s2_gui_scene.py [A|C|orig]   （按空格开始）
"""
import sys
sys.dont_write_bytecode = True
import os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../step1_5"))
import numpy as np
import Sofa, Sofa.Core
import common
import s2_puncture as S
from m3_gui_scene import SmoothMonitor


class StatePrinter(Sofa.Core.Controller):
    def __init__(self, s, **kw):
        super().__init__(**kw)
        self.s, self.k, self.punctured = s, 0, False

    def onAnimateEndEvent(self, e):
        self.k += 1
        s, root = self.s, self.s["root"]
        lam = np.array(root.csolver.constraintForces.value)
        z = S.L - (np.asarray(s["base"].position.value)[0][2] - S.LEN)
        p_base = np.asarray(s["base"].position.value)[0][:3]
        F, _, _ = common.needle_wrench([s["pmo"], s["body"]], lam, s["dt"], p_base)
        n = len(S.coupling_points(s["body"]))
        if (n > 0) != self.punctured:
            self.punctured = n > 0
            print(f"[状态] 第 {self.k} 步，压深 {z*1e3:+.2f} mm：{'刺穿 → 插入状态' if n else '耦合点已全部删除 → 回到刺穿前状态'}",
                  flush=True)
        if self.k % 20 == 0:
            state = "插入" if n else ("接触（压凹）" if F[2] > 1e-6 else "未接触")   # 有约束行不等于有接触力（DistanceFilter 2 cm 内就建行）
            print(f"[针] 第 {self.k} 步：压深 {z*1e3:+6.2f} mm，轴向力 {F[2]:.3f} N，侧向 {np.hypot(F[0], F[1]):.3f} N，"
                  f"耦合点 {n} 个，状态：{state}", flush=True)


def add_camera(root, look_at=(S.L / 2, S.L / 2, 0.05), azimuth=25.0, elevation=20.0, distance=0.17,
               resolution=(1280, 960)):
    """斜视相机：从 look_at 出发，朝 -y 方向偏 azimuth 度、仰角 elevation 度、距离 distance（m）的位置看向 look_at，
    世界 z 轴朝上。另设窗口分辨率和深色纯色背景。只影响显示。"""
    from scipy.spatial.transform import Rotation
    a, e = np.radians(azimuth), np.radians(elevation)
    c = np.asarray(look_at, float)
    pos = c + distance * np.array([np.sin(a) * np.cos(e), -np.cos(a) * np.cos(e), np.sin(e)])
    f = (c - pos) / np.linalg.norm(c - pos)                   # 视线方向（相机局部 -z）
    r = np.cross(f, [0.0, 0.0, 1.0]); r /= np.linalg.norm(r)  # 右（局部 x）
    u = np.cross(r, f)                                        # 上（局部 y）
    q = Rotation.from_matrix(np.c_[r, u, -f]).as_quat()       # [x, y, z, w]，和 SOFA 的四元数顺序相同
    root.addObject("RequiredPlugin", name="Sofa.Component.Setting")
    root.addObject("ViewerSetting", resolution=list(resolution))
    root.addObject("BackgroundSetting", color=[0.12, 0.12, 0.15, 1.0])
    cam = root.addObject("InteractiveCamera", name="camera", position=pos.tolist(), orientation=q.tolist(),
                         lookAt=c.tolist(), distance=distance, widthViewport=resolution[0], heightViewport=resolution[1])
    # GLFW 界面启动时会把相机重置成"看全整个场景"（包括很长的针），所以动画开始后的前几步再设一次
    root.addObject(CameraSetter(cam, pos, q, c, distance, name="cameraSetter"))


class CameraSetter(Sofa.Core.Controller):
    def __init__(self, cam, pos, q, c, distance, steps=3, **kw):
        super().__init__(**kw)
        self.cam, self.v, self.n = cam, (pos.tolist(), q.tolist(), c.tolist(), float(distance)), steps

    def onAnimateBeginEvent(self, e):
        if self.n > 0:
            self.n -= 1
            self.cam.position.value, self.cam.orientation.value, self.cam.lookAt.value, self.cam.distance.value = self.v


SCHEMES = {"A": ("async", S.BETA, S.THRESHOLD), "C": ("ldl", S.BETA, S.THRESHOLD), "orig": ("async", 0.0, S.THRESHOLD)}


def createScene(root):
    args = [a for a in sys.argv[1:] if not a.endswith(".py")]
    scheme = args[0] if args else "A"
    linear, beta, thr = SCHEMES[scheme]
    s = S.build(root, dt=0.01, linear=linear, threshold=thr, rayleigh_stiffness=beta)
    tis = s["tis"]
    root.addObject("VisualStyle", displayFlags="showVisualModels")      # 耦合点由插件绘制（drawcollision），不受显示开关影响
    root.algo.drawPointsScale = 0.001
    tis["fem"].drawing = False
    vis = tis["surface"].addChild("Visual")            # 显示模型放在单独的子节点（第 1 步的教训）
    vis.addObject("OglModel", name="visual", src="@../tris", color=[0.9, 0.6, 0.6, 0.45])
    vis.addObject("IdentityMapping", input="@../dofs", output="@visual")
    common.add_needle_visual(root.Needle, S.LEN, radius=0.0004)
    add_camera(root)
    root.addObject(SmoothMonitor(name="monitor"))
    root.addObject(StatePrinter(s, name="state"))
    print(f"[场景] 第 2 步刺穿，方案 {scheme}：求解器 {linear}，阻尼 β = {beta}，W50，单点针尖，阈值 {thr} N，耦合点间距 {S.TIP_DIST*1e3:.0f} mm，"
          f"下插到 {S.DEPTH*1e3:.0f} mm 再退回", flush=True)
    return root
