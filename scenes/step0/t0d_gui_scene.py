"""0d GUI 目测场景：冻结配置下的组织 + 来回运动的纯运动学针 + W2 约束（21 行）。

冻结配置（0d）：FTCfix polar，ν = 0.45，E = 5 kPa，8 cm 立方体 n = 9（h = 10 mm），dt = 0.01 s，
EulerImplicitSolver（不加 Rayleigh）+ AsyncSparseLDLSolver（R3）+ LinearSolverConstraintCorrection，
FreeMotionAnimationLoop + BlockGaussSeidelConstraintSolver。
针沿轴向以 5 mm/s 往复运动（振幅 10 mm），带着组织里的 7 个约束点（完全粘住的针道）。

运行（在项目根目录）：
  export SOFA_ROOT=~/sofa/SOFA_v25.12.00_Linux
  export PYTHONPATH=$SOFA_ROOT/plugins/SofaPython3/lib/python3/site-packages:$SOFA_ROOT/plugins/STLIB/lib/python3/site-packages
  export LD_LIBRARY_PATH=$SOFA_ROOT/lib
  $SOFA_ROOT/bin/runSofa -l SofaPython3 scenes/step0/t0d_gui_scene.py
然后点 Animate（GLFW 界面按空格）。动画按真实时间播放（针以 5 mm/s 运动，全局时钟对齐）。
终端每 100 步打印计算时间、滞后和有效实时比；每 500 步打印平均 / 中位数 / p95 / 最大值。
参数（tools/run_gui.sh <场景> [参数...]）：R1 改用直接求解器作对比；debug 改用调试显示
（画出所有四面体和力学对象，绘图开销大，约 12 ms/帧，见 0d 结果 3）。
默认是轻量显示：只画组织外表面（Tetra2TriangleTopologicalMapping + OglModel）和针。
"""
import sys
sys.dont_write_bytecode = True
import os, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
import Sofa, Sofa.Core
import common
import t0d_timing as T


class FrameMonitor(Sofa.Core.Controller):
    """按真实时间播放并统计。

    节拍（全局时钟对齐，"固定时间步 + 追赶"）：第 k 步的目标开始时刻 = t0 + k·dt。
    提前了就等待；落后了就不等待，直接算下一步，这样偶尔超时的帧落下的时间可以被后面有余量的帧追回来。
    （第一版是"每帧不足 dt 就补足"，超时的帧落下的时间永远追不回来。）
    统计（每 500 步，去掉前 20 帧）：
      计算 = AnimateBegin→AnimateEnd 的纯计算时间；计算 + 绘制 = 相邻两帧之间扣除等待后的时间；
      吞吐量判据：计算 + 绘制的平均值 ≤ dt（长期能跟上真实时间）；
      抖动：滞后 = 实际开始时刻 - 目标开始时刻（≥ 0 表示落后），记录 p95 和最大值；
      有效实时比 = 仿真时间 / 墙钟时间。
    """

    def __init__(self, dt, **kw):
        super().__init__(**kw)
        self.dt, self.k, self.t0, self.t_begin, self.t_end = dt, 0, None, None, None
        self.steps, self.busy, self.lag = [], [], []

    def onAnimateBeginEvent(self, e):
        now = time.perf_counter()
        if self.t0 is None:
            self.t0 = now
        else:
            self.busy.append(now - self.t_begin)              # 上一帧：计算 + 绘制（不含等待）
        target = self.t0 + self.k * self.dt
        if now < target:
            time.sleep(target - now)
            now = time.perf_counter()
        self.lag.append(now - target)
        self.t_begin = now

    def onAnimateEndEvent(self, e):
        self.steps.append(time.perf_counter() - self.t_begin)
        self.k += 1
        if self.k % 100 == 0 and len(self.busy) >= 100:
            s = 1e3 * np.median(self.steps[-100:]); b = 1e3 * np.median(self.busy[-100:])
            wall = time.perf_counter() - self.t0
            print(f"[0d] 第 {self.k} 步：计算 {s:.2f} ms/步，计算 + 绘制 {b:.2f} ms（中位数）；"
                  f"当前滞后 {1e3*self.lag[-1]:.1f} ms；有效实时比 {self.k*self.dt/wall:.3f}", flush=True)
        if self.k % 500 == 0 and len(self.busy) > 20:
            st = 1e3 * np.array(self.steps[20:]); bu = 1e3 * np.array(self.busy[20:]); lg = 1e3 * np.array(self.lag[20:])
            q = lambda a: f"平均 {a.mean():.2f} / 中位数 {np.median(a):.2f} / p95 {np.percentile(a, 95):.2f} / 最大 {a.max():.2f} ms"
            wall = time.perf_counter() - self.t0
            print(f"[0d 统计] 前 {self.k} 步（去掉前 20 帧）：计算 {q(st)}；计算 + 绘制 {q(bu)}；"
                  f"计算 + 绘制超过 dt 的帧 {100*(bu > 1e3*self.dt).mean():.1f}%；"
                  f"滞后 p95 {np.percentile(lg, 95):.1f} / 最大 {lg.max():.1f} ms；有效实时比 {self.k*self.dt/wall:.3f}", flush=True)


def createScene(root):
    args = sys.argv[1:]
    route = "R1" if "R1" in args else "R3"
    debug = "debug" in args
    root, n_pts = T.build(9, route, "W2", root=root)
    # 往复运动：三角波，振幅 10 mm，速度 5 mm/s
    drv = root.driver
    base0 = np.asarray(root.Needle.base.position.value)[0][:3].copy()
    amp, v = 0.010, 0.005

    def traj(t):
        period = 2 * amp / v
        ph = t % period
        s = v * ph if ph < period / 2 else amp - v * (ph - period / 2)
        return base0 + np.array([0, 0, -s]), np.array([0, 0, 0, 1.0])
    drv.traj = traj
    tissue = root.TissueRoot.Tissue
    if debug:
        root.addObject("VisualStyle", displayFlags="showBehaviorModels showForceFields showInteractionForceFields showWireframe")
        root.Needle.Body.dofs.showObject = True
        root.Needle.Body.dofs.showObjectScale = 4
        tissue.Embedded.dofs.showObject = True
        tissue.Embedded.dofs.showObjectScale = 6
    else:
        root.addObject("VisualStyle", displayFlags="showVisualModels showBehaviorModels")
        tissue.fem.drawing = False                      # 不画每个四面体
        surf = tissue.addChild("Surface")               # 组织外表面
        surf.addObject("TriangleSetTopologyContainer", name="tris")
        surf.addObject("TriangleSetTopologyModifier")
        surf.addObject("Tetra2TriangleTopologicalMapping", input="@../topo", output="@tris")
        surf.addObject("OglModel", name="visual", src="@tris", color=[0.9, 0.6, 0.6, 0.5])
        surf.addObject("IdentityMapping", input="@../dofs", output="@visual")
        root.Needle.Body.dofs.showObject = True         # 针身 7 个点 + 针根：画成点（开销很小）
        root.Needle.Body.dofs.showObjectScale = 4
        root.Needle.Body.dofs.drawMode = 1
    root.addObject(FrameMonitor(root.dt.value, name="monitor"))
    print(f"[0d] GUI 场景：求解路线 {route}，显示 {'调试' if debug else '轻量'}，n = 9，约束点 {n_pts} 个（{3*n_pts} 行），"
          f"dt = {root.dt.value}", flush=True)
    return root
