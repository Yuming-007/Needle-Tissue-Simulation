"""0b-1 纯运动学刚性针的运动学测试（没有组织）。见 docs/plan/step0_foundation.md §5.4。

场景：FreeMotionAnimationLoop + BlockGaussSeidelConstraintSolver（第 1 步要用的配置，这里没有约束，
只为走一遍自由运动 / 自由位置的计算路径）；针 = common.add_kinematic_needle（Rigid3，没有求解器）
+ 21 个针身点 + 1 个针尖点（RigidMapping）；common.NeedleDriver 每步设定位姿。
两种驱动写法：mode="end"（pos = 指令(t_{n+1})，freeVel = 0）和 mode="start"（pos = 指令(t_n)，freeVel = 步内平均速度）。
第一版只设定针根：子节点的 free_position 落后一步（第一步还是初始化时的局部坐标，误差约 0.1 m），
原因见 common.NeedleDriver 的说明；现在控制器同时写入子节点。
第二版 mode="start" + 转动时第一步误差 5e-5 m：针初始化时是竖直的，而轨迹 T2 一开始就倾斜 30°，
SOFA 用映射的雅可比（基于上一次更新的姿态）重算子节点的 freeVel（FreeMotionAnimationLoop.cpp:411），
第一步的雅可比因此是错的；现在针初始化在轨迹的起始位姿（pose0）。

轨迹：
  T1(v)：针竖直向下（针轴 -z），以 v = 1/5/10 mm/s 插入 80 mm，停 0.5 s，再以 v 拔出 20 mm。
  T2：  针先绕 y 轴倾斜 30°，然后针根沿 x 以 2 mm/s 平移，同时绕（过针根的）y 轴以 0.1 rad/s 转动，持续 2 s。
每步结束后读：针根 position / free_position，针身点和针尖的 position / free_position，
和解析值 p(t) + R(t)·r_i 比较。约束求解用的是 free_position，所以它必须等于这一步结束时的指令位姿。
运行：python3 t0b1_kinematics.py → results/step0/t0b1_kinematics.json
"""
import sys
sys.dont_write_bytecode = True
import os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
import Sofa, Sofa.Core, Sofa.Simulation
import common

DT, LENGTH = 0.01, 0.10
P0 = np.array([0.04, 0.04, 0.10])
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../results/step0/t0b1_kinematics.json")


def traj_T1(v, depth=0.08, hold=0.5, back=0.02):
    t1 = depth / v; t2 = t1 + hold; t3 = t2 + back / v
    q = np.array([0, 0, 0, 1.0])

    def f(t):
        if t <= t1:
            s, vz = v * t, -v
        elif t <= t2:
            s, vz = depth, 0.0
        elif t <= t3:
            s, vz = depth - v * (t - t2), v
        else:
            s, vz = depth - back, 0.0
        return P0 + np.array([0, 0, -s]), q, np.array([0, 0, vz]), np.zeros(3)
    return f, t3


def traj_T2(vx=0.002, w=0.1, tilt=np.radians(30), T=2.0):
    def f(t):
        tt = min(t, T)
        p = P0 + np.array([vx * tt, 0, 0])
        q = common.quat_from_axis_angle([0, 1, 0], tilt + w * tt)
        on = 1.0 if t < T else 0.0
        return p, q, np.array([vx * on, 0, 0]), np.array([0, w * on, 0])
    return f, T


def rot_angle(q1, q2):
    """两个四元数之间的转角（rad）。"""
    d = abs(np.dot(q1 / np.linalg.norm(q1), q2 / np.linalg.norm(q2)))
    return 2 * np.arccos(min(1.0, d))


def run(traj, t_end, mode):
    root = common.make_root(dt=DT, loop="free")
    nd, base, body, tip, local = common.add_kinematic_needle(root, length=LENGTH, pose0=traj(0.0)[:2])
    drv = root.addObject(common.NeedleDriver(base, traj, DT, children=[(body, local), (tip, local[-1:])],
                                             mode=mode, name="driver"))
    Sofa.Simulation.init(root)
    nsteps = int(round(t_end / DT)) + 20          # 结束后再多走 20 步（静止）
    err = {k: 0.0 for k in ["base_pos", "base_free", "base_rot", "base_free_rot",
                            "body_pos", "body_free", "tip_pos", "tip_free"]}
    t_step = []
    for k in range(nsteps):
        t0 = time.perf_counter()
        Sofa.Simulation.animate(root, DT)
        t_step.append(time.perf_counter() - t0)
        p1, q1, _, _ = traj((k + 1) * DT)         # 这一步结束时的指令位姿
        R1 = common.quat_to_mat(q1)
        X1 = p1 + local @ R1.T                    # 针身点的解析位置
        b = np.asarray(base.position.value)[0]; bf = np.asarray(base.free_position.value)[0]
        if mode == "end":
            err["base_pos"] = max(err["base_pos"], np.linalg.norm(b[:3] - p1))
            err["base_rot"] = max(err["base_rot"], rot_angle(b[3:], q1))
        else:                                     # mode="start"：步末 position 仍是 指令(t_n)
            p0, q0, _, _ = traj(k * DT)
            err["base_pos"] = max(err["base_pos"], np.linalg.norm(b[:3] - p0))
            err["base_rot"] = max(err["base_rot"], rot_angle(b[3:], q0))
        err["base_free"] = max(err["base_free"], np.linalg.norm(bf[:3] - p1))
        err["base_free_rot"] = max(err["base_free_rot"], rot_angle(bf[3:], q1))
        bp = body.position.array(); bfp = body.free_position.array()
        tp = tip.position.array()[0]; tfp = tip.free_position.array()[0]
        if mode == "end":
            err["body_pos"] = max(err["body_pos"], np.abs(bp - X1).max())
            err["tip_pos"] = max(err["tip_pos"], np.abs(tp - X1[-1]).max())
        err["body_free"] = max(err["body_free"], np.abs(bfp - X1).max())
        err["tip_free"] = max(err["tip_free"], np.abs(tfp - X1[-1]).max())
    t_now = float(root.time.value)
    return dict(err=err, steps=nsteps, time_drift=t_now - nsteps * DT,
                t_step_median_ms=1e3 * float(np.median(t_step)))


if __name__ == "__main__":
    out = []
    cases = [("T1", v) for v in [1e-3, 5e-3, 10e-3]] + [("T2", None)]
    for mode in ["end", "start"]:
        for name, v in cases:
            traj, t_end = traj_T1(v) if name == "T1" else traj_T2()
            r = run(traj, t_end, mode)
            r.update(mode=mode, traj=name, v=v)
            out.append(r)
            e = r["err"]
            print(f"mode={mode:5s} {name} v={'-' if v is None else f'{v*1e3:.0f}mm/s':>6s} steps={r['steps']:5d}: "
                  f"base pos {e['base_pos']:.1e} rot {e['base_rot']:.1e} | base free {e['base_free']:.1e} rot {e['base_free_rot']:.1e} | "
                  f"body pos {e['body_pos']:.1e} free {e['body_free']:.1e} | tip pos {e['tip_pos']:.1e} free {e['tip_free']:.1e} | "
                  f"time drift {r['time_drift']:.1e} s | {r['t_step_median_ms']:.3f} ms/step", flush=True)
    json.dump(out, open(RES, "w"), indent=1)
