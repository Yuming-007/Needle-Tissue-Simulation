"""项目公共工具：SOFA 环境、组织立方体构建。

用法：在脚本开头 `import common`（scenes/ 在 sys.path 中）。
约定（0c）：SI 单位（m, kg, s, N, Pa）。
"""
import sys
sys.dont_write_bytecode = True
import os

SOFA_ROOT = os.path.expanduser("~/sofa/SOFA_v25.12.00_Linux")
os.environ.setdefault("SOFA_ROOT", SOFA_ROOT)

import numpy as np
import Sofa
import Sofa.Core
import Sofa.Simulation
import SofaRuntime

SofaRuntime.importPlugin("Sofa.Component")

# 本项目的修补插件（plugins/NeedleSimFixes，编译到仓库外）
NEEDLESIM_PLUGIN_DIR = os.path.expanduser("~/sofa/needle_build/NeedleSimFixes")
SofaRuntime.PluginRepository.addFirstPath(NEEDLESIM_PLUGIN_DIR)
SofaRuntime.importPlugin("NeedleSimFixes")


def _vec3(v):
    return np.array(v if np.ndim(v) else [v, v, v], dtype=float)


def grid_nodes(L, n, origin=(0.0, 0.0, 0.0)):
    """RegularGridTopology 的节点顺序：x 最快，其次 y，最后 z。返回 (nx*ny*nz, 3) 坐标。

    L、n 可以是标量（立方体）或三元组（长方体：各方向尺寸 / 节点数）。
    """
    L3, n3, o = _vec3(L), _vec3(n).astype(int), _vec3(origin)
    sx, sy, sz = (np.linspace(o[i], o[i] + L3[i], n3[i]) for i in range(3))
    z, y, x = np.meshgrid(sz, sy, sx, indexing="ij")
    return np.stack([x.ravel(), y.ravel(), z.ravel()], axis=1)


def add_tissue(parent, L=0.08, n=9, E=5000.0, nu=0.3, rho=1000.0, method="polar",
               name="Tissue", origin=(0.0, 0.0, 0.0), ff_type="FTCfix"):
    """8 cm 立方体（默认），n×n×n 节点的规则网格切成四面体，共旋线弹性。

    parent 节点下不放求解器，由调用者决定（静力 / 动力）。
    返回 (tissue 节点, MechanicalObject, 力场)。
    """
    o, L3, n3 = _vec3(origin), _vec3(L), _vec3(n).astype(int)
    topo = parent.addChild(name + "Grid")
    topo.addObject("RegularGridTopology", name="hexa", n=n3.tolist(),
                   min=o.tolist(), max=(o + L3).tolist())
    t = parent.addChild(name)
    t.addObject("TetrahedronSetTopologyContainer", name="topo",
                position="@../%sGrid/hexa.position" % name)
    t.addObject("TetrahedronSetTopologyModifier", name="modifier")
    t.addObject("Hexa2TetraTopologicalMapping", input="@../%sGrid/hexa" % name,
                output="@topo", swapping=False)
    mo = t.addObject("MechanicalObject", name="dofs", template="Vec3d")
    # 质量：静力测试不用，但动力测试需要；密度 × 体积
    t.addObject("MeshMatrixMass", name="mass", massDensity=rho, topology="@topo")
    cls = {"FTCfix": "FastTetrahedralCorotationalForceFieldFixed",
           "FTC": "FastTetrahedralCorotationalForceField", "TFEM": "TetrahedronFEMForceField"}[ff_type]
    ff = t.addObject(cls, name="fem", youngModulus=E, poissonRatio=nu, method=method, topology="@topo")
    return t, mo, ff


def add_static_solver(node, linear="cg", newton_iters=50, abs_tol=1e-20):
    """静力求解：StaticSolver + NewtonRaphsonSolver + 线性求解器。

    linear: "cg"（不组装矩阵）/ "ldl"（SparseLDLSolver）/ "eigen"（EigenSimplicialLDLT）。
    原版 FTC 的组装矩阵有 bug，只能用 "cg"；修复后的 FTCfix 三种都可以（见 docs/plan/step0_results.md）。
    abs_tol：残差范数（N）的绝对停止阈值；按问题的力的量级设（例如 1e-9 × 典型力）。
    阈值过严（低于舍入误差下限）时 Newton 会跑满 newton_iters 次。
    """
    node.addObject("NewtonRaphsonSolver", name="newton",
                   maxNbIterationsNewton=newton_iters,
                   absoluteResidualStoppingThreshold=abs_tol,
                   relativeInitialStoppingThreshold=1e-14,
                   relativeSuccessiveStoppingThreshold=1e-14,
                   relativeEstimateDifferenceThreshold=1e-14,
                   absoluteEstimateDifferenceThreshold=1e-16,
                   warnWhenLineSearchFails=False, warnWhenDiverge=False)
    node.addObject("StaticSolver", name="static", newtonSolver="@newton")
    if linear == "cg":
        node.addObject("CGLinearSolver", name="linsolver", iterations=5000,
                       tolerance=1e-24, threshold=1e-30)
    elif linear == "eigen":
        node.addObject("EigenSimplicialLDLT", name="linsolver",
                       template="CompressedRowSparseMatrixMat3x3d")
    else:
        node.addObject("SparseLDLSolver", name="linsolver",
                       template="CompressedRowSparseMatrixMat3x3d")


def internal_force(root, node, constraints):
    """读当前构形下的节点内力（只用于静力测试的最后一步，会拆掉约束）。

    投影约束会把受约束自由度上的力清零，所以先移除 constraints，再让 Newton 迭代 0 次：
    只计算一次内力、不移动节点。返回 (N,3) 数组：f = 内力（外载为 0 时，反力 = -f）。
    """
    for c in constraints:
        node.removeObject(c)
    root.newton.maxNbIterationsNewton = 0
    mo = node.dofs
    x0 = mo.position.array().copy()
    Sofa.Simulation.animate(root, root.dt.value)
    if np.all(np.isfinite(x0)):  # 发散（NaN）的对照工况跳过检查
        assert np.abs(mo.position.array() - x0).max() == 0.0
    return mo.force.array().copy()


def make_root(gravity=(0.0, 0.0, 0.0), dt=0.01):
    root = Sofa.Core.Node("root")
    root.gravity = list(gravity)
    root.dt = dt
    root.addObject("DefaultAnimationLoop")
    return root
