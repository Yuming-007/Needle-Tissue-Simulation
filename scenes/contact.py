"""针–组织接触的场景构建（基于 InfinyTech3D CollisionAlgorithm + ConstraintGeometry v25.12）。

组件组合参考插件官方场景 `CollisionAlgorithm/scenes/NeedleInsertion.py:136-205`（只读，`docs/code/01`），
差别（第 1 步）：
  - 针是纯运动学的（common.add_kinematic_needle + NeedleDriver），不是梁 + 弹簧；
  - 刺穿阈值不设（默认无穷大，InsertionAlgorithm.cpp:59-65），关闭针身碰撞和插入，只保留针尖–表面的单边约束；
  - 约束方向用 SecondDirection（表面法向取反），和官方场景一致（不用 ContactDirection，见 step1_tenting.md §1.4）。
"""
import sys
sys.dont_write_bytecode = True
import os
import numpy as np
import SofaRuntime
import common

PLUGIN_DIR = os.path.expanduser("~/sofa/reference_install/lib")
_loaded = False


def load_plugins():
    global _loaded
    if not _loaded:
        SofaRuntime.PluginRepository.addFirstPath(PLUGIN_DIR)
        SofaRuntime.importPlugin("CollisionAlgorithm")
        SofaRuntime.importPlugin("ConstraintGeometry")
        _loaded = True


def add_tissue_with_surface(root, n=9, L=0.08, E=5000.0, nu=0.45, rho=1000.0, linear="ldl", swapping=True,
                            flip_normals=False):
    """组织（动力学，FTCfix polar）+ 底面固定 + 约束修正 + 表面几何（插件）。

    swapping=True（默认）：Hexa2Tetra 相邻六面体的切分方向交替，网格在中心节点处镜像对称。
    swapping=False 时，中心节点受竖直点载荷会沿对角线侧向漂移约 30%（第 1 步发现），接触点会离开节点。

    返回 dict：node、dofs、surface 节点、表面几何 / 法向组件、体几何、底面节点编号、网格节点坐标 X。
    """
    tr = root.addChild("TissueRoot")
    common.add_dynamic_solver(tr, linear=linear)
    t, mo, ff = common.add_tissue(tr, L=L, n=n, E=E, nu=nu, rho=rho, swapping=swapping)
    X = common.grid_nodes(L, n)
    bot = np.where(np.isclose(X[:, 2], 0))[0]
    t.addObject("FixedProjectiveConstraint", name="bottom", indices=bot.tolist())
    t.addObject("LinearSolverConstraintCorrection", linearSolver="@../linsolver")
    vol = t.addObject("TetrahedronGeometry", name="geom_tetra", mstate="@dofs", topology="@topo", draw=False)
    s = t.addChild("Surface")
    s.addObject("TriangleSetTopologyContainer", name="tris")
    s.addObject("TriangleSetTopologyModifier")
    # flip_normals=False（插件 A1 用）：顶面法向指向组织内部（第 1 步实测），配合 SecondDirection（取反）使用。
    # flip_normals=True（SOFA 自带接触 A2 用）：法向朝外，点–三角形接触按正面判断。
    s.addObject("Tetra2TriangleTopologicalMapping", input="@../topo", output="@tris", flipNormals=flip_normals)
    smo = s.addObject("MechanicalObject", name="dofs", template="Vec3d", position="@../topo.position")
    sg = s.addObject("TriangleGeometry", name="geom_tri", mstate="@dofs", topology="@tris", draw=False)
    nh = s.addObject("PhongTriangleNormalHandler", name="normals", geometry="@geom_tri")
    s.addObject("AABBBroadPhase", name="aabb", thread=1, nbox=[2, 2, 2], method=2)
    s.addObject("IdentityMapping", input="@../dofs", output="@dofs", isMechanical=True)
    return dict(node=t, dofs=mo, fem=ff, surface=s, surf_dofs=smo, surf_geom=sg, normals=nh, vol_geom=vol,
                bottom=bot, X=X, L=L, n=n)


def add_needle_geometry(nd, body, tip, n_body):
    """给纯运动学针的子节点加上插件的几何：针尖 PointGeometry，针身 EdgeGeometry + EdgeNormalHandler。"""
    tip_node = tip.getContext()
    tg = tip_node.addObject("PointGeometry", name="geom_tip", mstate="@dofs")
    body_node = body.getContext()
    body_node.addObject("EdgeSetTopologyContainer", name="edges", position="@dofs.position",
                        edges=[[i, i + 1] for i in range(n_body - 1)])
    bg = body_node.addObject("EdgeGeometry", name="geom_body", mstate="@dofs", topology="@edges")
    body_node.addObject("EdgeNormalHandler", name="needleEdges", geometry="@geom_body")
    return tg, bg


def add_tip_contact(root, tip_geom, body_geom, tissue, distance=0.01):
    """针尖–表面的单边接触（刺穿前）。第 1 步：不刺穿、不插入、不检测针身。"""
    algo = root.addObject("InsertionAlgorithm", name="algo",
                          tipGeom=tip_geom.getLinkPath(), surfGeom=tissue["surf_geom"].getLinkPath(),
                          shaftGeom=body_geom.getLinkPath(), volGeom=tissue["vol_geom"].getLinkPath(),
                          enableShaftCollision=False, enableInsertion=False)
    root.addObject("DistanceFilter", algo="@algo", distance=distance)
    root.addObject("SecondDirection", name="contactDir", handler=tissue["normals"].getLinkPath())
    cu = root.addObject("ConstraintUnilateral", name="tipContact", input="@algo.collisionOutput",
                        directions="@contactDir")
    return algo, cu
