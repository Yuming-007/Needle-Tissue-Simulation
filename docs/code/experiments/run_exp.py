"""Headless experiment on the official v25.12 NeedleInsertion.py scene (read-only use).

Usage: python3 run_exp.py <frictionCoeff> <tolerance> <maxIt> <nsteps> <out.json> [insertion_only]
Records per step: GS iterations/error, #constraints, needle base axial force (from the
RestShapeSprings coupling), and determinant / diagonal statistics of the 3x3 diagonal
blocks of the constraint compliance matrix W.
"""
import sys, os, json, importlib.util
sys.dont_write_bytecode = True  # never write __pycache__ into the read-only scene directory
import numpy as np
import Sofa, SofaRuntime
from Sofa import SofaConstraintSolver  # noqa: F401  (enables solver.W())

SCENE_DIR = os.path.expanduser("~/sofa/reference_sources/CollisionAlgorithm/scenes")
PLUGIN_DIR = os.path.expanduser("~/sofa/reference_install/lib")

mu, tol, maxit, nsteps, out = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
DT = float(sys.argv[6]) if len(sys.argv) > 6 else None

SofaRuntime.PluginRepository.addFirstPath(PLUGIN_DIR)
spec = importlib.util.spec_from_file_location("needle_scene", os.path.join(SCENE_DIR, "NeedleInsertion.py"))
scene = importlib.util.module_from_spec(spec); spec.loader.exec_module(scene)


def find(node, cls):
    for o in node.objects:
        if o.getClassName() == cls:
            return o
    return None


class Probe(Sofa.Core.Controller):
    def __init__(self, *a, **k):
        Sofa.Core.Controller.__init__(self, *a, **k)
        self.solver = k["solver"]; self.root = k["root"]
        self.rows = []; self.cur = {}

    def onBuildConstraintSystemEndEvent(self, e):
        W = np.array(self.solver.W())
        n = W.shape[0] if W.ndim == 2 else 0
        dets, diags = [], []
        for i in range(0, n - n % 3, 3):
            B = W[i:i + 3, i:i + 3]
            dets.append(abs(np.linalg.det(B))); diags.extend(np.diag(B).tolist())
        self.cur = dict(n=int(n),
                        minDet=float(min(dets)) if dets else None,
                        maxDet=float(max(dets)) if dets else None,
                        nDetBelowEps=int(sum(d <= np.finfo(float).eps for d in dets)),
                        minDiag=float(min(diags)) if diags else None,
                        maxDiag=float(max(diags)) if diags else None)

    def onAnimateEndEvent(self, e):
        base = self.root.Needle.needleBase.mstate_base.position.value[0]
        master = self.root.NeedleBaseMaster.mstate_baseMaster.position.value[0]
        f = 1e8 * (np.array(master[:3]) - np.array(base[:3]))  # spring force on needle base
        tip = self.root.Needle.tipCollision.mstate_tip.position.value[0]
        r = dict(t=round(float(self.root.time.value), 4),
                 it=int(self.solver.currentIterations.value), err=float(self.solver.currentError.value),
                 nC=int(self.solver.currentNumConstraints.value),
                 Fz=float(f[2]), Fx=float(f[0]), tipZ=float(tip[2]))
        r.update(self.cur); self.cur = {}
        self.rows.append(r)


root = Sofa.Core.Node("root")
scene.createScene(root)
if DT: root.dt.value = DT
solver = find(root, "BlockGaussSeidelConstraintSolver")
solver.tolerance.value = tol; solver.maxIterations.value = maxit
solver.computeConstraintForces.value = True
ci = find(root, "ConstraintInsertion"); ci.frictionCoeff.value = mu
probe = root.addObject(Probe(name="probe", solver=solver, root=root))

old = os.getcwd(); os.chdir(SCENE_DIR)  # ReadState uses a relative path
Sofa.Simulation.initRoot(root)
for _ in range(nsteps):
    Sofa.Simulation.animate(root, root.dt.value)
os.chdir(old)
json.dump(dict(mu=mu, tol=tol, maxit=maxit, dt=root.dt.value, rows=probe.rows), open(out, "w"), indent=0)
print("done", out, len(probe.rows))
