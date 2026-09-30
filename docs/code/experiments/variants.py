"""Variants on the official scene: needle integrator order and constraint-correction type.
usage: variants.py <gridN> <needleOrder: 1|2> <cc: lin|pre> <out.json>"""
import sys, os, json, time, importlib.util, numpy as np
sys.dont_write_bytecode = True
import Sofa, SofaRuntime
SCENE_DIR=os.path.expanduser("~/sofa/reference_sources/CollisionAlgorithm/scenes")
SofaRuntime.PluginRepository.addFirstPath(os.path.expanduser("~/sofa/reference_install/lib"))
N, order, cc, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4]
spec=importlib.util.spec_from_file_location("s",os.path.join(SCENE_DIR,"NeedleInsertion.py")); s=importlib.util.module_from_spec(spec); spec.loader.exec_module(s)
s.g_gelRegularGridParameters["n"]=[N,N,N]
root=Sofa.Core.Node("root"); s.createScene(root)
f=lambda node,cls:[o for o in node.objects if o.getClassName()==cls][0]
f(root.Needle,"EulerImplicitSolver").firstOrder.value = (order==1)
tipd = 0.25/(N-1)                       # constraint spacing = element size (avoid over-constraint)
f(root,"InsertionAlgorithm").tipDistThreshold.value = tipd
solver=f(root,"BlockGaussSeidelConstraintSolver"); solver.tolerance.value=1e-8; solver.maxIterations.value=2000
f(root,"ConstraintInsertion").frictionCoeff.value = 0.0
if cc=="pre":
    vol=root.Volume; old=f(vol,"LinearSolverConstraintCorrection"); vol.removeObject(old)
    vol.addObject("PrecomputedConstraintCorrection", rotations=True, recompute=True)
vol=root.Volume; box=[o for o in vol.objects if o.getClassName()=="BoxROI"][0]
rows=[]
class P(Sofa.Core.Controller):
    def onAnimateBeginEvent(self,e): self.t0=time.perf_counter()
    def onAnimateEndEvent(self,e):
        x=np.array(vol.mstate_gel.position.value); x0=np.array(vol.mstate_gel.rest_position.value); idx=np.array(box.indices.value,dtype=int)
        base=root.Needle.needleBase.mstate_base.position.value[0]; m=root.NeedleBaseMaster.mstate_baseMaster.position.value[0]
        nd=np.array(root.Needle.mstate.position.value)[:,:3]
        rows.append(dict(t=round(root.time.value,3), ms=(time.perf_counter()-self.t0)*1e3, it=int(solver.currentIterations.value),
            nC=int(solver.currentNumConstraints.value), needleFz=float(1e8*(m[2]-base[2])),
            tissueFz=float(1e6*(x0[idx]-x[idx]).sum(axis=0)[2]), tipZ=float(root.Needle.tipCollision.mstate_tip.position.value[0][2]), tipX=float(root.Needle.tipCollision.mstate_tip.position.value[0][0]), tipY=float(root.Needle.tipCollision.mstate_tip.position.value[0][1]),
            needleFinite=bool(np.isfinite(nd).all()), needleLen=float(np.linalg.norm(np.diff(nd,axis=0),axis=1).sum())))
root.addObject(P(name="p"))
here=os.path.dirname(os.path.abspath(__file__)); os.chdir(SCENE_DIR)
t0=time.perf_counter(); Sofa.Simulation.initRoot(root); tinit=time.perf_counter()-t0
for _ in range(100): Sofa.Simulation.animate(root, root.dt.value)
os.chdir(here); json.dump(dict(N=N,order=order,cc=cc,tinit_s=tinit,rows=rows),open(out,"w")); print("done",out)
