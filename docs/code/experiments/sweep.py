"""Sweeps on the official v25.12 scene: constraint spacing and mesh resolution.
usage: sweep.py <gridN> <tipDist_m> <mu> <tol> <maxit> <out.json>"""
import sys, os, json, time, importlib.util, numpy as np
sys.dont_write_bytecode = True
import Sofa, SofaRuntime
SCENE_DIR=os.path.expanduser("~/sofa/reference_sources/CollisionAlgorithm/scenes")
SofaRuntime.PluginRepository.addFirstPath(os.path.expanduser("~/sofa/reference_install/lib"))
N, tipd, mu, tol, maxit, out = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
spec=importlib.util.spec_from_file_location("s",os.path.join(SCENE_DIR,"NeedleInsertion.py")); s=importlib.util.module_from_spec(spec); spec.loader.exec_module(s)
s.g_gelRegularGridParameters["n"]=[N,N,N]
root=Sofa.Core.Node("root"); s.createScene(root)
f=lambda cls:[o for o in root.objects if o.getClassName()==cls][0]
solver=f("BlockGaussSeidelConstraintSolver"); solver.tolerance.value=tol; solver.maxIterations.value=maxit
f("ConstraintInsertion").frictionCoeff.value=mu; f("InsertionAlgorithm").tipDistThreshold.value=tipd
vol=root.Volume; box=[o for o in vol.objects if o.getClassName()=="BoxROI"][0]
rows=[]
class P(Sofa.Core.Controller):
    def onAnimateBeginEvent(self,e): self.t0=time.perf_counter()
    def onAnimateEndEvent(self,e):
        x=np.array(vol.mstate_gel.position.value); x0=np.array(vol.mstate_gel.rest_position.value); idx=np.array(box.indices.value,dtype=int)
        base=root.Needle.needleBase.mstate_base.position.value[0]; m=root.NeedleBaseMaster.mstate_baseMaster.position.value[0]
        rows.append(dict(t=round(root.time.value,3), ms=(time.perf_counter()-self.t0)*1e3, it=int(solver.currentIterations.value),
            err=float(solver.currentError.value), nC=int(solver.currentNumConstraints.value),
            needleFz=float(1e8*(m[2]-base[2])), tissueFz=float(1e6*(x0[idx]-x[idx]).sum(axis=0)[2]),
            tipZ=float(root.Needle.tipCollision.mstate_tip.position.value[0][2])))
root.addObject(P(name="p"))
here=os.path.dirname(os.path.abspath(__file__)); os.chdir(SCENE_DIR); Sofa.Simulation.initRoot(root)
ntet=len(vol.TetraContainer.tetrahedra.value)
for _ in range(100): Sofa.Simulation.animate(root, root.dt.value)
os.chdir(here); json.dump(dict(N=N,tipd=tipd,mu=mu,tol=tol,ntet=ntet,nnodes=len(vol.mstate_gel.position.value),rows=rows),open(out,"w"))
print("done",out)
