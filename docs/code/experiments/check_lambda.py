import sys, os, importlib.util, numpy as np
sys.dont_write_bytecode = True  # never write __pycache__ into the read-only scene directory
import Sofa, SofaRuntime
SCENE_DIR=os.path.expanduser("~/sofa/reference_sources/CollisionAlgorithm/scenes")
SofaRuntime.PluginRepository.addFirstPath(os.path.expanduser("~/sofa/reference_install/lib"))
spec=importlib.util.spec_from_file_location("s",os.path.join(SCENE_DIR,"NeedleInsertion.py")); s=importlib.util.module_from_spec(spec); spec.loader.exec_module(s)
dt=float(sys.argv[1])
root=Sofa.Core.Node("root"); s.createScene(root); root.dt.value=dt
solver=[o for o in root.objects if o.getClassName()=="BlockGaussSeidelConstraintSolver"][0]
solver.computeConstraintForces.value=True
class P(Sofa.Core.Controller):
    def onAnimateBeginEvent(self,e):
        # stored lambda (motion space) on tip mstate, read at the same moment as the puncture check
        lam=np.array(root.Needle.tipCollision.mstate_tip.findData('lambda').value) if root.Needle.tipCollision.mstate_tip.findData('lambda') else None
        self.lam=lam
    def onAnimateEndEvent(self,e):
        t=root.time.value
        if 0.66<=t<=0.73:
            f=np.array(solver.constraintForces.value)
            base=root.Needle.needleBase.mstate_base.position.value[0]; m=root.NeedleBaseMaster.mstate_baseMaster.position.value[0]
            print(f"t={t:.2f} spring_Fz={1e8*(m[2]-base[2]):8.3f}  GS f={np.round(f,4)}")
root.addObject(P(name="p"))
os.chdir(SCENE_DIR); Sofa.Simulation.initRoot(root)
for _ in range(int(0.73/dt)+1): Sofa.Simulation.animate(root, dt)
