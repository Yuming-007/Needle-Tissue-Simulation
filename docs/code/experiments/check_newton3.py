import sys, os, importlib.util, numpy as np
sys.dont_write_bytecode = True  # never write __pycache__ into the read-only scene directory
import Sofa, SofaRuntime
SCENE_DIR=os.path.expanduser("~/sofa/reference_sources/CollisionAlgorithm/scenes")
SofaRuntime.PluginRepository.addFirstPath(os.path.expanduser("~/sofa/reference_install/lib"))
spec=importlib.util.spec_from_file_location("s",os.path.join(SCENE_DIR,"NeedleInsertion.py")); s=importlib.util.module_from_spec(spec); spec.loader.exec_module(s)
dt=float(sys.argv[1]); needle_first_order = sys.argv[2]=="1"
root=Sofa.Core.Node("root"); s.createScene(root); root.dt.value=dt
[o for o in root.Needle.objects if o.getClassName()=="EulerImplicitSolver"][0].firstOrder.value=needle_first_order
solver=[o for o in root.objects if o.getClassName()=="BlockGaussSeidelConstraintSolver"][0]; solver.computeConstraintForces.value=True
vol=root.Volume; box=[o for o in vol.objects if o.getClassName()=="BoxROI"][0]
class P(Sofa.Core.Controller):
    def onAnimateEndEvent(self,e):
        t=root.time.value
        if 0.60<=t<=0.715:
            x=np.array(vol.mstate_gel.position.value); x0=np.array(vol.mstate_gel.rest_position.value)
            idx=np.array(box.indices.value,dtype=int)
            Fb=1e6*(x0[idx]-x[idx]).sum(axis=0)
            base=root.Needle.needleBase.mstate_base.position.value[0]; m=root.NeedleBaseMaster.mstate_baseMaster.position.value[0]
            tip=root.Needle.tipCollision.mstate_tip.position.value[0]
            f=solver.constraintForces.value
            print(f"t={t:.3f} needleSpringFz={1e8*(m[2]-base[2]):9.3f} GSf={(f[0] if len(f) else 0):9.3f} gelBottomFz={Fb[2]:10.3f} indent={(-0.1-tip[2])*1000:6.2f}mm")
root.addObject(P(name="p"))
os.chdir(SCENE_DIR); Sofa.Simulation.initRoot(root)
for _ in range(int(0.715/dt)+1): Sofa.Simulation.animate(root, dt)
