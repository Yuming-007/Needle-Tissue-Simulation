import sys, os, json, importlib.util, numpy as np
sys.dont_write_bytecode = True  # never write __pycache__ into the read-only scene directory
import Sofa, SofaRuntime
SCENE_DIR=os.path.expanduser("~/sofa/reference_sources/CollisionAlgorithm/scenes")
SofaRuntime.PluginRepository.addFirstPath(os.path.expanduser("~/sofa/reference_install/lib"))
scene_file, nsteps, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
spec=importlib.util.spec_from_file_location("s",os.path.join(SCENE_DIR,scene_file)); s=importlib.util.module_from_spec(spec); spec.loader.exec_module(s)
root=Sofa.Core.Node("root"); s.createScene(root)
solver=[o for o in root.objects if "ConstraintSolver" in o.getClassName()][0]
rows=[]
class P(Sofa.Core.Controller):
    def onAnimateEndEvent(self,e):
        g=np.array(root.Volume.mstate_gel.position.value); g0=np.array(root.Volume.mstate_gel.rest_position.value); d=np.linalg.norm(g-g0,axis=1)
        rows.append(dict(t=float(root.time.value), nC=int(solver.currentNumConstraints.value), it=int(solver.currentIterations.value), gelMaxDisp=float(d.max())))
root.addObject(P(name="p")); os.chdir(SCENE_DIR); Sofa.Simulation.initRoot(root)
for _ in range(nsteps): Sofa.Simulation.animate(root, root.dt.value)
os.chdir(os.path.dirname(os.path.abspath(__file__))); json.dump(rows, open(out,"w")); print("done", len(rows))
