import sys, os, json, importlib.util, time, numpy as np
sys.dont_write_bytecode = True
import Sofa, SofaRuntime
SCENE_DIR=os.path.expanduser("~/sofa/reference_sources/CollisionAlgorithm/scenes")
SofaRuntime.PluginRepository.addFirstPath(os.path.expanduser("~/sofa/reference_install/lib"))
spec=importlib.util.spec_from_file_location("s",os.path.join(SCENE_DIR,"NeedleInsertion.py")); s=importlib.util.module_from_spec(spec); spec.loader.exec_module(s)
root=Sofa.Core.Node("root"); s.createScene(root)
solver=[o for o in root.objects if o.getClassName()=="BlockGaussSeidelConstraintSolver"][0]
ci=[o for o in root.objects if o.getClassName()=="ConstraintInsertion"][0]
here=os.path.dirname(os.path.abspath(__file__)); os.chdir(SCENE_DIR); Sofa.Simulation.initRoot(root)
for _ in range(80): Sofa.Simulation.animate(root, root.dt.value)   # needle inserted, coupling points exist
def run(n):
    for _ in range(n): Sofa.Simulation.animate(root, root.dt.value)
    return dict(t=round(root.time.value,3), nC=int(solver.currentNumConstraints.value),
                gel=float(np.abs(np.array(root.Volume.mstate_gel.position.value)).sum()))
out=os.path.join(here,"fork_result_%d.json")
t0=time.time(); pids=[]
for k,mu in enumerate([0.002, 0.002, 0.2]):          # branch 0 and 1 identical, branch 2 different
    pid=os.fork()
    if pid==0:
        ci.frictionCoeff.value=mu; r=run(20); r["mu"]=mu
        json.dump(r,open(out%k,"w")); os._exit(0)
    pids.append(pid)
for p in pids: os.waitpid(p,0)
tfork=time.time()-t0
parent=run(20)   # parent continues from the snapshot, unaffected by children
print(json.dumps(dict(parent=parent, children=[json.load(open(out%k)) for k in range(3)], wall_s=round(tfork,2))))
