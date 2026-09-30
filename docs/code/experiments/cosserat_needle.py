import sys, os
sys.dont_write_bytecode = True
EX="/tmp/claude-1000/-home-yuming-sofa-needle-project/982afcea-fafe-44ba-8bab-2e129c167ce6/scratchpad/cosserat_ex"
sys.path.insert(0, os.path.expanduser("~/sofa/resources/install_cosserat/lib/python3/site-packages")); sys.path.insert(0, EX)
import numpy as np, Sofa, SofaRuntime
SofaRuntime.PluginRepository.addFirstPath(os.path.expanduser("~/sofa/resources/install_cosserat/lib"))
import importlib.util
spec=importlib.util.spec_from_file_location("ns",os.path.join(EX,"NeedleInsertion.py")); ns=importlib.util.module_from_spec(spec); spec.loader.exec_module(ns)
SofaRuntime.importPlugin("Sofa.Component"); SofaRuntime.importPlugin("Sofa.GL.Component")
root=Sofa.Core.Node("root"); ns.createScene(root); Sofa.Simulation.initRoot(root)
anim=[o for o in root.objects if o.getClassName()=="Animation" or type(o).__name__=="Animation"][0]
rb=anim.rigidBaseMO; step_cm=float(sys.argv[1]) if len(sys.argv)>1 else 0.05
for k in range(int(sys.argv[2]) if len(sys.argv)>2 else 300):
    with rb.rest_position.writeable() as p: p[0][0]+=step_cm
    Sofa.Simulation.animate(root, root.dt.value)
    if k%20==0 or anim.inside and k%10==0:
        n=len(anim.constraintPts.position.value) if anim.inside else 0
        f=anim.generic.constraintForces.value
        print(f"k={k:3d} base_x={rb.rest_position.value[0][0]:6.2f}cm inside={anim.inside} nCP={n} |f|max={(np.abs(f).max() if len(f) else 0):.3f}")
