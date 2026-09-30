import sys, os, importlib.util
sys.dont_write_bytecode = True
import Sofa, SofaRuntime
SCENE_DIR=os.path.expanduser("~/sofa/reference_sources/CollisionAlgorithm/scenes")
SofaRuntime.PluginRepository.addFirstPath(os.path.expanduser("~/sofa/reference_install/lib"))
spec=importlib.util.spec_from_file_location("s",os.path.join(SCENE_DIR,"NeedleInsertion.py")); s=importlib.util.module_from_spec(spec); spec.loader.exec_module(s)
root=Sofa.Core.Node("root"); s.createScene(root); os.chdir(SCENE_DIR); Sofa.Simulation.initRoot(root)
for k in range(100):
    Sofa.Simulation.animate(root, root.dt.value)
    if k in (50,70,99): print(f"t={root.time.value:.2f} Volume (K,U)={root.Volume.computeEnergy()}  Needle (K,U)={root.Needle.computeEnergy()}")
