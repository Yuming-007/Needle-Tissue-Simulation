#!/usr/bin/env bash
# 用 runSofa 图形界面打开一个本项目的 Python 场景。
# 用法（在任意目录）：tools/run_gui.sh scenes/step0/t0d_gui_scene.py [场景参数...]
# 无界面测试：RUNSOFA_OPTS="-g batch -n 300" tools/run_gui.sh <场景>
# 换成简单的 GLFW 界面（没有 ImGui 面板）：RUNSOFA_OPTS="-g glfw" tools/run_gui.sh <场景>
# 用 NVIDIA 独显绘制（双显卡笔记本）：RUNSOFA_GPU=nvidia tools/run_gui.sh <场景>
set -e
PROJECT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export SOFA_ROOT="$HOME/sofa/SOFA_v25.12.00_Linux"
export PYTHONPATH="$SOFA_ROOT/plugins/SofaPython3/lib/python3/site-packages:$SOFA_ROOT/plugins/STLIB/lib/python3/site-packages"
export LD_LIBRARY_PATH="$SOFA_ROOT/lib"
export PYTHONDONTWRITEBYTECODE=1
if [ "${RUNSOFA_GPU:-}" = "nvidia" ]; then
  export __NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia __GL_SYNC_TO_VBLANK=0
fi
SCENE="${1:?请给出场景文件，例如 scenes/step0/t0d_gui_scene.py}"; shift
cd "$PROJECT"
if [ $# -gt 0 ]; then
  exec "$SOFA_ROOT/bin/runSofa" ${RUNSOFA_OPTS} -l SofaPython3 "$SCENE" --argv "$@"
else
  exec "$SOFA_ROOT/bin/runSofa" ${RUNSOFA_OPTS} -l SofaPython3 "$SCENE"
fi
