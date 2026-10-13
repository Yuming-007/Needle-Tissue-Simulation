#!/usr/bin/env bash
# 远程（SSH）时在游戏本屏幕上运行 GUI 场景并截图 / 录屏。流程固定为：解锁 → 运行 + 截图 / 录制 → 立刻锁回（用户 2026-10-13 授权，见 CLAUDE.md）。
# 无论成功、出错还是超时，退出时都会关闭 runSofa 并锁屏；另有看门狗：超过 MAX_SECONDS（默认 60）强制结束。
#
# 用法（项目目录下）：
#   tools/remote_capture.sh shot  <输出前缀> "<秒1> <秒2> ..." <场景> [场景参数...]
#       窗口出现后第 秒1、秒2… 秒各截一张窗口图 → <输出前缀>_<秒>.png
#   tools/remote_capture.sh video <输出.mp4> <最长秒数> <结束标记> <场景> [场景参数...]
#       录 SOFA 窗口（30 帧/秒，h264，原始分辨率，不缩放），日志里出现 <结束标记> 后再录 1.5 s 结束，最多 <最长秒数> 秒
# 场景自动开始（runSofa -a）；窗口大小由场景里的 ViewerSetting 决定。终端输出存在 <输出>.log。
set -u
MODE=$1; OUT=$2; shift 2
export DISPLAY=:0 XAUTHORITY=/run/user/1000/gdm/Xauthority DBUS_SESSION_BUS_ADDRESS=unix:path=/run/user/1000/bus
PROJECT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MAX_SECONDS=${MAX_SECONDS:-60}
SESSION=$(loginctl list-sessions --no-legend | awk '$4=="seat0"{print $1; exit}')

cleanup() {
  pkill -x runSofa 2>/dev/null
  [ -n "${FPID:-}" ] && kill -INT "$FPID" 2>/dev/null && wait "$FPID" 2>/dev/null
  loginctl lock-session "$SESSION"
  [ -n "${WATCHDOG:-}" ] && kill "$WATCHDOG" 2>/dev/null
  sleep 1
  echo "已锁屏：$(gdbus call --session --dest org.gnome.ScreenSaver --object-path /org/gnome/ScreenSaver --method org.gnome.ScreenSaver.GetActive)"
}
trap cleanup EXIT
trap 'exit 1' INT TERM
( sleep "$MAX_SECONDS"; pkill -x runSofa; loginctl lock-session "$SESSION" ) & WATCHDOG=$!

loginctl unlock-session "$SESSION"; xset dpms force on; xset s reset; sleep 1.5

if [ "$MODE" = shot ]; then TIMES=$1; shift; else MAXREC=$1; STOP=$2; shift 2; fi
SCENE=$1; shift
LOG="${OUT%.*}.log"; [ "$MODE" = shot ] && LOG="${OUT}.log"
cd "$PROJECT"
RUNSOFA_OPTS="-g glfw -a" tools/run_gui.sh "$SCENE" "$@" > "$LOG" 2>&1 &

TITLE="SOFA - $SCENE"
W=""
for i in $(seq 1 150); do
  W=$(xwininfo -root -tree | grep -F "\"$TITLE\": (\"$TITLE\"" | head -1 | awk '{print $1}')
  [ -n "$W" ] && break; sleep 0.2
done
[ -z "$W" ] && { echo "没找到 SOFA 窗口"; exit 1; }
sleep 0.5
read -r WW HH X Y <<< "$(xwininfo -id "$W" | awk '/Absolute upper-left X/{x=$4}/Absolute upper-left Y/{y=$4}/Width/{w=$2}/Height/{h=$2}END{print w, h, x, y}')"
WW=$((WW / 2 * 2)); HH=$((HH / 2 * 2))
echo "窗口 $W：${WW}x${HH}+${X}+${Y}"
T0=$(date +%s.%N)

if [ "$MODE" = shot ]; then
  for t in $TIMES; do
    while (( $(echo "$(date +%s.%N) - $T0 < $t" | bc) )); do sleep 0.05; done
    ffmpeg -v error -y -f x11grab -video_size "${WW}x${HH}" -i ":0.0+$X,$Y" -frames:v 1 "${OUT}_${t}.png"
    echo "截图 ${OUT}_${t}.png（$(grep -E '^\[针\]' "$LOG" | tail -1 | awk -F'：' '{print $1}')）"
  done
else
  ffmpeg -y -loglevel error -f x11grab -framerate 30 -video_size "${WW}x${HH}" -i ":0.0+$X,$Y" \
         -c:v libx264 -pix_fmt yuv420p -preset veryfast -t "$MAXREC" "$OUT" < /dev/null &
  FPID=$!
  while kill -0 "$FPID" 2>/dev/null && ! grep -qF "$STOP" "$LOG"; do sleep 0.3; done
  sleep 1.5
  kill -INT "$FPID" 2>/dev/null; wait "$FPID" 2>/dev/null; FPID=""
  echo "录制结束，$(echo "$(date +%s.%N) - $T0" | bc | cut -d. -f1) s → $OUT"
fi
