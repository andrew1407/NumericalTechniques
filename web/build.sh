#!/usr/bin/env bash
# Compile every lab in sripts/ unchanged to an ES6 WebAssembly module in web/wasm/.
# Needs emscripten (source ~/emsdk/emsdk_env.sh). -O0 matches the labs' own `c++ file.cpp`
# build and keeps their undefined behaviour (uninitialised accumulators) as it is natively.
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p web/wasm
for src in sripts/*-lab.cpp; do
  lab=$(basename "$src" -lab.cpp)
  exports=_main
  case $lab in
    2) exports=_main,_M ;;       # float M[4][5]: the system's coefficients
    5) exports=_main,_X,_Y ;;    # float X[5], Y[5]: the interpolation table
  esac
  em++ -O0 -w "$src" -o "web/wasm/lab$lab.mjs" \
    -sMODULARIZE -sEXPORT_ES6 -sEXPORT_NAME=createLab -sINVOKE_RUN=0 \
    -sENVIRONMENT=web,worker -sEXPORTED_FUNCTIONS=$exports \
    -sEXPORTED_RUNTIME_METHODS=callMain,HEAPF32
  echo "built web/wasm/lab$lab.mjs"
done
