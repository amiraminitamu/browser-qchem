#!/usr/bin/env bash
set -euo pipefail
project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
for tool in cmake emcmake emcc; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "Missing $tool. See the README setup instructions." >&2
        exit 1
    fi
done
emcmake cmake -S "$project_dir" -B "$project_dir/build/wasm" \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_POLICY_VERSION_MINIMUM=3.5
cmake --build "$project_dir/build/wasm" --parallel "${CMAKE_BUILD_PARALLEL_LEVEL:-4}"
echo "Ready. Run: python3 -m http.server 8000 --bind 127.0.0.1 --directory build/wasm/dist"
