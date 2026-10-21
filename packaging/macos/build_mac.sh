#!/usr/bin/env bash
# Builds "Mol3D Bridge.app" and a .dmg on macOS (Apple Silicon or Intel). Needs Python 3.10-3.13.
set -euo pipefail
cd "$(dirname "$0")/../.."
PY=""
for v in 3.12 3.13 3.11 3.10 3.14 3.9; do command -v python$v >/dev/null && { PY=python$v; break; }; done
PY=${PY:-python3}
$PY -c "import sys; assert sys.version_info >= (3,9)" || { echo "Python 3.9+ needed"; exit 1; }
[ -d .buildenv ] && .buildenv/bin/python -c "import sys; assert sys.version_info >= (3,9)" 2>/dev/null || rm -rf .buildenv
[ -d .buildenv ] || $PY -m venv .buildenv
.buildenv/bin/python -m pip install --upgrade pip
.buildenv/bin/python -m pip install ".[build]"
.buildenv/bin/python packaging/build.py
echo "DMG written to dist/"
