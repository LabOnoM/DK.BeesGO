#!/usr/bin/env bash
# Local check: job "build-and-test" of the former GitHub Actions workflow .github/workflows/ci.yaml, removed on
# 2026-10-13 (LabOnoM does not use GitHub Actions). Run before every merge: bash scripts/ci-local.sh
# Needs Node 20 or >= 22 (source/package.json "engines"; the workflow used Node 20) and network access for npm.
# The build writes source/public/ and source/quartz/.quartz-cache/ (not committed). npx marks
# source/quartz/bootstrap-cli.mjs executable, which git shows as a mode change; it is harmless.
set -euo pipefail
cd "$(dirname "$0")/../source"

done_steps=()
step() { printf '\n== %s\n' "$1"; current="$1"; }
pass() { done_steps+=("$current"); }
summary() {
  status=$?
  printf '\n== Summary\n'
  for s in "${done_steps[@]}"; do echo "PASS  $s"; done
  if [ "$status" -ne 0 ]; then echo "FAIL  ${current:-setup} (exit $status)"; else echo "All ${#done_steps[@]} steps passed."; fi
}
trap summary EXIT

# --- Job "build-and-test" ---
step "build-and-test: npm install (Node $(node --version))"
npm install
pass

step "build-and-test: chmod +x ./quartz/bootstrap-cli.mjs"
chmod +x ./quartz/bootstrap-cli.mjs
pass

step "build-and-test: npx quartz build --bundleInfo"
npx quartz build --bundleInfo
pass
