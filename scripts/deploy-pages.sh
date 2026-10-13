#!/usr/bin/env bash
# Deploy the Quartz site to GitHub Pages by hand: job "deploy" of the former GitHub Actions workflow
# .github/workflows/ci.yaml, removed on 2026-10-13 (LabOnoM does not use GitHub Actions), which ran after
# "build-and-test" on every push to main. When to run: after merging into main, from an up-to-date main.
#
# The workflow published with actions/upload-pages-artifact and actions/deploy-pages, which work only inside GitHub
# Actions. This script publishes the same files as the whole content of the gh-pages branch instead.
# ONE-TIME OWNER STEP (a repository setting; this script does not change it): after the first run, set
# Settings > Pages > Build and deployment > Source to "Deploy from a branch", branch gh-pages, folder "/ (root)".
# Until then https://www.bs-gou.com/DK.BeesGO/ keeps serving the last deployment made by GitHub Actions.
# Credentials: the owner's own git access to github.com (for example gh auth login); nothing is read from a file
# in the repository.
set -euo pipefail
cd "$(dirname "$0")/.."

# Guards: the workflow deployed what was pushed to main, so deploy only a clean, pushed main.
# (File-mode changes are ignored: npx marks source/quartz/bootstrap-cli.mjs executable on every build.)
branch=$(git rev-parse --abbrev-ref HEAD)
[ "$branch" = main ] || { echo "FAIL  deploy from main, not $branch" >&2; exit 1; }
git fetch --quiet origin main
[ "$(git rev-parse HEAD)" = "$(git rev-parse origin/main)" ] || { echo "FAIL  local main differs from origin/main" >&2; exit 1; }
git -c core.fileMode=false diff --quiet HEAD -- || { echo "FAIL  tracked files have uncommitted changes" >&2; exit 1; }
sha=$(git rev-parse HEAD)
remote=$(git remote get-url origin)

# needs: build-and-test
bash scripts/ci-local.sh

# "Build Quartz"
cd source
npm install
npx quartz build

# "Copy raw HTML folder"
cp -r raw_html/* public/ || true

# "Upload GitHub Pages artifact" and "Deploy to GitHub Pages": source/public becomes the whole gh-pages branch.
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
cp -r public/. "$tmp/"
touch "$tmp/.nojekyll"   # serve the files as built; no Jekyll processing on a branch build
git -C "$tmp" init --quiet -b gh-pages
git -C "$tmp" add -A
git -C "$tmp" commit --quiet -m "Deploy $sha"
git -C "$tmp" push --quiet --force "$remote" gh-pages
echo "OK    pushed the site built from $sha to gh-pages"
echo "      Pages source must be 'Deploy from a branch: gh-pages' for https://www.bs-gou.com/DK.BeesGO/ to serve it."
