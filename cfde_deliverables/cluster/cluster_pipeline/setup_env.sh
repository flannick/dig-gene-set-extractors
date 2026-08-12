#!/usr/bin/env bash
# Cluster env for the scRNA->cNMF->programs pipeline. Python 3.10 user-space env.
# Uses MICROMAMBA (tiny, low-memory C++ solver, no root) to AVOID the classic-conda repodata OOM
# (the original 'Killed') AND the missing-libmamba problem. The env is created in your user conda
# envs dir ($HOME/.conda/envs) so the downstream run scripts can still `conda activate gsx310`.
# Compiled stack comes from conda-forge (PREBUILT) -> no source builds against the node's old gcc.
set -uo pipefail
ENV=${ENV:-gsx310}
TOOLCHAIN=${TOOLCHAIN:-/humgen/diabetes2/users/ryank/software/geneset-extractor-dev}
ENVDIR="$HOME/.conda/envs/$ENV"

echo "=== 1. conda (for base path + later activation) ==="
command -v conda >/dev/null 2>&1 || { echo "No conda on PATH (need it for activation). Load anaconda, then re-run."; exit 1; }
source "$(conda info --base)/etc/profile.d/conda.sh"

echo "=== 2. bootstrap micromamba (user-space; no root) ==="
MM=""
command -v micromamba >/dev/null 2>&1 && MM=micromamba
[ -z "$MM" ] && [ -x "$HOME/bin/micromamba" ] && MM="$HOME/bin/micromamba"
if [ -z "$MM" ]; then
  echo "  downloading micromamba ..."
  ( cd "$HOME" && curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj bin/micromamba ) \
    || { echo "micromamba download failed (no internet on this node? run on login node dig-ae-dev-03)"; exit 1; }
  MM="$HOME/bin/micromamba"
fi
echo "  micromamba: $MM"

echo "=== 3. create py3.10 env + CORE numeric stack only (conda-forge ONLY; scanpy via pip later) ==="
# --override-channels -c conda-forge => ignore defaults/pkgs-main/pkgs-r repodata (big memory saver).
# Keep the conda solve SMALL: just python + core numeric. scanpy/anndata/cnmf go in via pip (step 5)
# because their dep tree (numba, igraph, leidenalg, ...) is what blows up the solver's memory.
if [ -x "$ENVDIR/bin/python" ]; then
  echo "  (env exists at $ENVDIR — skipping create)"
else
  # MINIMAL solve: just Python 3.10 (near-zero solver memory). The whole numeric/scanpy stack
  # is installed via pip WHEELS in step 5 (wheels work on 3.10 -> no conda solve, no source build).
  "$MM" create -y -p "$ENVDIR" --override-channels -c conda-forge python=3.10 \
    || { echo ""; echo "micromamba create FAILED. If 'Killed' => this session's memory is too small."; \
         echo "  NOTE: your 'ish -l h_vmem=16' has NO UNIT — request '-l h_vmem=32g' (with the g),"; \
         echo "        OR just run on the LOGIN NODE (no cap):  exit; cd /humgen/diabetes2/users/gage/CFDE; bash setup_env.sh"; \
         echo "        (env is on shared FS, so compute-node jobs use it afterward.)"; exit 1; }
fi

echo "=== 4. activate + HARD-STOP unless python>=3.10 ==="
conda activate "$ENV" 2>/dev/null || conda activate "$ENVDIR" || { echo "activate failed for $ENV / $ENVDIR"; exit 1; }
python - <<'PY' || { echo "ABORT: env is not Python >=3.10"; exit 1; }
import sys
assert sys.version_info[:2] >= (3,10), "need >=3.10"
print("python", sys.version.split()[0])
PY

echo "=== 5a. conda-forge: numeric + scanpy stack (PREBUILT; h5py/numba bundled -> no source builds) ==="
# h5py/numba/scanpy must come from conda-forge: pip would source-build h5py (needs HDF5 dev libs absent here).
"$MM" install -y -p "$ENVDIR" --override-channels -c conda-forge \
  "numpy>=1.23,<2" "pandas>=1.5" "scipy>=1.10" "scikit-learn>=1.2" cython h5py numba scanpy anndata \
  || echo "WARN: conda-forge install failed — send the error"

echo "=== 5b. pip: cnmf (pip-only wheel; deps already satisfied by conda) ==="
python -m pip install -U pip
python -m pip install "cnmf" || echo "WARN: cnmf install failed"
# NOTE: the pipeline is SELF-CONTAINED (scrna_cnmf_programs.py uses cnmf + anndata directly).
# There is NO external 'geneset-extractors' toolchain to install.

echo "=== 6. VERIFY (send me this block) ==="
python -c "import sys; print('python', sys.version.split()[0])"
python - <<'PY' 2>&1 || echo "VERIFY-FAIL: a dep is missing"
import importlib
for m in ("cnmf","scanpy","anndata","sklearn","numpy","pandas","scipy"):
    importlib.import_module(m); print(m, "OK")
# confirm the cNMF API we rely on actually exists on THIS install (no assumptions):
from cnmf import cNMF
for meth in ("prepare","factorize","combine","consensus"):
    assert hasattr(cNMF, meth), f"cNMF.{meth} MISSING — tell Claude"
print("cNMF API (prepare/factorize/combine/consensus): OK")
PY
echo "=== done. Activate later with: conda activate $ENV ==="
