#!/usr/bin/env bash
# Reproduce project 01 (crash insurance and the G10 carry premium) from a fresh
# clone of this repository, given the private input data.
#
# The script builds a clean virtual environment from requirements.lock with the
# commands in README.md, runs the test suite, and then runs the estimation
# pipeline in the order of the reproduction blocks in
# projects/01_fx_crash_insurance_carry/reports/:
#
#   1. the Stage 1 audit of the 23 September 2026 retrieval (data_audit.md);
#   2. the one-month composite smile calibration, 30 July 2010 to 31 August 2026
#      (smile_calibration.md, e1.md);
#   3. the one-month Fenics smile calibration, 31 January 2007 to 31 August 2026
#      (e1.md);
#   4. E1 and E4 (e1.md);
#   5. Stage 4: portfolio returns, E2 and E3 (stage4.md);
#   6. E5, from the Cboe VX settlement files already on disk (stage4.md);
#   7. the moment intervals (stage4.md);
#   8. the vanna-volga calibration (robustness.md);
#   9. the robustness grid (robustness.md);
#  10. the three-month calibration, with the 24 September 2026 retrieval as a
#      supplementary raw root (robustness.md);
#  11. the three-month estimation (robustness.md).
#
# It stops at the first command that fails. It never runs an acquisition
# script and needs no network access after the packages are installed.
#
# Input data
#
# The repository contains no LSEG data, and the LSEG data must come from the
# reader's own LSEG Workspace licence. I obtained mine under a student licence
# provided by my university, which covers only me, permits individual study
# and research and does not permit redistribution, so I cannot supply them. A
# reader with an entitlement acquires them in the separate LSEG environment
# (requirements-lseg.lock; see the project README), with LSEG Workspace open
# and the reader's own App Key in ~/.lseg/app_key:
#
#   .venv-lseg/bin/python scripts/acquire_lseg_fx.py --retrieval-date 2026-09-23
#   .venv-lseg/bin/python scripts/acquire_lseg_fx.py --retrieval-date 2026-09-24 \
#       --skip-search \
#       --only AUD3M= CAD3M= CHF3M= EUR3M= GBP3M= JPY3M= NOK3M= NZD3M= SEK3M= \
#              AUD3MD= CAD3MD= CHF3MD= EUR3MD= GBP3MD= JPY3MD= NOK3MD= NZD3MD= \
#              SEK3MD= USD3MD= USD3MOIS= USDSROIS3M=
#
# The second retrieval holds the three-month forwards and rates of robustness
# variant 2 (data plan) and, like mine, makes no metadata search. The Cboe VX
# files come from Cboe's public historical data:
#
#   .venv/bin/python scripts/acquire_cboe_vx.py --retrieval-date 2026-09-24
#
# The script expects this layout, which Git ignores:
#
#   data/private/lseg/2026-09-23/   raw/, manifest.json, metadata.csv, requests.jsonl
#   data/private/lseg/2026-09-24/   raw/ with the three-month forwards and rates
#   data/private/cboe/2026-09-24/   raw/ with one settlement file per VX contract
#
# The option --retrieval-date only names the folder. The requested range is
# fixed in acquire_lseg_fx.py (1 January 1995 to 22 September 2026) and in
# acquire_cboe_vx.py (contracts from January 2007 to December 2026), so a later
# retrieval with these commands covers the same dates and lands where this
# script expects it. Its results can still differ from those reported if the
# provider has revised quotes since my retrieval or if the reader's
# entitlements differ.
#
# Usage, from the root of a fresh clone:
#
#   PYTHON=python3.13 scripts/reproduce.sh [PRIVATE_DIR]
#
# PRIVATE_DIR, if given, is a directory holding lseg/ and cboe/ in the layout
# above; the script links data/private/lseg and data/private/cboe to it and
# only reads from it. Without PRIVATE_DIR the inputs must already be in
# data/private/. The environment variable PYTHON selects the interpreter that
# builds the environment (default python3); the project needs Python 3.11 or
# later.
#
# I produced the reported results with Python 3.13.9 on macOS (Apple silicon),
# in an environment whose numpy, scipy, pandas and matplotlib were Anaconda
# builds of the versions pinned in requirements.lock, with numpy linked to
# OpenBLAS. The PyPI numpy that this script installs is linked to Apple's
# Accelerate library instead, so outputs can differ from mine in the last
# digits.
#
# On 26 September 2026 I ran this script from a fresh clone of commit 448b5a1,
# with Python 3.13.9 and my private data. All 308 tests passed and every step
# completed. The six summaries behind the reports (e1_summary.md,
# stage4_summary.md, e5_summary.md, moments_summary.md, robustness_summary.md
# and tenor3m_summary.md) were byte-identical to mine, as were the smile
# inputs, the E1, Stage 4, moment and three-month series and five of the eight
# audit tables. In the other numerical outputs every number agreed with mine to
# within 1e-9 relative plus 1e-12 absolute, apart from 212 values of a
# vanna-volga density diagnostic that differed by at most 3e-9 relative and
# changed no flag; the calibration logs differed only in the order and timing
# of lines. Three audit files differed for reasons unrelated to the estimates:
# the coverage table and the audit report now list the 21 three-month
# instruments that I added to the catalogue on 24 September, which the
# 23 September retrieval does not contain, and the conventions table carries
# the corrected citation of the delta conventions.
#
# Outputs go to data/private/audit/2026-09-23/, data/private/results/2026-09-23/
# and data/private/results/2026-09-24/. The calibration logs are written next
# to the results, as in my own runs; the other step logs, the interpreter
# version, the commit and the installed package versions go to
# data/private/reproduce/. Every output is LSEG-derived and stays under
# data/private/. The script stops rather than overwrite an existing .venv,
# existing outputs or existing links. On an eight-core Apple M2 the whole run
# took 14 minutes, of which the vanna-volga calibration took between five and
# six and the other calibrations about four and a half.

set -Eeuo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

PYTHON="${PYTHON:-python3}"
PRIVATE=data/private
LOGS="$PRIVATE/reproduce"
RESULTS="$PRIVATE/results/2026-09-23"
STEP="setup"

fail() {
  echo "reproduce.sh: $*" >&2
  exit 1
}

trap 'printf "\nreproduce.sh: stopped at step \"%s\"; its log is under %s\n" "$STEP" "$ROOT/$PRIVATE" >&2' ERR

# run NAME LOGFILE COMMAND...: run one step, sending its output to LOGFILE.
run() {
  STEP="$1"
  local log="$2"
  shift 2
  local start=$SECONDS
  printf '%-34s' "$STEP"
  "$@" > "$log" 2>&1
  echo "done in $((SECONDS - start)) s"
}

if [[ $# -gt 1 ]]; then
  fail "usage: scripts/reproduce.sh [PRIVATE_DIR]"
fi

# Refuse to overwrite anything.
if [[ -e .venv ]]; then
  fail ".venv already exists; start from a fresh clone or remove it"
fi
for out in audit/2026-09-23 results/2026-09-23 results/2026-09-24 reproduce; do
  if [[ -e "$PRIVATE/$out" ]]; then
    fail "$PRIVATE/$out already exists; move it aside so that nothing is overwritten"
  fi
done

# Check the inputs, in PRIVATE_DIR if one was given, before changing anything.
inputs="$PRIVATE"
if [[ $# -eq 1 ]]; then
  [[ -d "$1" ]] || fail "$1 is not a directory"
  inputs="$(cd "$1" && pwd)"
  for d in lseg cboe; do
    if [[ -e "$PRIVATE/$d" || -L "$PRIVATE/$d" ]]; then
      fail "$PRIVATE/$d already exists; remove it or omit PRIVATE_DIR"
    fi
  done
fi
for d in lseg/2026-09-23/raw lseg/2026-09-24/raw cboe/2026-09-24/raw; do
  [[ -d "$inputs/$d" ]] || fail "missing input $inputs/$d (see the header of this script)"
done
for f in lseg/2026-09-23/manifest.json lseg/2026-09-23/metadata.csv lseg/2026-09-23/requests.jsonl; do
  [[ -f "$inputs/$f" ]] || fail "missing input $inputs/$f (see the header of this script)"
done

command -v "$PYTHON" > /dev/null || fail "interpreter $PYTHON not found; set PYTHON to Python 3.11 or later"
"$PYTHON" -c 'import sys; sys.exit(sys.version_info < (3, 11))' \
  || fail "$PYTHON is older than Python 3.11; set PYTHON to a suitable interpreter"

# Link the inputs into data/private/ if PRIVATE_DIR was given. The pipeline
# only reads them; within data/private/ it writes only to audit/, results/
# and reproduce/.
if [[ $# -eq 1 ]]; then
  mkdir -p "$PRIVATE"
  ln -s "$inputs/lseg" "$PRIVATE/lseg"
  ln -s "$inputs/cboe" "$PRIVATE/cboe"
fi

mkdir -p "$LOGS"

# Environment (README.md).
run "virtual environment" "$LOGS/venv.log" "$PYTHON" -m venv .venv
run "install requirements.lock" "$LOGS/pip_requirements.log" .venv/bin/pip install -r requirements.lock
run "install the package" "$LOGS/pip_package.log" .venv/bin/pip install -e . --no-deps
{
  .venv/bin/python --version
  git rev-parse HEAD 2> /dev/null || echo "not a Git checkout"
  git status --short --untracked-files=no 2> /dev/null || true
  .venv/bin/pip freeze
} > "$LOGS/environment.txt"

# Tests.
run "test suite" "$LOGS/pytest.log" .venv/bin/python -m pytest -q

# Pipeline.
PY=.venv/bin/python
run "audit" "$LOGS/audit.log" "$PY" scripts/audit_fx_panel.py --retrieval-date 2026-09-23
mkdir -p "$RESULTS"
run "calibration, 1M composite" "$RESULTS/calibrate.log" \
  "$PY" scripts/calibrate_smiles.py --start 2010-07-30 --end 2026-08-31
run "calibration, 1M Fenics" "$RESULTS/calibrate_fn.log" \
  "$PY" scripts/calibrate_smiles.py --contributor FN --start 2007-01-31 --end 2026-08-31
run "E1 and E4" "$LOGS/e1.log" "$PY" scripts/estimate_e1.py
run "Stage 4" "$LOGS/stage4.log" "$PY" scripts/estimate_stage4.py
run "E5" "$LOGS/e5.log" "$PY" scripts/estimate_e5.py
run "moment intervals" "$LOGS/moments.log" "$PY" scripts/estimate_moments.py
run "vanna-volga calibration" "$LOGS/calibrate_vv.log" "$PY" scripts/calibrate_vv.py
run "robustness grid" "$LOGS/robustness.log" "$PY" scripts/robustness.py
run "calibration, 3M" "$RESULTS/calibrate_3m.log" \
  "$PY" scripts/calibrate_smiles.py --tenor 3M --extra-raw-root data/private/lseg/2026-09-24/raw
run "three-month estimation" "$LOGS/estimate_3m.log" "$PY" scripts/estimate_3m.py

STEP="finished"
echo "All steps completed. Outputs are in $PRIVATE/audit/ and $PRIVATE/results/."
