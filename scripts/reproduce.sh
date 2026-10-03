#!/usr/bin/env bash
# Reproduce the project on crash insurance and the G10 carry premium from a fresh
# clone of this repository, given the private input data.
#
# The script builds a clean virtual environment from requirements.lock with the
# commands in docs/development.md, runs the test suite, and then runs the estimation
# pipeline: the steps of the reproduction blocks in reports/, in their order,
# followed by the checks that the paper reports in its Section 4.7:
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
#  11. the three-month estimation (robustness.md);
#  12. to 14., added on 27 September 2026: the post hoc regime attribution (e1.md), the model validation of
#      the moment code (model_validation.md, from synthetic models only) and the
#      sharp identification of option-implied moments (stage4.md), which takes
#      about 45 minutes on eight cores;
#  15. added on 27 September 2026: the out-of-sample test of the tail
#      hypotheses with the 5-delta quotes (stage4.md), about 15 minutes;
#  16. added on 27 September 2026: the sign and magnitude check of HML_FX
#      against Verdelhan's public portfolios (paper, Section 4.6);
#  17. added on 27 September 2026: the E1 bootstrap intervals recomputed with
#      the arch package as an independent implementation (paper, Section 4.6).
#  18. and 19., added on 28 September 2026: the long at-the-money sample
#      (stage4.md, about 20 seconds) and the design checks (e1.md and
#      stage4.md; SOFR OIS, the named-broker contributor, delta-method
#      intervals and the moment code's closed-form checks, about 2 minutes on
#      eight cores);
#  20. and 21., added on 2 October 2026: the design items completed then
#      (E3 by regime, stale Fenics butterflies in the extended sample and the
#      9,999-draw comparison; stage4.md and robustness.md), about 1 minute,
#      and the descriptive account of 5 August 2024 (stage4.md), about
#      30 seconds;
#  22. added on 3 October 2026: the checks behind statements of the paper
#      and the theory notes that no earlier step produced (the hypotheses of
#      R5(c) and R8 on the calibrated smiles, a denser grid of means in the
#      identified sets, and the 5-delta breakdown comparison on the same
#      currency-months; stage4.md), about 1 minute.
#
# Before the tests, the C++ kernel is built with cpp/build.sh, which needs a
# C++17 compiler. scripts/benchmark_kernel.py is not run, because it rewrites
# the timing block of the tracked cpp/README.md.
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
# (requirements-lseg.lock; see docs/development.md), with LSEG Workspace open
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
# variant 2 (data plan) and, like mine, makes no metadata search. A third
# retrieval, of 27 September 2026, holds the one-month 5-delta Fenics quotes of
# the out-of-sample test (step 15):
#
#   .venv-lseg/bin/python scripts/acquire_lseg_fx.py --retrieval-date 2026-09-27 \
#       --skip-search \
#       --only EUR1MR5=FN EUR1MB5=FN GBP1MR5=FN GBP1MB5=FN AUD1MR5=FN AUD1MB5=FN \
#              NZD1MR5=FN NZD1MB5=FN JPY1MR5=FN JPY1MB5=FN CHF1MR5=FN CHF1MB5=FN \
#              CAD1MR5=FN CAD1MB5=FN NOK1MR5=FN NOK1MB5=FN SEK1MR5=FN SEK1MB5=FN
#
# An optional fourth retrieval, of 1 October 2026, checks the coverage of the
# long at-the-money sample's series before 1995 (step 18 reads its manifest if
# present, and its results do not depend on it):
#
#   .venv-lseg/bin/python scripts/acquire_lseg_fx.py --retrieval-date 2026-10-01 \
#       --start 1970-01-01 --end 1995-01-31 --skip-search \
#       --only <CCY>= <CCY>1M= <CCY>1MO=   (for each of the nine currencies)
#
# The Cboe VX files come from Cboe's public historical data:
#
#   .venv/bin/python scripts/acquire_cboe_vx.py --retrieval-date 2026-09-24
#
# Verdelhan's public currency portfolios (step 16) are downloaded with
#
#   .venv/bin/python scripts/acquire_verdelhan.py --retrieval-date 2026-09-27
#
# The script expects this layout, which Git ignores:
#
#   data/private/lseg/2026-09-23/   raw/, manifest.json, metadata.csv, requests.jsonl
#   data/private/lseg/2026-09-24/   raw/ with the three-month forwards and rates
#   data/private/lseg/2026-09-27/   raw/ with the one-month 5-delta quotes (step 15)
#   data/private/lseg/2026-10-01/   optional: spot, forwards and ATM before 1995 (step 18)
#   data/private/verdelhan/2026-09-27/   CurrencyPortfolios.xls and its manifest (step 16)
#   data/private/cboe/2026-09-24/   raw/ with one settlement file per VX contract
#
# The option --retrieval-date only names the folder. The requested range
# defaults in acquire_lseg_fx.py to 1 January 1995 to 22 September 2026
# (--start and --end change it, as for the coverage retrieval above) and defaults in
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
# PRIVATE_DIR, if given, is a directory holding lseg/, cboe/ and verdelhan/ in
# the layout above; the script links those folders of data/private/ to it and
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
# On 3 October 2026 I ran this script, with every step listed above, from a
# fresh clone of the repository in a new environment, with Python 3.13.9 and
# my private data. All tests passed (442, with the 50-digit check of the
# kernel skipped because mpmath is not in requirements.lock) and every step
# completed. Of the 50 result files, 23 were byte-identical to mine, and in
# the others every number agreed with mine to within 1e-9 relative plus 1e-12
# absolute, apart from run times, 216 values of a vanna-volga density
# diagnostic that differed by at most 4.2e-10 and changed no flag, and one
# bootstrap block length that differed in its sixteenth significant digit.
# Five of the eight audit tables were byte-identical. The other three differed
# for reasons unrelated to the estimates: the coverage table and the audit
# report list instruments added to the catalogue after the 23 September
# retrieval, which that retrieval does not contain, and the conventions table
# carries the corrected citation of the delta conventions. An earlier run on
# 26 September, of steps 1 to 11, agreed in the same way.
#
# Outputs go to data/private/audit/2026-09-23/, data/private/results/2026-09-23/
# and data/private/results/2026-09-24/. The calibration logs are written next
# to the results, as in my own runs; the other step logs, the interpreter
# version, the commit and the installed package versions go to
# data/private/reproduce/. Every output is LSEG-derived and stays under
# data/private/. The script stops rather than overwrite an existing .venv,
# existing outputs or existing links. On an eight-core Apple M2 the whole run
# took between 47 and 88 minutes, depending on the other load on the machine,
# most of it in the sharp identification (27 to 65 minutes); the calibrations
# took about six minutes together.

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
  for d in lseg cboe verdelhan; do
    if [[ -e "$PRIVATE/$d" || -L "$PRIVATE/$d" ]]; then
      fail "$PRIVATE/$d already exists; remove it or omit PRIVATE_DIR"
    fi
  done
fi
for d in lseg/2026-09-23/raw lseg/2026-09-24/raw lseg/2026-09-27/raw cboe/2026-09-24/raw verdelhan/2026-09-27; do
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
  ln -s "$inputs/verdelhan" "$PRIVATE/verdelhan"
fi

mkdir -p "$LOGS"

# Environment (docs/development.md).
run "virtual environment" "$LOGS/venv.log" "$PYTHON" -m venv .venv
run "install requirements.lock" "$LOGS/pip_requirements.log" .venv/bin/pip install -r requirements.lock
run "install the package" "$LOGS/pip_package.log" .venv/bin/pip install -e . --no-deps
{
  .venv/bin/python --version
  git rev-parse HEAD 2> /dev/null || echo "not a Git checkout"
  git status --short --untracked-files=no 2> /dev/null || true
  .venv/bin/pip freeze
} > "$LOGS/environment.txt"

# Tests. The C++ kernel is built first; its parity tests would otherwise build it (or skip without a compiler).
run "C++ kernel build" "$LOGS/cpp_build.log" bash cpp/build.sh
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
# Post hoc analyses added on 27 September 2026.
run "regime attribution (post hoc)" "$LOGS/regime_attribution.log" "$PY" scripts/estimate_regime_attribution.py
run "model validation (post hoc)" "$LOGS/model_validation.log" "$PY" scripts/validate_moments_models.py
run "sharp identification (post hoc)" "$LOGS/identification.log" "$PY" scripts/estimate_identification.py
run "5-delta out-of-sample test (post hoc)" "$LOGS/wing_test.log" "$PY" scripts/estimate_wing_test.py
run "comparison with published portfolios" "$LOGS/verdelhan.log" "$PY" scripts/compare_verdelhan.py
run "bootstrap check against arch" "$LOGS/bootstrap_arch.log" "$PY" scripts/check_bootstrap_arch.py
# Items of the design computed after the main results (research log, 27 and 28 September 2026).
run "long at-the-money sample" "$LOGS/long_atm.log" "$PY" scripts/estimate_long_atm.py
run "design checks" "$LOGS/design_checks.log" "$PY" scripts/estimate_design_checks.py
# Items of the design completed on 2 October 2026 (research log, 2 October 2026).
run "outstanding design items" "$LOGS/outstanding_items.log" "$PY" scripts/estimate_outstanding_items.py
run "5 August 2024 account" "$LOGS/august_2024.log" "$PY" scripts/describe_august_2024.py
# Checks behind reported statements (research log, 3 October 2026).
run "checks behind reported statements" "$LOGS/reported_conditions.log" "$PY" scripts/check_reported_conditions.py
# scripts/benchmark_kernel.py is not run here: it rewrites the timing block of the tracked cpp/README.md.

STEP="finished"
echo "All steps completed. Outputs are in $PRIVATE/audit/ and $PRIVATE/results/."
