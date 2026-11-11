#!/bin/bash
# Checks the logs against the numbers in the paper. Run after code/run_cm.sh, code/run_scans.sh,
# code/run_models.sh and code/run_misc.sh.
# Usage: code/check_all.sh   (from any directory)
set -eu
cd "$(dirname "$0")/.."
sage code/review_gonality.sage > logs/review_gonality.log
python3 code/review_manuscript_logs.py > logs/review_manuscript_logs.log
python3 code/check_large.py > logs/check_large.log
python3 code/check_small.py > logs/check_small.log
grep -q "ALL REDUCTION CHECKS PASSED" logs/review_manuscript_reduction.log
grep -q "DONE2" logs/recheck_lemmaE.log
grep -q "idx=4 failures with 200<N: \[210\] 1" logs/gon_general.log
grep -q "0 mismatches" logs/index_formula_check.log
echo ALL_CHECKS_PASSED
