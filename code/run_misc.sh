#!/bin/bash
# Section 5 (multiplicities of T_l^n(P) on X(1)), Section 7 (the lemma on gonality) and Section 10 (the index of
# the normalizer of Gamma_Delta(N), N <= 400) of the paper.
# Usage: code/run_misc.sh   (from any directory; JOBS can be set, default 6)
set -eu
cd "$(dirname "$0")/.."
JOBS=${JOBS:-6}
python3 code/recheck_lemmaE.py > logs/recheck_lemmaE.log
python3 code/gon_general.py > logs/gon_general.log
python3 code/aut_orders.py 13 200 "$JOBS" > logs/aut_orders_13_200.tsv 2> logs/aut_orders_13_200.err
python3 code/aut_orders.py 201 400 "$JOBS" > logs/aut_orders_201_400.tsv 2> logs/aut_orders_201_400.err
python3 code/index_formula_check.py logs/aut_orders_13_200.tsv logs/aut_orders_201_400.tsv > logs/index_formula_check.log
echo MISC_DONE
