#!/bin/bash
# Sections 8 and 9 of the paper and the remark on Momose's preprint: the explicit models of X_1(13), X_1(16),
# X_1(18), X_1(20) and X_{+-1,+-8}(21), the automorphism groups over finite fields (canon2.m) and their checks, the
# modular automorphisms from explicit matrices (normalizer.m), the descent lemma (recheck_c4.m), the curves of level
# 63, and the computations at levels 37 and 65.
# Usage: code/run_models.sh   (from any directory; MAGMA and JOBS can be set, default magma and 6 jobs at a time)
set -eu
cd "$(dirname "$0")"
MAGMA=${MAGMA:-magma}; JOBS=${JOBS:-6}
run() {  # run LOG ARGS...: one Magma job in the background, at most JOBS at a time
  while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do sleep 2; done
  local log=$1; shift
  "$MAGMA" -b "$@" < /dev/null > "$log" 2>&1 &
}
mkdir -p ../logs/models ../logs/canon2 ../logs/c4
rows() { grep -v '^#' table_small.tsv | tail -n +2; }  # Table 1: the lines "N gens i g entry" of table_small.tsv
# Automorphism groups over finite fields: the rows whose entry is a number (or m/k), except X_{+-1,+-8}(21)
# (explicit model below) and level 63 (Section 9); and X_1(N) for N = 21, 22, 24, 30.
while IFS=$'\t' read -r N D i g entry; do
  case "$entry" in *")"*) continue;; esac
  if [ "$N:$D" = "21:13,8" ] || [ "$N" = 63 ]; then continue; fi
  run ../logs/canon2/N${N}_D${D//,/_}.log N:=$N D:=$D canon2.m
done < <(rows)
for N in 21 22 24 30; do run ../logs/canon2/X1_N$N.log N:=$N D:=1 canon2.m; done
# The descent lemma: the 12 rows marked a) or b), and level 63 (Section 9).
while IFS=$'\t' read -r N D i g entry; do
  case "$entry" in *")"*) run ../logs/c4/N${N}_D${D//,/_}.log N:=$N D:=$D recheck_c4.m;; esac
done < <(rows)
run ../logs/c4/N63_D59_37.log N:=63 D:=59,37 recheck_c4.m
# Explicit models.
for s in t2 t3 t4 t_small; do run ../logs/models/$s.log $s.m; done
# Level 63 (Section 9), level 37 and level 65 (the remark on Momose's preprint).
for s in review_xdelta63 review_frobenius_fields x063 t37c t65; do run ../logs/models/$s.log $s.m; done
run ../logs/models/enum.log Nmax:=70 enum.m
wait
# Checks of the models over finite fields, and the modular automorphisms of all the curves of Table 1.
run ../logs/review_manuscript_reduction.log review_manuscript_reduction.m
: > ../logs/normalizer.log
while IFS=$'\t' read -r N D i g entry; do
  "$MAGMA" -b N:=$N D:=$D normalizer.m < /dev/null | grep '^NORMALIZER' >> ../logs/normalizer.log
done < <(rows)
wait
echo MODELS_DONE
