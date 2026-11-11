#!/bin/bash
# Section 6 of the paper: the local factors, the finite set of levels, the exact CM bound at these levels (two
# implementations), the checks for each Delta at the remaining 55 levels, and the comparison of the bound with
# the dimensions of the CM subspaces at 24 levels.
# Usage: code/run_cm.sh   (from any directory; MAGMA and JOBS can be set, default magma and 6 jobs at a time)
set -eu
cd "$(dirname "$0")/.."
MAGMA=${MAGMA:-magma}; JOBS=${JOBS:-6}
run() {  # run CMD: one shell command in the background, at most JOBS at a time
  while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do sleep 2; done
  bash -c "$1" &
}
mkdir -p logs/cm_exact logs/cm_exact_integer logs/perdelta logs/V
python3 code/cmlocal.py > logs/cmlocal.log
python3 code/theta_check.py > logs/theta_check.log
sage code/review_dfs_bounds.sage > logs/review_dfs_bounds.log
python3 code/cm_general_dfs.py > logs/cm_general_dfs.log     # writes logs/cm_general_candidates.txt
run "sage code/review_dfs_interval.sage > logs/review_dfs_interval.log 2>&1"
tmp=$(mktemp -d)
gcc -O2 -fopenmp -o "$tmp/cm_exact" code/cm_exact.c -lm
gcc -O2 -ftrapv -o "$tmp/cm_exact_integer" code/review_cm_exact_integer.c -lm
split -d -n r/8 logs/cm_general_candidates.txt "$tmp/part_"
for p in "$tmp"/part_0?; do
  b=$(basename "$p")
  run "'$tmp/cm_exact' < '$p' > logs/cm_exact/$b.out 2> logs/cm_exact/$b.err"
  run "'$tmp/cm_exact_integer' < '$p' > logs/cm_exact_integer/$b.out 2> logs/cm_exact_integer/$b.err"
done
for N in $(python3 code/levels.py exact); do
  run "python3 code/cmV_table.py $N > '$tmp/table_$N.txt' && cd code && '$MAGMA' -b N:=$N T:='$tmp/table_$N.txt' perdelta.m < /dev/null > ../logs/perdelta/N$N.log 2>&1"
done
for c in 216:6 243:9 256:16 288:8 324:6 384:8 400:4 432:6 486:9 512:8 576:4 648:6 729:3 768:4 864:4 972:3 1024:4 1152:2 1296:2 1728:2 625:4 375:4 675:4 1250:2; do
  N=${c%%:*}; M=${c#*:}
  run "cd code && '$MAGMA' -b N:=$N MAXORD:=$M B:=200 recheck_V.m < /dev/null > ../logs/V/N$N.log 2>&1"
done
for N in $(seq 11 200); do
  run "cd code && '$MAGMA' -b N:=$N MAXORD:=1000 B:=200 recheck_V.m < /dev/null > ../logs/V/N$N.log 2>&1"
done
wait
rm -rf "$tmp"
cat logs/cm_exact/part_0?.out | sed -n 's/^FAIL N=\([0-9]*\) .*/\1/p' | sort -n > logs/cm_exact_fails.txt
python3 code/recheck_V_compare.py > logs/recheck_V_compare.log
echo CM_DONE
