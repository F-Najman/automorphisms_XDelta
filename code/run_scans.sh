#!/bin/bash
# The scans of Sections 7 and 8 of the paper (xdelta.m and xdelta_small.m), and the check that Magma's Subgroups
# finds all the subgroups Delta (recheck_subgroups.m and recheck_subgroups.py).
# Usage: code/run_scans.sh   (from any directory; MAGMA and JOBS can be set, default magma and 6 jobs at a time)
set -eu
cd "$(dirname "$0")"
MAGMA=${MAGMA:-magma}; JOBS=${JOBS:-6}
mkdir -p ../logs/xdelta ../logs/x0
python3 levels.py x0 | tr ' ' '\n' > ../logs/x0_levels.txt
run() {  # run LOG ARGS...: one Magma job in the background, at most JOBS at a time
  while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do sleep 2; done
  local log=$1; shift
  "$MAGMA" -b "$@" < /dev/null > "$log" 2>&1 &
}
desc() { tr ' ' '\n' | sort -rn; }  # the largest levels first
for N in $(python3 levels.py x0 | desc); do run ../logs/x0/N$N.log N:=$N L:=37 B:=60 MAXIDX:=1 xdelta_small.m; done
for N in $(python3 levels.py idx4 | desc); do run ../logs/xdelta/small4_N$N.log N:=$N L:=37 B:=60 MAXIDX:=4 xdelta_small.m; done
for N in $(python3 levels.py nonsqfree | desc); do run ../logs/xdelta/small_N$N.log N:=$N L:=37 B:=60 MAXIDX:=3 xdelta_small.m; done
for N in $(echo $(python3 levels.py sqfree) $(python3 levels.py small) | desc); do
  run ../logs/xdelta/N$N.log N:=$N L:=37 B:=60 xdelta.m
done
run ../logs/subgroups_magma.log recheck_subgroups.m
wait
cd ..
python3 code/recheck_subgroups.py > logs/recheck_subgroups.log
echo SCANS_DONE
