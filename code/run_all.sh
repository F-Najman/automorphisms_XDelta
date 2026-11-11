#!/bin/bash
# Runs every computation of the paper and then checks the logs (code/check_all.sh).
# Usage: code/run_all.sh   (from any directory; MAGMA and JOBS can be set, default magma and 6 jobs at a time)
set -eu
cd "$(dirname "$0")"
mkdir -p ../logs
./versions.sh > ../logs/versions.txt
./run_cm.sh
./run_scans.sh
./run_models.sh
./run_misc.sh
./check_all.sh
