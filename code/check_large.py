"""Checks the numbers of Sections 6 and 7 of the paper (levels N > 200) against the logs.

Run from the repository root, after code/run_cm.sh and code/run_scans.sh:
    python3 code/check_large.py
Every check is an assertion; the script ends with ALL CHECKS PASSED.
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from levels import EXACT, IDX4, x0_levels, squarefree


def numbers(path):
    return [int(x) for x in Path(path).read_text().split()]


def subgroup(n, gens):
    seen, todo = {1}, [1]
    while todo:
        a = todo.pop()
        for b in gens + [-1]:
            c = a * b % n
            if c not in seen:
                seen.add(c)
                todo.append(c)
    return frozenset(seen)


def ints(s):
    return [int(x) for x in re.findall(r"-?\d+", s)]


# --- Section 6: the constants of the uniform bound ---------------------------------------------------------------
bounds = Path("logs/review_dfs_bounds.log").read_text()
val = {k: float(v) for k, v in re.findall(r"^([^:\n]+?)(?: upper)?(?: <=)?: ([0-9.e+-]+)$", bounds, re.M)}
assert str(val["Finite range max f"]).startswith("0.9742") and "ALL CHECKS PASSED" in bounds
assert str(val["Finite range max f(p,e)/f(p,1), e>=2"]).startswith("0.6856")
small = [float(x) for x in re.findall(r"Small-prime factor \d <= ([0-9.e+-]+)", bounds)]
assert len(small) == 4 and small[0] * small[1] * small[2] * small[3] < 16
assert all(s < b for s, b in zip(small, [3.786, 2.491, 1.454, 1.146]))
assert val["Uniform majorant when some exponent >=60"] < 0.04
assert val["Majorant for discarded q^e > 10^40, e<60"] < 5e-13
theta = Path("logs/theta_check.log").read_text()
r = [float(x) for x in re.search(r"R\(2000,2\) = ([0-9.e+-]+)\s+R\(11,26\) = ([0-9.e+-]+)", theta).groups()]
assert max(r) < 0.64 and "theta(11,1) = 0.9742" in theta
print("Uniform bound: theta(11,1) = 0.9742..., max theta(q,e)/theta(q,1) = 0.6856... on the finite range,",
      "max{R(2000,2), R(11,26)} < 0.64, product of the factors for p <= 7 below 16, tail below 5e-13.")

# --- Section 6: the finite set of levels and the exact CM bound --------------------------------------------------
candidates = numbers("logs/cm_general_candidates.txt")
assert candidates == sorted(set(candidates)) and len(candidates) == 2370285 and min(candidates) == 2
assert max(candidates) == 266378112000
print("Finite set of levels:", len(candidates), "integers, the largest", max(candidates))
assert "DFS ALL CHECKS PASSED" in Path("logs/review_dfs_interval.log").read_text()
print("The interval version of the search gives the same set.")

fails = {}
for name in ["cm_exact", "cm_exact_integer"]:
    checked, fs = 0, []
    for p in sorted(Path("logs", name).glob("part_*.out")):
        text = p.read_text()
        m = re.search(r"SUMMARY checked=(\d+) fail=(\d+)", text)
        assert m, p
        f = [int(x) for x in re.findall(r"^FAIL N=(\d+)", text, re.M)]
        assert len(f) == int(m[2]), p
        checked += int(m[1])
        fs += f
    assert checked == len(candidates), name
    fails[name] = sorted(fs)
assert fails["cm_exact"] == fails["cm_exact_integer"] == numbers("logs/cm_exact_fails.txt")
fails = fails["cm_exact"]
large = [n for n in fails if n > 200]
assert len(fails) == 134 and len(large) == 58
assert [n for n in large if squarefree(n)] == [213, 285, 429]
assert [n for n in large if not squarefree(n)] == EXACT
print("Exact CM bound: 134 failures (the same for both implementations), 58 of them above 200: 213, 285, 429 and the",
      "55 levels of the paper.")

# Checks for each Delta at the 55 levels.
flagged = []
clear = []
for n in EXACT:
    log = Path(f"logs/perdelta/N{n}.log").read_text()
    m = re.search(r"PERDELTA N=\d+ subgroups=\d+ fail=(\d+)", log)
    assert m, n
    fs = re.findall(r"PERDELTA_FAIL N=\d+ gens=\[([^\]]*)\] idx=(\d+)", log)
    assert len(fs) == int(m[1]), n
    flagged += [(n, subgroup(n, ints(gs)), int(idx)) for gs, idx in fs]
    if not fs:
        clear.append(n)
assert len(flagged) == 63 and {i for _, _, i in flagged} == {2, 3, 4}
assert len({n for n, _, _ in flagged}) == 26 and len(clear) == 29 and 6912 in clear
print("Checks for each Delta: 63 subgroups with i_Delta in {2, 3, 4} left at 26 levels; the other 29 levels clear.")
assert "VIOLATIONS (V_D < proxy): []" in Path("logs/recheck_V_compare.log").read_text()
assert sorted(int(p.stem[1:]) for p in Path("logs/V").glob("N*.log")) == list(range(11, 201)) + sorted(
    [216, 243, 256, 288, 324, 384, 400, 432, 486, 512, 576, 648, 729, 768, 864, 972, 1024, 1152, 1296, 1728, 625, 375,
     675, 1250])
assert all("DONE" in p.read_text() for p in Path("logs/V").glob("N*.log"))
print("The CM bound is at least the dimension of the CM subspace at all levels 11 <= N <= 200 and at 24 levels above 200.")

# --- Section 7: the scans ----------------------------------------------------------------------------------------
direct = {}
done = {}
for p in list(Path("logs/xdelta").glob("*.log")) + list(Path("logs/x0").glob("*.log")):
    text = p.read_text()
    m = re.search(r"SUMMARY N=(\d+) ok=(\d+) fail=(\d+)", text)
    assert m and "DONE" in text, p
    done[p.parent.name + "/" + p.name] = int(m[1]), int(m[3])
    for n, gs, idx, genus, cm, status in re.findall(
            r"DELTA N=(\d+) gens=\[([^\]]*)\] idx=(\d+) genus=(\d+) cmproxy=(\d+)[^\n]* (OK|FAIL)\n", text):
        if status == "OK":
            direct[(int(n), subgroup(int(n), ints(gs)))] = int(idx), int(genus), int(cm)
# a) 200 < N <= 500: all Delta if N is square-free, i_Delta <= 3 otherwise.
for n in range(201, 501):
    name = f"xdelta/N{n}.log" if squarefree(n) else f"xdelta/small_N{n}.log"
    assert done[name] == (n, 0), n
# b) the 22 levels, i_Delta <= 4; c) the 113 levels, X_0(N).
for n in IDX4:
    assert done[f"xdelta/small4_N{n}.log"] == (n, 0), n
assert numbers("logs/x0_levels.txt") == x0_levels()
for n in x0_levels():
    assert done[f"x0/N{n}.log"] == (n, 0), n
print("Scans: all the X_Delta(N) of parts a), b) and c) of the proposition on the scans satisfy a condition",
      f"({len(x0_levels())} levels in c)).")
# The 63 subgroups left over at the 55 levels are covered by the scans of parts a) and b).
for n, d, idx in flagged:
    assert (n, d) in direct and direct[(n, d)][0] == idx, (n, sorted(d))
print("Each of the 63 subgroups is covered by a scan, as sets of residues.")
# The numbers quoted after the proposition.
full = subgroup(6912, [x for x in range(1, 6912) if x % 2 and x % 3])
assert direct[(6912, full)] == (1, 1081, 236)
at4608 = sorted(direct[(n, d)][1:] for n, d, idx in flagged if n == 4608)
assert at4608 == [(1409, 392), (1409, 428)]
print("N = 6912: g(X_0(6912)) = 1081 and CM bound 236; N = 4608: the two subgroups have genus 1409 and CM bounds 392, 428.")
for line in Path("logs/recheck_subgroups.log").read_text().splitlines():
    if line.startswith("compared"):
        assert line.endswith("mismatches: []"), line
        print("Magma's Subgroups:", line)
print("ALL CHECKS PASSED")
