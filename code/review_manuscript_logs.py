"""Checks of the automorphism groups over finite fields of Section 8.1 of the paper (the logs of canon2.m).

Run from the repository root, after code/run_models.sh:
    python3 code/review_manuscript_logs.py
For every canon2.m run behind a number in Table 1 or behind the curves X_1(N), N in {21,22,24,30}, it asserts:
  * the run finished (RESULT line);
  * there are (g-2)(g-3)/2 quadrics, the model over F_ell has genus g, and the reduction is not hyperelliptic;
  * the point count agrees with ell + 1 - tr T_ell (canon2.m prints this but does not assert it);
  * g <= 10, so every m with phi(m) <= 4g^2 is at most 1680, below the cutoff 2000 in canon2.m.
It also recomputes two constants of the proof of the lemma bounding theta(q, e) in Section 6.
Table 1 is read from code/table_small.tsv.
"""
from pathlib import Path
import math
import re
import sys

from sympy import totient



def subgroup(n, generators):
    seen = {1}
    todo = [1]
    while todo:
        a = todo.pop()
        for b in [-1] + generators:
            c = a * b % n
            if c not in seen:
                seen.add(c)
                todo.append(c)
    return frozenset(seen)


# Table 1 rows whose last column is a number (computed as in Section 8.1).
numbered = {}
for row in Path("code/table_small.tsv").read_text().splitlines():
    fields = row.split("\t")
    if row.startswith("#") or not fields[0].isdigit():
        continue
    n, last = int(fields[0]), fields[4]
    if "," not in last:
        numbered[(n, subgroup(n, [int(x) for x in fields[1].split(",")]))] = last
assert len(numbered) == 30, len(numbered)

# Curves not computed with canon2.m: X_{+-1,+-8}(21) (hyperelliptic, explicit model) and level 63 (Section 9).
direct = {(21, subgroup(21, [13, 8])), (63, subgroup(63, [59, 37]))}
assert direct <= set(numbered)

logs = sorted(Path("logs/canon2").glob("*.log"))

covered = set()
x1_levels = set()
max_g = 0
for p in logs:
    text = p.read_text()
    m = re.search(r"RESULT N=(\d+) Delta=\[([^\]]*)\] g=(\d+) ell=(\d+) k=(\d+) : \|Aut\(X_Fbar\)\| = (\d+)", text)
    assert m, p
    n, delta, g = int(m.group(1)), frozenset(int(x) for x in m.group(2).split(",")), int(m.group(3))
    assert f"genus of model over F_ell: {g}" in text, p
    assert "HYPERELLIPTIC" not in text, p
    quads = re.findall(r"number of quadrics: (\d+) ", text)
    assert quads == [str((g - 2) * (g - 3) // 2)], (p, quads)
    counts = re.findall(r"#C\(F_ell\) = (\d+) \(predicted (\d+)\)", text)
    assert len(counts) == 1 and counts[0][0] == counts[0][1], p
    max_g = max(max_g, g)
    if delta == frozenset({1, n - 1}):
        x1_levels.add(n)
    else:
        key = (n, delta)
        assert key in numbered, p
        aut = int(m.group(6))
        expected = numbered[key]
        assert expected == str(aut) or expected.startswith(f"{aut}/"), (p, expected, aut)
        covered.add(key)

assert covered == set(numbered) - direct, set(numbered) - direct - covered
assert x1_levels == {21, 22, 24, 30}, x1_levels
print(f"canon2 runs: {len(logs)} logs, {len(covered)} Table 1 curves and X_1(N) for N in {sorted(x1_levels)}.")
print("All runs: (g-2)(g-3)/2 quadrics, model genus = g, no hyperelliptic reduction, point count = ell + 1 - tr T_ell.")

cutoff = max(m for m in range(1, 2 * (4 * max_g**2) ** 2 + 1) if totient(m) <= 4 * max_g**2)
assert cutoff < 2000
print(f"max genus {max_g}; largest m with phi(m) <= {4 * max_g**2} is {cutoff} < 2000.")

sys.path.insert(0, "code")
from cmlocal import Lam, psi_pe


def f(q, e):
    return Lam(q, e) / psi_pe(q, e) * (1 + e * math.log(q) / 2)


print("f(11,1) =", f(11, 1))
print("f(11,2)/f(11,1) =", f(11, 2) / f(11, 1))
print("REVIEW LOG CHECKS COMPLETED")
