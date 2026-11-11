# Automorphisms of the modular curves X_Δ(N): code and logs

Code and logs for the computations in

> F. Najman, *Automorphisms of the modular curves X_Δ(N)* (in preparation).

For a subgroup Δ of (Z/NZ)^× containing −1, let X_Δ(N) be the modular curve of Γ_Δ(N), where
Γ_1(N) ⊆ Γ_Δ(N) ⊆ Γ_0(N). The paper proves that every automorphism of a curve X_Δ(N) of genus at least 2 is
modular, that is, induced by the normalizer of Γ_Δ(N) in PSL₂(R), except for five curves. The proof needs
computations for finitely many levels:

* a bound for the dimension of the CM part of J_Δ(N), checked at the levels of an explicit finite set;
* checks of the conditions of Proposition 4.1 for all the curves of level N ≤ 200 and some of larger level
  (the *scans*);
* the automorphism groups of explicit models of some curves of small level.

This repository contains all of these computations, their logs and scripts that check every number in the paper
against the logs. Statement numbers refer to the version of the paper of November 2026.

The results of Section 10 on the normalizer of Γ_Δ(N) are also formalized in Lean 4, see
<https://github.com/F-Najman/normalizer-lean>.

## Requirements

* Magma. The logs were produced with V2.29-4.
* SageMath (10.7). The Python scripts are run with Sage's `python3` (3.11), which provides `sympy`, `mpmath` and
  `cypari2`.
* gcc (15.2), with `__int128`.
* bash and GNU coreutils.

## Running

Run the drivers from the repository root. Each one changes into `code/` where Magma needs it: the Magma scripts
`load "lib.m"`, so they are run from `code/`. The other scripts are run from the repository root.

    code/run_all.sh      # everything below, in this order, then code/check_all.sh

    code/run_cm.sh       # Section 6: the CM bound
    code/run_scans.sh    # Sections 7 and 8: the scans
    code/run_models.sh   # Sections 8 and 9, Remark 9.3: explicit models, descent, level 63
    code/run_misc.sh     # Sections 5, 7 and 10
    code/check_all.sh    # checks the logs against the numbers in the paper

`MAGMA=/path/to/magma` selects the Magma binary and `JOBS=n` the number of jobs run at the same time (default 6).
The drivers write into `logs/` and overwrite the logs there. The four drivers were run at the same time, with
`JOBS` = 3, 5, 2 and 4, on a 32-thread AMD EPYC 9175F that was also running other jobs. `run_cm.sh` took about 50
minutes, `run_scans.sh` 70 minutes (most of it for the scan of level 4608 with i_Δ ≤ 4), `run_models.sh` 35
minutes and `run_misc.sh` 15 minutes.

A single computation is run directly. Parameters are passed as `Name:=value`, for example

    cd code && magma -b N:=63 L:=37 B:=60 xdelta.m < /dev/null    # the scan of level 63

Each Magma script states its parameters in its first lines, and the drivers show the values used.

`MANIFEST.sha256` lists the SHA-256 checksums of all the other files. Check it with `sha256sum -c MANIFEST.sha256`.

## The computations of the paper

| paper | what is computed | code | logs |
|---|---|---|---|
| §5, before Lemma 5.4 | the support and the multiplicities of T_ℓ^n(P) on X(1) for elliptic, CM and generic points | `recheck_lemmaE.py` | `recheck_lemmaE.log` |
| §6.1, after Proposition 6.3 | the bound B(χ) of Proposition 6.3 against the dimension of the subspace of S_2(N, χ) killed by T_q for the primes q ≤ 200 inert in F, for 11 ≤ N ≤ 200 (there are no cusp forms on Γ_1(N) for N ≤ 10 and N = 12) and 24 levels above 200 | `recheck_V.m`, `cmV.py`, `recheck_V_compare.py` | `V/`, `recheck_V_compare.log` |
| §6.2, before Lemma 6.6 | Λ_p(1) and Λ_p(2) for p = 3, 5, 7, 11, and Λ_p(e)/ψ(p^e) for p ≤ 7 and small e | `cmlocal.py` | `cmlocal.log` |
| Lemma 6.7 | θ(11, 1) = 0.9742…, θ(q, e)/θ(q, 1) ≤ 0.6856… for 11 ≤ q < 2000 and 2 ≤ e ≤ 25 (interval arithmetic), and R(2000, 2), R(11, 26) | `review_dfs_bounds.sage`, `theta_check.py` | `review_dfs_bounds.log`, `theta_check.log` |
| Proposition 6.8 | the set 𝒮 (2 370 285 levels, the largest 266 378 112 000), and the constants 3.786, 2.491, 1.454, 1.146, 0.04 and 5·10⁻¹³ of the proof; a second search with interval arithmetic | `cm_general_dfs.py`, `review_dfs_bounds.sage`, `review_dfs_interval.sage` | `cm_general_candidates.txt`, `cm_general_dfs.log`, `review_dfs_bounds.log`, `review_dfs_interval.log` |
| Proposition 6.9 | max_χ B(χ) at every level of 𝒮: 134 failures, 58 of them above 200; a second implementation in integer arithmetic | `cm_exact.c`, `review_cm_exact_integer.c` | `cm_exact/`, `cm_exact_integer/`, `cm_exact_fails.txt` |
| Proposition 6.10 | the left-hand side of (6.2) for every Δ at the 55 levels | `perdelta.m`, `cmV_table.py`, `cmV.py` | `perdelta/` |
| §7, before Proposition 7.1 | the scan of a level: the conditions of Proposition 4.1 for every Δ | `xdelta.m`, `xdelta_small.m` | |
| Proposition 7.1 | a) the scans of 200 < N ≤ 500; b) the scans with i_Δ ≤ 4 at 22 levels; c) the scans of X_0(N) at 113 levels; and the numbers for N = 6912 and N = 4608 after it | `xdelta.m`, `xdelta_small.m`, `levels.py` | `xdelta/N*.log`, `xdelta/small_N*.log`, `xdelta/small4_N*.log`, `x0/`, `x0_levels.txt` |
| Lemma 7.2 | the levels 200 < N < 2743 where the inequalities fail (105 levels in b), N = 210 in c)) | `gon_general.py`, `review_gonality.sage` | `gon_general.log`, `review_gonality.log` |
| proof of Theorem 1.1 for N > 200 | the 63 subgroups of Proposition 6.10 are covered by Proposition 7.1, as sets of residues | `check_large.py` | `check_large.log` |
| §§6–8 | Magma's `Subgroups` finds all the subgroups Δ at the levels used | `recheck_subgroups.m`, `recheck_subgroups.py` | `subgroups_magma.log`, `recheck_subgroups.log` |
| §8 | the scans of 13 ≤ N ≤ 200: 1449 curves of genus at least 2, of which 1362 satisfy a condition; the 37 curves X_0(N), the 8 curves X_1(N) and the 42 curves of Table 1 that do not | `xdelta.m` | `xdelta/N13.log`, …, `xdelta/N200.log`, `check_small.log` |
| §8.1 | Aut of the models of X_1(13), X_1(16), X_1(18), X_1(20) and X_{±1,±8}(21) over Q̄ | `t2.m`, `t_small.m`, `t3.m`, `t4.m`, `hyp.m` | `models/` |
| §8.1, Table 1 | Aut over finite fields: canonical models, point counts and the field of definition of the automorphisms | `canon2.m`, `review_manuscript_reduction.m`, `review_manuscript_logs.py` | `canon2/`, `review_manuscript_reduction.log`, `review_manuscript_logs.log` |
| §8.1, Table 1 | the modular automorphisms from explicit matrices (in the normalizer of Γ_0(N)) | `normalizer.m` | `normalizer.log` |
| Lemma 8.2, Table 1 | conditions a) and b) of the descent lemma for the 12 curves marked a) or b) | `recheck_c4.m` | `c4/` |
| Proposition 9.1 | condition a) of Lemma 8.2 at level 63; Aut of the models of X and Y over F_{127²} and the field of definition of their automorphisms; the 48 modular automorphisms; the 24 automorphisms of X_0(63) preserving the cusps | `recheck_c4.m`, `review_xdelta63.m`, `review_frobenius_fields.m`, `normalizer.m`, `x063.m` | `c4/N63_D59_37.log`, `models/`, `normalizer.log` |
| Remark 9.3 | the order 21 of the cuspidal class on X_{Δ37}(37) modulo 3, 5, 7, 11, 13; the Atkin–Lehner involutions not normalizing Γ_Δ(N), square-free N ≤ 70 | `t37c.m`, `enum.m`, `t65.m` | `models/` |
| §10, Table 2 | the index of the normalizer of Γ_Δ(N) for the 4397 curves with N ≤ 400, by testing the matrices σH, and the formula of Theorem 10.6 | `aut_orders.py`, `index_formula_check.py` | `aut_orders_13_200.tsv`, `aut_orders_201_400.tsv`, `index_formula_check.log` |
| Tables 1 and 2 | every entry, against the logs above | `check_small.py`, `table_small.tsv`, `table_order.tsv` | `check_small.log` |

The logs are in `logs/`. `code/check_large.py` (Sections 6 and 7) and `code/check_small.py` (Sections 8–10 and
the tables) assert every number of the paper that comes from a computation, and `code/check_all.sh` runs them and
the other checks.

### Reading the logs

* The scans print one line per Δ:
  `DELTA N=… gens=[…] idx=i_Δ genus=g_Δ cmproxy=… C1=… C3_l=… C3p_l=… C2=… OK|FAIL`. Here `cmproxy` is the upper
  bound dim W for g_CM(Δ) of Section 7, `C1` is the condition (CM), `C3_l` and `C3p_l` are the least ℓ for which
  conditions a) and b) of Proposition 4.1 hold (0 if there is none), and `C2` is condition c), with the gonality
  bound (4.2). The line says `OK` if one of a), b) and c) holds.
* `recheck_c4.m` prints a line `C4CHECK` for every Δ' ⊋ Δ. `coprime_any_p` is the least prime p ∤ N at which the
  characteristic polynomials of ⟨p⟩⁻¹T_p² on J^H and on A are coprime (0 if there is none), and `cm'(bound)` and
  `cm_complement` are upper bounds for the parts of J^H and A with CM by each field. Condition a) of Lemma 8.2
  holds if `coprime_any_p` is non-zero and no field occurs on both sides (`overlapD=[]`). For b), dim S is at most
  the sum of both bounds over the fields in `overlapD`. The labels `FULL`, `NCM+C1` and `FAIL` at the end of the
  line are not the conditions of Lemma 8.2; `code/check_small.py` checks a) and b) as just described.
* `canon2.m` prints a line `RESULT` with #Aut(X_{F̄_ℓ}) and the number |B_Δ| of automorphisms coming from
  diamond operators and Atkin–Lehner operators.

The names of the scripts that start with `review_` or `recheck_` come from the independent checks and second
implementations that were written while the computations were reviewed. The paper cites them under these names.

## Verification

The logs in `logs/` were produced by the drivers above from a clean checkout of this repository, in November 2026,
with Magma V2.29-4, SageMath 10.7, Python 3.11 and gcc 15.2 (`logs/versions.txt`). `code/check_all.sh` passes on
them.

The computations were first run in September 2026, with the same scripts, before they were collected in this
repository. Apart from timings, the new logs are identical to the logs of those runs wherever those exist, and the
two implementations of Proposition 6.9 give the same 134 levels. The runs of `recheck_V.m` at the levels N ≤ 200,
and of `normalizer.m` for all the rows of Table 1 (the earlier runs covered 15 of them), have no earlier logs.

## Layout

| path | contents |
|---|---|
| `code/lib.m`, `code/hyp.m` | Magma helpers: characters trivial on Δ, q-expansion bases of S_2(Γ_Δ(N)), relations; hyperelliptic models |
| `code/xdelta.m`, `code/xdelta_small.m` | the scan of a level (`xdelta_small.m` only for i_Δ ≤ `MAXIDX`) |
| `code/levels.py` | the lists of levels of Sections 6–8 |
| `code/cmlocal.py`, `code/cmV.py`, `code/cmV_table.py` | the local factors of the CM bound, B(χ) for one type of χ, and the table of B(χ) for all types at a level |
| `code/cm_general_dfs.py` | the search for the set 𝒮 of Proposition 6.8 |
| `code/cm_exact.c`, `code/review_cm_exact_integer.c` | Proposition 6.9, in floating point and in integer arithmetic |
| `code/perdelta.m` | the check of (6.2) for every Δ at one level |
| `code/recheck_V.m`, `code/recheck_V_compare.py` | B(χ) against the dimensions of the subspaces killed by inert T_q |
| `code/review_dfs_bounds.sage`, `code/review_dfs_interval.sage`, `code/theta_check.py` | Lemma 6.7 and Proposition 6.8 with interval arithmetic, and the constants of the proof of Lemma 6.7 |
| `code/gon_general.py`, `code/review_gonality.sage` | Lemma 7.2 |
| `code/recheck_subgroups.m`, `code/recheck_subgroups.py` | the number of subgroups Δ found by Magma against Birkhoff's formula |
| `code/canon2.m`, `code/review_manuscript_reduction.m`, `code/review_manuscript_logs.py` | Aut of a canonical model over a finite field, and its checks |
| `code/normalizer.m` | the modular automorphisms of X_Δ(N) coming from the normalizer of Γ_0(N) |
| `code/recheck_c4.m` | the descent lemma |
| `code/t2.m`, `code/t3.m`, `code/t4.m`, `code/t_small.m` | the explicit models of §8.1 |
| `code/review_xdelta63.m`, `code/review_frobenius_fields.m`, `code/x063.m` | level 63 |
| `code/t37c.m`, `code/enum.m`, `code/t65.m` | levels 37 and 65 (Remark 9.3) |
| `code/aut_orders.py`, `code/index_formula_check.py` | Section 10 |
| `code/table_small.tsv`, `code/table_order.tsv` | Tables 1 and 2 of the paper |
| `code/check_large.py`, `code/check_small.py`, `code/check_all.sh` | the checks of the logs against the paper |
| `code/run_*.sh`, `code/versions.sh`, `code/make_manifest.sh` | the drivers, the versions of the programs, the manifest |
| `logs/` | the outputs, named as in the table above |
