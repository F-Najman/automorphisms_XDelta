"""The lists of levels used in Sections 6-8 of the paper.

Usage: python3 code/levels.py NAME, where NAME is one of
  small      13 <= N <= 200 (scans of all X_Delta(N), Section 8)
  sqfree     the square-free N with 200 < N <= 500 (scans of all X_Delta(N))
  nonsqfree  the other N with 200 < N <= 500 (scans of the X_Delta(N) with i_Delta <= 3)
  idx4       the 22 levels of part b) of the proposition on the scans (scans with i_Delta <= 4)
  x0         the 113 levels of part c) of the same proposition (scans of X_0(N))
  exact      the 55 levels of the proposition on the exact CM bound (checks for each Delta)
The lists idx4 and exact are copied from the paper; code/check_large.py checks that they are the right ones.
"""
import sys
from math import isqrt

EXACT = [207, 216, 225, 240, 243, 252, 256, 279, 288, 300, 320, 324, 336, 351,
         360, 363, 375, 378, 384, 392, 400, 432, 441, 480, 486, 504, 512, 540,
         567, 576, 600, 640, 648, 675, 720, 756, 768, 784, 800, 864, 900, 972,
         1024, 1080, 1152, 1296, 1440, 1536, 1728, 2304, 2592, 3072, 3456, 4608, 6912]
IDX4 = [216, 225, 256, 288, 320, 384, 432, 512, 576, 648, 675,
        768, 864, 1024, 1152, 1296, 1536, 1728, 2304, 2592, 3456, 4608]


def factor(n):
    f, p = {}, 2
    while p * p <= n:
        while n % p == 0:
            f[p] = f.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def psi(n):
    r = n
    for p in factor(n):
        r = r // p * (p + 1)
    return r


def least_prime_not_dividing(n):
    ell = 2
    while n % ell == 0 or any(ell % q == 0 for q in range(2, isqrt(ell) + 1)):
        ell += 1
    return ell


def squarefree(n):
    return all(e == 1 for e in factor(n).values())


def x0_levels():
    # N > 500 not in IDX4, and either in EXACT or 7 psi(N) <= 1600 (ell + 1); the second condition fails for
    # N >= 2743 by the lemma on gonality, and max(EXACT) = 6912.
    return [n for n in range(501, max(EXACT) + 1) if n not in IDX4
            and (n in EXACT or 7 * psi(n) <= 1600 * (least_prime_not_dividing(n) + 1))]


assert len(EXACT) == 55 and len(IDX4) == 22 and len(x0_levels()) == 113

LISTS = {
    "small": lambda: list(range(13, 201)),
    "sqfree": lambda: [n for n in range(201, 501) if squarefree(n)],
    "nonsqfree": lambda: [n for n in range(201, 501) if not squarefree(n)],
    "idx4": lambda: IDX4,
    "x0": x0_levels,
    "exact": lambda: EXACT,
}

if __name__ == "__main__":
    print(" ".join(map(str, LISTS[sys.argv[1]]())))
