"""Exact check of the reservoir construction for G = q U[c-H, c+H] + (1-q) T, T any law on [0, R].

Thin positions (K of them, K in {floor m, floor m + 1}, m = n(1-q)) take T-values Q((j+W)/K), j = 0..K-1, W ~ U[0,1).
Fat positions take c + Z, Z with density (1 + beta z/H)/(2H), beta = 3C/((n-K)H), C = n mu + m_E - S - Y - (n-K)c,
where Y ~ E is one extra independent variable on an interval of length L (L = 0: none).
Checked in rationals for discrete T (the worst case for the quantile grid):
  |S - K mT| <= range(T); E_W[S] = K mT; |beta| <= 1; n - K >= 3; E[C] = 0 (fat marginal exactly uniform).
"""
import random
from fractions import Fraction as Fr
from math import floor

random.seed(20261029)


def grid_pieces(atoms, k):
    """Pieces of W in [0,1) on which S(W) = sum_j Q((j+W)/k) is constant: list of (length, S)."""
    cum, cdf = Fr(0), []
    for x, p in atoms:
        cum += p
        cdf.append((cum, x))
    breaks = {Fr(0), Fr(1)}
    for u, _ in cdf[:-1]:
        for j in range(k):
            w = u * k - j
            if 0 < w < 1:
                breaks.add(w)
    b = sorted(breaks)
    out = []
    for lo, hi in zip(b, b[1:]):
        mid = (lo + hi) / 2
        s = Fr(0)
        for j in range(k):
            v = (j + mid) / k
            s += next(x for u, x in cdf if v < u)
        out.append((hi - lo, s))
    return out


def trial():
    R = Fr(random.randint(1, 20))
    na = random.randint(1, 5)
    xs = [Fr(random.randint(0, 40), 2) * R / 20 for _ in range(na)]
    ws = [random.randint(1, 9) for _ in range(na)]
    atoms = sorted((x, Fr(w, sum(ws))) for x, w in zip(xs, ws))
    rng_T = atoms[-1][0] - atoms[0][0]
    mT = sum(x * p for x, p in atoms)
    H = Fr(random.randint(1, 40), 40) * R / 2
    c = H + Fr(random.randint(0, 40), 40) * (R - 2 * H) if R > 2 * H else R / 2
    q = Fr(random.randint(1, 19), 20)
    mu = q * c + (1 - q) * mT
    # extra independent variable with a discrete law E on [0, L] (L = 0 means no extra variable)
    L = Fr(random.randint(0, 3)) * R / 2
    ne = random.randint(1, 4)
    ex = [Fr(random.randint(0, 8), 8) * L for _ in range(ne)]
    ew = [random.randint(1, 9) for _ in range(ne)]
    E = [(x, Fr(w, sum(ew))) for x, w in zip(ex, ew)]
    mE = sum(x * p for x, p in E)
    Nf = (3 * (2 * R + L) / H + 1) / q
    N = -(-Nf.numerator // Nf.denominator)
    bad = 0
    for n in (N, N + 1, N + 3):
        m = n * (1 - q)
        k0 = floor(m)
        ks = [(k0, (k0 + 1) - m), (k0 + 1, m - k0)]
        EC = Fr(0)
        for k, pk in ks:
            if pk == 0:
                continue
            if n - k < 3:
                bad += 1
            pieces = grid_pieces(atoms, k) if k > 0 else [(Fr(1), Fr(0))]
            ES = sum(l * s for l, s in pieces)
            if ES != k * mT:
                bad += 1
            for l, s in pieces:
                if abs(s - k * mT) > rng_T:
                    bad += 1
                for y, py in E:
                    C = n * mu + mE - s - y - (n - k) * c
                    if abs(C) > 2 * R + L:
                        bad += 1
                    beta = 3 * C / ((n - k) * H)
                    if abs(beta) > 1:
                        bad += 1
                    EC += pk * l * py * C
        if EC != 0:
            bad += 1
    return bad


if __name__ == "__main__":
    trials, fails = 1000, 0
    for _ in range(trials):
        fails += trial() > 0
    print("trials", trials, "failures", fails)
