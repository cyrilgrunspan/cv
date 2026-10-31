#!/usr/bin/env python3
"""
Cagnotte critique d'un état de la course (Grunspan & Pérez-Marco,
« Is the game worth the candle? Optimal stopping in a race with sunk investment »).

Usage :  python3 cagnotte.py a h q z [tolérance]

Calcule un encadrement [ℓ, u] de la cagnotte critique v*_z(a, h) (Définition 3.17) :
la cagnotte au-dessus de laquelle, à l'état (a, h), continuer (Wait) vaut mieux
qu'abandonner (Quit). Pour (a, h) = (0, 0) on retrouve v*(q, z) (Lemme 3.18(ii)).

Méthode (Algorithme 1 adapté à un état quelconque) : pour une cagnotte v, la valeur
V_v de la course est encadrée par deux inductions rétrogrades sur le triangle
a + h ≤ N̄, qui diffèrent seulement sur la dernière diagonale (abandon forcé pour
la borne inférieure, borne P*_z(a, h)·(a + max(v, 0)) de la valeur pour la borne supérieure, analogue de U*) ;
on continue à (a, h) ssi q·V_v(a+1, h) + p·V_v(a, h+1) − q > 0 (Lemme 3.13 et 3.19(iii)).
Une bissection en v localise le changement de signe de chaque borne, et N̄ est
augmenté jusqu'à ce que l'encadrement ait la largeur demandée.
"""
import sys


def win(a, h, z):
    return a > h and h >= z


def pstar(a, h, q, z):
    """P*_z(a, h) : probabilité de gagner d'un challenger qui n'abandonne jamais (formule (3.2))."""
    from math import comb
    p = 1.0 - q
    lam = q / p
    if win(a, h, z):
        return 1.0
    if h >= z:
        return lam ** (h - a + 1)
    m = z - h
    s = 0.0
    sb = 0.0
    for k in range(0, z - a + 1):
        b = comb(k + m - 1, k) * q ** k * p ** m
        sb += b
        s += b * lam ** (z + 1 - a - k)
    return s + (1.0 - sb)


def wait_value(a0, h0, q, z, v, N, optimistic):
    """q·V(a0+1,h0) + p·V(a0,h0+1) − q, V étant la valeur tronquée (Définition 3.14
    adaptée : gains au lieu de la forme linéarisée) sur le triangle a + h ≤ N."""
    p = 1.0 - q
    # dernière diagonale a + h = N
    V = [0.0] * (N + 1)
    for a in range(N + 1):
        h = N - a
        if win(a, h, z):
            V[a] = max(a + v, 0.0)
        elif optimistic:
            V[a] = pstar(a, h, q, z) * (a + max(v, 0.0))
        else:
            V[a] = 0.0
    n0 = a0 + h0 + 1
    for n in range(N - 1, n0 - 1, -1):
        W = [0.0] * (n + 1)
        for a in range(n + 1):
            h = n - a
            if win(a, h, z):
                W[a] = max(a + v, q * V[a + 1] + p * V[a] - q, 0.0)
            else:
                W[a] = max(0.0, q * V[a + 1] + p * V[a] - q)
        V = W
    return q * V[a0 + 1] + p * V[a0] - q


def critical_bound(a, h, q, z, N, optimistic, tol):
    """Changement de signe en v de wait_value, par bissection (géométrique si v > 0)."""
    f = lambda v: wait_value(a, h, q, z, v, N, optimistic)
    lo, hi = -float(a), 1.0            # v*_z(x) ≥ −a (Lemme 3.18(ii))
    if f(lo) > 0:
        return lo, lo
    while f(hi) <= 0:
        hi *= 4.0
        if hi > 1e300:
            return float("inf"), float("inf")
    # invariant : f(lo) ≤ 0 < f(hi)
    for _ in range(200):
        if hi > 0 and lo > 0:
            if hi / lo <= 1 + tol / 4:
                break
            m = (lo * hi) ** 0.5
        else:
            if hi - lo <= tol / 4:
                break
            m = 0.5 * (lo + hi)
        if f(m) > 0:
            hi = m
        else:
            lo = m
    return lo, hi


def critical_prize(a, h, q, z, tol=1e-6, N=None, verbose=False):
    """Encadrement [ℓ, u] de v*_z(a, h), de largeur relative (ou absolue si ℓ ≤ 0) ≤ tol."""
    if not (0 < q < 0.5):
        raise ValueError("q doit être dans (0, 1/2)")
    if z < 1 or a < 0 or h < 0:
        raise ValueError("z ≥ 1, a ≥ 0, h ≥ 0")
    if win(a, h, z):
        raise ValueError("Win est disponible à l'état (a, h) : la cagnotte critique n'y est pas définie")
    if N is None:
        N = max(2 * z + 1, a + h + 2, 8)
    while True:
        _, u = critical_bound(a, h, q, z, N, False, tol)  # abandon forcé : valeur ≤ vraie, seuil ≥ v* → borne sup
        l, _ = critical_bound(a, h, q, z, N, True, tol)   # borne optimiste : valeur ≥ vraie, seuil ≤ v* → borne inf
        if verbose:
            print(f"N̄ = {N:4d} : [{l:.10g}, {u:.10g}]", file=sys.stderr)
        width = (u / l - 1) if l > 0 else (u - l)
        if width <= tol or N > 20000:
            return l, u, N
        N = int(N * 1.5) + 1


if __name__ == "__main__":
    if len(sys.argv) < 5:
        print(__doc__)
        sys.exit(1)
    a, h, q, z = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
    tol = float(sys.argv[5]) if len(sys.argv) > 5 else 1e-6
    l, u, N = critical_prize(a, h, q, z, tol, verbose=True)
    print(f"v*_{z}(({a}, {h})) ∈ [{l:.10g}, {u:.10g}]   (q = {q}, troncature N̄ = {N})")
