"""EXPLORATORY, NOT PREREGISTERED: Mie + kappa-Koehler forward check of the H1 magnitude.

For each confirmation pollution site, a lognormal fine mode is calibrated so that the modelled
fine-mode Angstrom exponent (alpha_f, 500 nm) equals the site's observed median alpha_f at the site's
median RH. The dry particles are then grown with kappa-Koehler theory (no Kelvin term) and the expected
change of alpha_f and FMF per +10 % RH is averaged over the site's observed hourly RH distribution.

usage: python mie_check.py SITE_LIST HOURLY_DIR OUT_DIR
"""
import sys, itertools
from pathlib import Path
import numpy as np, pandas as pd
from scipy.optimize import brentq

M_WATER = 1.333 + 0j


def qext(m, x):
    """Extinction efficiency for homogeneous spheres (Bohren & Huffman BHMIE), vectorised over x."""
    x = np.atleast_1d(np.asarray(x, float))
    nstop = np.floor(x + 4 * x ** (1 / 3) + 2).astype(int)
    N = nstop.max()
    y = m * x
    nmx = int(max(N, np.abs(y).max())) + 16
    D = np.zeros((nmx + 1, x.size), complex)
    for n in range(nmx, 0, -1):
        D[n - 1] = n / y - 1 / (D[n] + n / y)
    psi0, psi1 = np.cos(x), np.sin(x)
    chi0, chi1 = -np.sin(x), np.cos(x)
    xi1 = psi1 - 1j * chi1
    q = np.zeros(x.size)
    with np.errstate(all="ignore"):
        for n in range(1, N + 1):
            psi = (2 * n - 1) / x * psi1 - psi0
            chi = (2 * n - 1) / x * chi1 - chi0
            xi = psi - 1j * chi
            a = ((D[n] / m + n / x) * psi - psi1) / ((D[n] / m + n / x) * xi - xi1)
            b = ((m * D[n] + n / x) * psi - psi1) / ((m * D[n] + n / x) * xi - xi1)
            q += np.where(n <= nstop, (2 * n + 1) * (a + b).real, 0.0)
            psi0, psi1, chi0, chi1 = psi1, psi, chi1, chi
            xi1 = psi1 - 1j * chi1
    return 2 / x**2 * q


def selftest():
    # Rayleigh limit and a textbook value (m = 1.5, x = 1: Qext = Qsca ~ 0.2153)
    m, x = 1.5 + 0j, 0.01
    ray = 8 / 3 * x**4 * abs((m**2 - 1) / (m**2 + 2)) ** 2
    assert abs(qext(m, x)[0] / ray - 1) < 1e-3, (qext(m, x), ray)
    assert abs(qext(m, 1.0)[0] - 0.2153) < 5e-4, qext(m, 1.0)
    assert abs(qext(1.33 + 0j, 2000.0)[0] - 2) < 0.02
    print("Mie self-test passed:", round(float(qext(m, 1.0)[0]), 5))


R_DRY = np.logspace(np.log10(0.005), np.log10(1.5), 300)          # um, dry radius grid
LNR = np.log(R_DRY)


def growth(rh, kappa):
    aw = rh / 100.0
    return (1 + kappa * aw / (1 - aw)) ** (1 / 3)


def tau(rg, lnsig, m_dry, kappa, rh, lam, r=R_DRY):
    """Extinction (arbitrary units) of a lognormal number distribution (median rg, ln sigma) at wavelength lam (um)."""
    n = np.exp(-0.5 * ((np.log(r) - np.log(rg)) / lnsig) ** 2)
    gf = growth(rh, kappa)
    m = (m_dry + (gf**3 - 1) * M_WATER) / gf**3                      # volume-weighted mixing rule
    rw = r * gf
    return np.trapezoid(qext(m, 2 * np.pi * rw / lam) * np.pi * rw**2 * n, np.log(r))


def alpha_f(rg, lnsig, m_dry, kappa, rh):
    t1, t2 = tau(rg, lnsig, m_dry, kappa, rh, 0.48), tau(rg, lnsig, m_dry, kappa, rh, 0.52)
    return -np.log(t2 / t1) / np.log(0.52 / 0.48)


_COARSE = {}


def site_check(rh_obs, af_med, fmf_med, lnsig, m_dry, kappa, kappa_c, rh_grid):
    rh_med = float(np.median(rh_obs))
    f = lambda lg: alpha_f(np.exp(lg), lnsig, m_dry, kappa, rh_med) - af_med
    try:
        rg = np.exp(brentq(f, np.log(0.01), np.log(0.4), xtol=1e-6))
    except ValueError:
        return None
    af = np.array([alpha_f(rg, lnsig, m_dry, kappa, h) for h in rh_grid])
    tf = np.array([tau(rg, lnsig, m_dry, kappa, h, 0.5) for h in rh_grid])
    # coarse mode: volume median radius 2.5 um, ln sigma 0.65, dry index 1.53+0.003i -> number median radius
    rc_n = 2.5 * np.exp(-3 * 0.65**2)
    rcoarse = np.logspace(np.log10(0.1), np.log10(20.0), 300)
    key = (kappa_c, rh_grid[0], rh_grid[-1], rh_grid.size)
    if key not in _COARSE:
        _COARSE[key] = np.array([tau(rc_n, 0.65, 1.53 + 0.003j, kappa_c, h, 0.5, rcoarse) for h in rh_grid])
    tc = _COARSE[key]
    i_med = np.argmin(np.abs(rh_grid - rh_med))
    scale = tf[i_med] * (1 - fmf_med) / fmf_med / tc[i_med]               # match observed median FMF
    fmf = tf / (tf + scale * tc)
    d_af = np.interp(rh_obs, rh_grid[1:-1], (af[2:] - af[:-2]) / (rh_grid[2:] - rh_grid[:-2])) * 10
    d_fmf = np.interp(rh_obs, rh_grid[1:-1], (fmf[2:] - fmf[:-2]) / (rh_grid[2:] - rh_grid[:-2])) * 10
    return dict(rg_dry=rg, rh_med=rh_med, d_alpha_f=float(d_af.mean()), d_FMF=float(d_fmf.mean()),
                af_at_40=float(np.interp(40, rh_grid, af)), af_at_85=float(np.interp(85, rh_grid, af)))


def main():
    selftest()
    sites, hourly, out = pd.read_csv(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    out.mkdir(parents=True, exist_ok=True)
    poll = sites[sites["class"] == "pollution"].site
    rh_grid = np.arange(5.0, 97.5, 1.0)
    rows = []
    grid = dict(kappa=[0.1, 0.2, 0.3, 0.4, 0.6], lnsig=[0.4, 0.5, 0.6],
                n_dry=[1.45, 1.53, 1.60], kappa_c=[0.0, 0.1])
    for s in poll:
        d = pd.read_csv(hourly / f"{s}.csv")
        d = d[d.RH.between(1, 99)].dropna(subset=["alpha_f", "FMF"])
        rh = d.RH.clip(upper=96).values
        for kappa, lnsig, n_dry, kc in itertools.product(*grid.values()):
            r = site_check(rh, d.alpha_f.median(), d.FMF.median(), lnsig, n_dry + 0.01j, kappa, kc, rh_grid)
            if r:
                rows.append(dict(site=s, kappa=kappa, lnsig=lnsig, n_dry=n_dry, kappa_c=kc, **r))
        print("done", s, flush=True)
    R = pd.DataFrame(rows); R.to_csv(out / "mie_site_results.csv", index=False)
    S = R.groupby(["kappa", "lnsig", "n_dry", "kappa_c"])[["d_alpha_f", "d_FMF", "rg_dry"]].mean().reset_index()
    S.to_csv(out / "mie_summary.csv", index=False)
    base = S[(S.lnsig == 0.5) & (S.n_dry == 1.53) & (S.kappa_c == 0.0)]
    print(base.round(4).to_string(index=False))
    print(S.groupby("kappa")[["d_alpha_f", "d_FMF"]].agg(["min", "max"]).round(4))


if __name__ == "__main__":
    main()
