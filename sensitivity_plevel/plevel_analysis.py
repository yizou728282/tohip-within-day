"""DEVIATION / SUBSTITUTE for preregistered sensitivity #6 (layer RH instead of 2 m RH).

Uses the frozen load_site / keep_days / fit / meta functions of confirm_analysis.py (prereg-v2) unchanged.
For pollution and marine sites, the primary model is refitted with four RH predictors on one common set of
hours (2021-03 to 2024-12, all four RH values present and both pressure levels at least 100 m above the
site, so the comparison is not confounded by sampling):
  RH2m_om  - Open-Meteo archive 2 m RH (the preregistered predictor)
  RH2m_gfs - GFS 2 m RH
  RH925    - GFS 925 hPa RH
  RH850    - GFS 850 hPa RH
Choices made before running: the 100 m clearance rule, the common-hours design, MIN_DAYS and MIN_H unchanged.

usage: python plevel_analysis.py FROZEN_SCRIPTS_DIR SITE_LIST HOURLY_DIR RH_FULL_DIR GFS_DIR OUT_DIR
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

sys.path.insert(0, sys.argv[1])
import confirm_analysis as ca

VARIANTS = ["RH2m_om", "RH2m_gfs", "RH925", "RH850"]
TESTS = {"H1a": ("pollution", "alpha_f"), "H1b": ("pollution", "FMF"), "H2": ("marine", "AE")}


def main():
    sites = ca.read_sites(sys.argv[2])
    hourly, rhdir, gfs, out = Path(sys.argv[3]), sys.argv[4], Path(sys.argv[5]), Path(sys.argv[6])
    out.mkdir(parents=True, exist_ok=True)
    elev = pd.read_csv(sys.argv[2]).set_index("site").elev
    rows, cover = [], []
    for s in sites[sites.cls.isin(["pollution", "marine"])].itertuples():
        g = gfs / f"{s.site_id}.csv"
        if not g.exists():
            cover.append(dict(site=s.site_id, cls=s.cls, status="NO_GFS_FILE")); continue
        d = ca.load_site(hourly / f"{s.site_id}.csv", rhdir, s.lon)
        G = pd.read_csv(g, parse_dates=["hour"])
        d = d.merge(G, on="hour", how="inner").rename(columns={"RH": "RH2m_om"})
        ok = (d.z925 - elev[s.site_id] >= 100) & (d.z850 - elev[s.site_id] >= 100)
        for v in VARIANTS:
            ok &= d[v].between(1, 99)
        d = d[ok].copy()
        for var in ca.VARS:
            e = ca.keep_days(d.dropna(subset=[var]), ca.MIN_H)
            nd = e.day.nunique()
            cover.append(dict(site=s.site_id, cls=s.cls, var=var, n_days=nd, status="ok" if nd >= ca.MIN_DAYS else "lt_min_days"))
            if nd < ca.MIN_DAYS:
                continue
            for v in VARIANTS:
                e2 = e.assign(RH=e[v])
                r = ca.fit(e2, var, ["RH"], ["PW"], 10)
                sd = e2.RH.sub(e2.groupby("day").RH.transform("mean")).std()
                rows.append(dict(site=s.site_id, cls=s.cls, var=var, rh=v, beta=r["beta"][0], se=r["se"][0],
                                 lo=r["lo"][0], hi=r["hi"][0], n_days=r["n_days"], n_hours=r["n_hours"], sd_within=sd))
        print("done", s.site_id, flush=True)
    R = pd.DataFrame(rows); R.to_csv(out / "plevel_site_results.csv", index=False)
    pd.DataFrame(cover).to_csv(out / "plevel_coverage.csv", index=False)
    M = []
    for (c, var, v), g in R.groupby(["cls", "var", "rh"]):
        m = ca.meta(g)
        M.append(dict(cls=c, var=var, rh=v, k=m["k"], mu=m["mu"], se=m["se"], lo=m["ci95"][0], hi=m["ci95"][1],
                      I2=m["I2"], n_sig=m["n_sig"], p_less=m["p_less"], p_greater=m["p_greater"],
                      median_sd_within=g.sd_within.median()))
    M = pd.DataFrame(M); M.to_csv(out / "plevel_meta.csv", index=False)
    H = []
    for v in VARIANTS:
        get = lambda c, var: M[(M.cls == c) & (M["var"] == var) & (M.rh == v)].iloc[0]
        for h, (c, var) in TESTS.items():
            m = get(c, var)
            H.append(dict(rh=v, test=h, k=m.k, mu=m.mu, lo=m.lo, hi=m.hi, p=m.p_less if h in ("H1a", "H2") else m.p_greater))
        a, b = get("pollution", "AE"), get("marine", "AE")
        diff, se = a.mu - b.mu, float(np.hypot(a.se, b.se)); df = min(a.k, b.k) - 1
        H.append(dict(rh=v, test="H3", k=f"{a.k}/{b.k}", mu=diff, lo=np.nan, hi=np.nan, p=stats.t.sf(diff / se, df)))
    H = pd.DataFrame(H); H.to_csv(out / "plevel_hypotheses.csv", index=False)
    print(M.round(4).to_string(index=False)); print(H.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
