"""Preregistered confirmatory analysis (TOHIP within-day humidity response, tag prereg-v1).

usage: python confirm_analysis.py SITE_LIST HOURLY_DIR RH_FULL_DIR OUT_DIR
  SITE_LIST   prereg/confirmation_sites.csv (or discovery/sites.csv for a dry run)
  HOURLY_DIR  per-site hourly tables written by fetch_hourly.py
  RH_FULL_DIR per-site full hourly RH series written by fetch_rh.py
Nothing here may be tuned after the preregistration timestamp; any change is a logged deviation.
"""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats, sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sitelist import read_sites

B, SEED, MIN_DAYS, MIN_H, MIN_SUB_DAYS, MIN_K = 1000, 42, 200, 3, 30, 8
VARS = ["AE", "alpha_f", "FMF"]
EQ = {"AE": 0.010, "alpha_f": 0.010, "FMF": 0.004}          # equivalence bounds, per +10 % RH
rng = np.random.default_rng(SEED)


# ---------------------------------------------------------------- site level
def load_site(f, rh_dir, lon):
    d = pd.read_csv(f, parse_dates=["hour"])
    full = pd.read_csv(Path(rh_dir) / f.name, parse_dates=["hour"]).set_index("hour").RH
    d["RH_lag7"] = (d.hour - pd.Timedelta(days=7)).map(full)
    d["RH_prev1h"] = (d.hour - pd.Timedelta(hours=1)).map(full)
    d["AOD"] = d.AOD500.fillna(d.AOD440)
    d = d[(d.AOD > 0) & d.RH.between(1, 99) & d.PW.notna()].copy()
    d["lnAOD"] = np.log(d.AOD)
    d["x_k"] = d.RH / (100 - d.RH)                        # kappa-Koehler water-activity term
    d["rise"] = (d.RH > d.RH_prev1h).astype(float)
    d["RHxrise"] = d.RH * d.rise
    lt = d.hour + pd.to_timedelta(lon / 15, unit="h")      # local solar time
    d["lhour"], d["day"], d["month"] = lt.dt.hour, lt.dt.date.astype(str), lt.dt.month
    d["dayRH"] = d.groupby("day").RH.transform("mean")
    d["aod_jump"] = d.AOD > 2 * d.groupby("day").AOD.transform("median")
    return d


def keep_days(e, min_h):
    return e[e.day.map(e.groupby("day").size()) >= min_h]


def demean_by(M, g):
    G = sparse.csr_matrix((np.ones(len(g)), (np.arange(len(g)), g)))
    cnt = np.asarray(G.sum(0)).ravel()
    return M - G @ ((G.T @ M) / cnt[:, None])


def fit(e, ycol, xcols, ctrls, scale):
    """Day FE + month x local-hour FE + controls; slope(s) of xcols.

    Frisch-Waugh-Lovell: y and x are day-demeaned and residualised on (controls + FE dummies);
    the day-clustered bootstrap (B draws of days with replacement) then resamples the final
    regression on those residuals, holding the control coefficients fixed.
    """
    e = e.dropna(subset=[ycol] + xcols + ctrls)
    if e.day.nunique() < 2:
        return None
    g = pd.factorize(e.day)[0]
    key = pd.factorize(e.month.astype(str) + "_" + e.lhour.astype(str))[0]
    H = np.zeros((len(e), key.max() + 1)); H[np.arange(len(e)), key] = 1; H = H[:, 1:]
    Z = demean_by(np.column_stack([e[ctrls].values, H]) if ctrls else H, g)
    YX = demean_by(e[[ycol] + xcols].values.astype(float), g)
    YX = YX - Z @ np.linalg.lstsq(Z, YX, rcond=None)[0]
    ry, RX = YX[:, 0], YX[:, 1:]
    k = RX.shape[1]; nd = g.max() + 1
    XtX = np.zeros((nd, k, k)); Xty = np.zeros((nd, k))
    np.add.at(XtX, g, RX[:, :, None] * RX[:, None, :]); np.add.at(Xty, g, RX * ry[:, None])
    beta = np.linalg.solve(XtX.sum(0), Xty.sum(0))
    C = rng.multinomial(nd, np.full(nd, 1 / nd), size=B).astype(float)
    bb = np.linalg.solve(np.einsum("bd,dij->bij", C, XtX), (C @ Xty)[:, :, None])[:, :, 0]
    return dict(beta=beta * scale, se=bb.std(0, ddof=1) * scale,
                lo=np.percentile(bb, 2.5, 0) * scale, hi=np.percentile(bb, 97.5, 0) * scale,
                n_days=nd, n_hours=len(e))


def site_rows(d, site, cls):
    out = []

    def add(var, spec, r, idx=0):
        if r is None:
            return
        out.append(dict(site=site, cls=cls, var=var, spec=spec, beta=r["beta"][idx], se=r["se"][idx],
                        lo=r["lo"][idx], hi=r["hi"][idx], n_days=r["n_days"], n_hours=r["n_hours"]))

    def sub(var, spec, s, x=("RH",), ctrls=("PW",), scale=10):
        if s.day.nunique() >= MIN_SUB_DAYS:
            add(var, spec, fit(s, var, list(x), list(ctrls), scale))

    for var in VARS:
        e = keep_days(d.dropna(subset=[var]), MIN_H)
        if e.day.nunique() < MIN_DAYS:
            out.append(dict(site=site, cls=cls, var=var, spec="EXCLUDED_lt_min_days", n_days=e.day.nunique()))
            continue
        add(var, "primary", fit(e, var, ["RH"], ["PW"], 10))
        add(var, "sens_lnAOD", fit(e, var, ["RH"], ["PW", "lnAOD"], 10))
        add(var, "sens_xk", fit(e, var, ["x_k"], ["PW"], 1))
        sub(var, "sens_sza70", keep_days(e[e.SZA <= 70], MIN_H))
        sub(var, "sens_5h", keep_days(e, 5))
        sub(var, "sens_aodjump", keep_days(e[~e.aod_jump], MIN_H))
        sub(var, "hiRH", e[e.dayRH >= 75])
        sub(var, "loRH", e[e.dayRH < 60])
        add(var, "placebo_lag7", fit(e.dropna(subset=["RH_lag7"]), var, ["RH_lag7"], ["PW"], 10))
        r = fit(e.dropna(subset=["RH_prev1h"]), var, ["RH", "RHxrise"], ["PW"], 10)
        add(var, "expl_fall", r, 0)                 # slope when RH is falling
        add(var, "expl_rise_minus_fall", r, 1)      # extra slope when RH is rising
    return out


def between_day(d, site, cls):
    """S5: daily means; Spearman of residualised Y vs RH under three control sets."""
    rows = []
    D = d.groupby("day")[["AE", "alpha_f", "FMF", "RH", "PW", "lnAOD"]].mean()
    for var in VARS:
        x = D.dropna(subset=[var, "RH", "PW", "lnAOD"])
        if len(x) < 50:
            continue

        def res(y, cols):
            M = np.column_stack([np.ones(len(x))] + [x[c].values for c in cols])
            return y - M @ np.linalg.lstsq(M, y, rcond=None)[0]

        rows.append(dict(site=site, cls=cls, var=var, n_days=len(x),
                         rho_raw=stats.spearmanr(x[var], x.RH)[0],
                         rho_PW=stats.spearmanr(res(x[var].values, ["PW"]), res(x.RH.values, ["PW"]))[0],
                         rho_PW_lnAOD=stats.spearmanr(res(x[var].values, ["PW", "lnAOD"]), res(x.RH.values, ["PW"]))[0]))
    return rows


# ---------------------------------------------------------------- meta analysis
def reml_tau2(y, v, it=500):
    """REML between-site variance by fixed-point iteration (Viechtbauer 2005)."""
    t2 = max(0.0, np.var(y, ddof=1) - v.mean())
    for _ in range(it):
        w = 1 / (v + t2); mu = (w * y).sum() / w.sum()
        new = max(0.0, (w**2 * ((y - mu)**2 - v)).sum() / (w**2).sum() + 1 / w.sum())
        if abs(new - t2) < 1e-14:
            break
        t2 = new
    return t2


def meta(df):
    df = df.dropna(subset=["beta", "se"]); df = df[df.se > 0]
    k = len(df)
    if k < 2:
        return dict(k=k)
    y, v = df.beta.values, df.se.values**2
    t2 = reml_tau2(y, v); w = 1 / (v + t2); mu = (w * y).sum() / w.sum()
    q = (w * (y - mu)**2).sum() / (k - 1); se = np.sqrt(max(q, 1e-300) / w.sum())     # HKSJ
    wf = 1 / v; Q = (wf * (y - (wf * y).sum() / wf.sum())**2).sum()
    tq = lambda a: stats.t.ppf(a, k - 1)
    return dict(k=k, mu=mu, se=se, tau2=t2, I2=max(0.0, (Q - (k - 1)) / Q) if Q > 0 else 0.0,
                ci95=(mu - tq(.975) * se, mu + tq(.975) * se), ci90=(mu - tq(.95) * se, mu + tq(.95) * se),
                p_less=stats.t.cdf(mu / se, k - 1), p_greater=stats.t.sf(mu / se, k - 1),
                n_sig=int(((df.lo > 0) | (df.hi < 0)).sum()), descriptive_only=bool(k < MIN_K))


def jsonable(o):
    if hasattr(o, "item"):
        return o.item()
    return str(o)


def main():
    sites, hourly, rhdir, out = read_sites(sys.argv[1]), Path(sys.argv[2]), sys.argv[3], Path(sys.argv[4])
    out.mkdir(parents=True, exist_ok=True)
    srows, brows = [], []
    for s in sites.itertuples():
        f = hourly / f"{s.site_id}.csv"
        if not f.exists():
            srows.append(dict(site=s.site_id, cls=s.cls, spec="MISSING_DATA")); continue
        d = load_site(f, rhdir, s.lon)
        srows += site_rows(d, s.site_id, s.cls); brows += between_day(d, s.site_id, s.cls)
        print("site done:", s.site_id, flush=True)
    S = pd.DataFrame(srows); S.to_csv(out / "site_results.csv", index=False)
    BD = pd.DataFrame(brows); BD.to_csv(out / "between_day.csv", index=False)
    ok = S[S.beta.notna()] if "beta" in S else S.iloc[0:0]
    M = {key: meta(g) for key, g in ok.groupby(["cls", "var", "spec"])}
    pd.DataFrame([dict(cls=c, var=v, spec=s, **{k: (str(x) if isinstance(x, tuple) else x) for k, x in m.items()})
                  for (c, v, s), m in M.items()]).to_csv(out / "meta_results.csv", index=False)
    get = lambda c, v, s="primary": M.get((c, v, s), {"k": 0})
    eqv = lambda m, b: bool(m.get("k", 0) >= 2 and -b < m["ci90"][0] and m["ci90"][1] < b)
    H = {}
    # ---- primary hypotheses (one-sided, Holm over 4)
    for h, (m, side) in {"H1a": (get("pollution", "alpha_f"), "p_less"),
                         "H1b": (get("pollution", "FMF"), "p_greater"),
                         "H2": (get("marine", "AE"), "p_less")}.items():
        H[h] = dict(k=m.get("k"), mu=m.get("mu"), se=m.get("se"), ci95=m.get("ci95"), p=m.get(side),
                    descriptive_only=m.get("descriptive_only", True))
    a, b = get("pollution", "AE"), get("marine", "AE")
    if a.get("k", 0) >= 2 and b.get("k", 0) >= 2:
        diff = a["mu"] - b["mu"]; se = float(np.hypot(a["se"], b["se"])); df_ = min(a["k"], b["k"]) - 1
        H["H3"] = dict(diff=diff, se=se, df=df_, p=stats.t.sf(diff / se, df_),
                       descriptive_only=a["descriptive_only"] or b["descriptive_only"])
    order = sorted([h for h in ["H1a", "H1b", "H2", "H3"] if H.get(h, {}).get("p") is not None], key=lambda h: H[h]["p"])
    running = 0.0
    for i, h in enumerate(order):
        running = max(running, min(1.0, (len(order) - i) * H[h]["p"]))
        H[h]["p_holm"] = running
        H[h]["supported"] = bool(running < 0.05 and not H[h]["descriptive_only"])
    # ---- secondary
    H["S1"] = dict(k=a.get("k"), mu=a.get("mu"), ci95=a.get("ci95"), p=a.get("p_greater"))
    m = get("marine", "FMF"); H["S2"] = dict(k=m.get("k"), mu=m.get("mu"), ci95=m.get("ci95"), p=m.get("p_less"))
    H["S3"] = {f"{c}_{v}": dict(k=get(c, v).get("k"), mu=get(c, v).get("mu"), ci90=get(c, v).get("ci90"), bound=EQ[v],
                                equivalent=eqv(get(c, v), EQ[v]), descriptive_only=get(c, v).get("descriptive_only", True))
               for c in ["dust", "burning"] for v in VARS}
    # S3 sensitivity: "core" dust sites only (climatological median AE < 0.6), descriptive
    raw = pd.read_csv(sys.argv[1])
    core = set(raw.loc[(raw["class"] == "dust") & (raw["medAE"] < 0.6), "site"]) \
        if {"class", "medAE", "site"} <= set(raw.columns) else set()
    H["S3_dust_core"] = {"sites": sorted(core)}
    for v in VARS:
        m = meta(ok[ok.site.isin(core) & (ok["var"] == v) & (ok.spec == "primary")])
        H["S3_dust_core"][v] = dict(k=m.get("k"), mu=m.get("mu"), ci90=m.get("ci90"), bound=EQ[v],
                                    equivalent=eqv(m, EQ[v]), descriptive_only=True)
    H["S4"] = {}
    for c, v, sign in [("pollution", "alpha_f", -1), ("pollution", "FMF", +1), ("marine", "AE", -1)]:
        sel = lambda spec: ok[(ok.cls == c) & (ok["var"] == v) & (ok.spec == spec)].set_index("site")[["beta", "se"]]
        j = sel("hiRH").join(sel("loRH"), lsuffix="_hi", rsuffix="_lo", how="inner")
        D = pd.DataFrame(dict(beta=sign * (j.beta_hi - j.beta_lo), se=np.hypot(j.se_hi, j.se_lo), lo=np.nan, hi=np.nan))
        mm = meta(D); H["S4"][f"{c}_{v}"] = dict(k=mm.get("k"), mu=mm.get("mu"), ci95=mm.get("ci95"), p=mm.get("p_greater"))
    H["S5"] = {}
    cols = ["rho_raw", "rho_PW", "rho_PW_lnAOD"]
    for v in VARS:
        x = BD[BD["var"] == v] if len(BD) else BD
        if len(x) < 3:
            continue
        med = lambda xx: [float(np.median(np.abs(xx[c]))) for c in cols]
        boots = np.array([med(x.sample(len(x), replace=True, random_state=i)) for i in range(B)])
        r0, r1, r2 = med(x)
        H["S5"][v] = dict(n_sites=len(x), med_abs_rho=dict(zip(cols, (r0, r1, r2))),
                          ci95={c: tuple(np.percentile(boots[:, i], [2.5, 97.5])) for i, c in enumerate(cols)},
                          PW_reduces_le_20pct=bool(r1 >= 0.8 * r0), lnAOD_reduces_ge_40pct=bool(r2 <= 0.6 * r1))
    H["S6"] = {}
    for v in VARS:
        pl = ok[(ok["var"] == v) & (ok.spec == "placebo_lag7")]
        frac = float(((pl.lo > 0) | (pl.hi < 0)).mean()) if len(pl) else None
        H["S6"][v] = dict(frac_sites_sig=frac, frac_ok=bool(frac is not None and frac <= 0.10),
                          by_class={c: dict(k=get(c, v, "placebo_lag7").get("k"), ci90=get(c, v, "placebo_lag7").get("ci90"),
                                            equivalent=eqv(get(c, v, "placebo_lag7"), EQ[v])) for c in sorted(pl.cls.unique())})
    json.dump(H, open(out / "hypotheses.json", "w"), indent=2, default=jsonable)
    print(json.dumps(H, indent=1, default=jsonable))


if __name__ == "__main__":
    main()
