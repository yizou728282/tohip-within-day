import sys, glob
from pathlib import Path
import numpy as np, pandas as pd
rng = np.random.default_rng(42)
hourly, rhdir, sitesf, outf = sys.argv[1:5]
sites = pd.read_csv(sitesf).set_index("site_id")
B = 300
def demean(df, cols, by): return df[cols] - df.groupby(by)[cols].transform("mean")
def coef(y, X):
    X = np.column_stack([np.ones(len(y)), X]); b, *_ = np.linalg.lstsq(X, y, rcond=None); return b[1]
def fit(d, ycol, xcol, ctrls, fe):
    key = d["lhour"].astype(str) if fe == "h" else d["month"].astype(str) + "_" + d["lhour"].astype(str)
    H = pd.get_dummies(key, prefix="k", drop_first=True).astype(float)
    D = pd.concat([d[[ycol, xcol, "day"] + ctrls], H], axis=1)
    cols = [ycol, xcol] + ctrls + list(H.columns)
    W = demean(D, cols, "day"); y = W[ycol].values; X = W[[xcol] + ctrls + list(H.columns)].values
    est = coef(y, X)
    days = d["day"].values; ud = np.unique(days); idx = {k: np.where(days == k)[0] for k in ud}
    bs = [coef(y[ii], X[ii]) for ii in (np.concatenate([idx[k] for k in rng.choice(ud, len(ud))]) for _ in range(B))]
    lo, hi = np.percentile(bs, [2.5, 97.5]); return est*10, lo*10, hi*10
specs = [("h_PW", "RH", ["PW"], "h"), ("mh_PW", "RH", ["PW"], "mh"), ("mh_PW_AOD", "RH", ["PW", "lnAOD"], "mh"),
         ("PLACEBO_lag7", "RH_lag7", ["PW"], "mh")]
rows = []
for f in sorted(glob.glob(str(Path(hourly) / "*.csv"))):
    d = pd.read_csv(f, parse_dates=["hour"]); site = d.site_id.iloc[0]; lon = sites.loc[site, "lon"]
    full = pd.read_csv(Path(rhdir) / f"{site}.csv", parse_dates=["hour"]).set_index("hour").RH
    d["RH_lag7"] = (d.hour - pd.Timedelta(days=7)).map(full)
    d["AOD"] = d.AOD500.fillna(d.AOD440)
    d = d[(d.AOD > 0) & d.RH.between(1, 99) & d.PW.notna() & d.RH_lag7.notna()].copy()
    d["lnAOD"] = np.log(d.AOD)
    lt = d.hour + pd.to_timedelta(lon / 15, unit="h")
    d["lhour"] = lt.dt.hour; d["day"] = lt.dt.date.astype(str); d["month"] = lt.dt.month
    for ycol in ["AE", "alpha_f", "FMF"]:
        e = d.dropna(subset=[ycol]).copy(); e = e[e.day.map(e.groupby("day").size()) >= 3]
        if e.day.nunique() < 50: continue
        r = dict(site=site, regime=e.regime.iloc[0], var=ycol, n_days=e.day.nunique(),
                 sdRHw=float(demean(e, ["RH"], "day").RH.std()),
                 corr_lag=float(np.corrcoef(demean(e, ["RH"], "day").RH, demean(e, ["RH_lag7"], "day").RH_lag7)[0, 1]))
        for name, x, ctl, fe in specs:
            r[name], r[name+"_lo"], r[name+"_hi"] = fit(e, ycol, x, ctl, fe)
        dm = e.groupby("day").RH.transform("mean")
        for tag, m in [("hiRH", dm >= 75), ("loRH", dm < 60)]:
            s = e[m]
            if s.day.nunique() >= 30: r[tag], r[tag+"_lo"], r[tag+"_hi"] = fit(s, ycol, "RH", ["PW"], "mh")
        rows.append(r)
R = pd.DataFrame(rows); R.to_csv(outf, index=False)
def fmt(r, k):
    v = r.get(k)
    if v is None or pd.isna(v): return "."
    return f"{v:+.3f}" + ("*" if (r[k+'_lo'] > 0 or r[k+'_hi'] < 0) else " ")
cols = ["h_PW", "mh_PW", "mh_PW_AOD", "hiRH", "loRH", "PLACEBO_lag7"]
for v in ["AE", "alpha_f", "FMF"]:
    s = R[R["var"] == v].sort_values(["regime", "site"])
    print(f"\n===== {v}: within-day slope per +10% RH; * = 95% day-cluster CI excludes 0 =====")
    print(f"{'site':18s}{'regime':7s}{'days':>6s}{'sdRHw':>6s}" + "".join(f"{c:>13s}" for c in cols))
    for _, r in s.iterrows():
        print(f"{r.site:18s}{r.regime:7s}{r.n_days:6d}{r.sdRHw:6.1f}" + "".join(f"{fmt(r, c):>13s}" for c in cols))
    for c in ["mh_PW", "PLACEBO_lag7"]:
        print(f"  {c}: sig-neg {(s[c+'_hi']<0).sum()}/{len(s)}, sig-pos {(s[c+'_lo']>0).sum()}/{len(s)}; median |slope| {s[c].abs().median():.4f}")
    print("  regime median mh_PW:", s.groupby("regime").mh_PW.median().round(4).to_dict())
print("\nmedian within-day corr(RH dev, RH_lag7 dev):", round(R.corr_lag.median(), 3))
