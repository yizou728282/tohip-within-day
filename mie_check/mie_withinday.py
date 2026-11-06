"""EXPLORATORY, NOT PREREGISTERED: within-day slopes implied by the Mie + kappa-Koehler model.

The model alpha_f and FMF are evaluated at every observed hour's RH (dry state fixed per site, calibrated
as in mie_check.py), then regressed on RH within days (day-demeaned OLS, slope per +10 % RH). This mirrors
the observed estimator except for the month x hour fixed effects and the PW control. The ratio of
observed to model slopes estimates the attenuation; the FMF/alpha_f slope ratio does not depend on it.

usage: python mie_withinday.py SITE_LIST HOURLY_DIR SITE_RESULTS_CSV OUT_DIR
"""
import sys, itertools
from pathlib import Path
import numpy as np, pandas as pd
from scipy.optimize import brentq
sys.path.insert(0, str(Path(__file__).resolve().parent))
import mie_check as mc


def tables(af_med, fmf_med, rh_med, lnsig, m_dry, kappa, kappa_c, grid):
    f = lambda lg: mc.alpha_f(np.exp(lg), lnsig, m_dry, kappa, rh_med) - af_med
    rg = np.exp(brentq(f, np.log(0.01), np.log(0.4), xtol=1e-6))
    af = np.array([mc.alpha_f(rg, lnsig, m_dry, kappa, h) for h in grid])
    tf = np.array([mc.tau(rg, lnsig, m_dry, kappa, h, 0.5) for h in grid])
    rc_n = 2.5 * np.exp(-3 * 0.65**2); rcoarse = np.logspace(np.log10(0.1), np.log10(20.0), 300)
    key = ("wd", kappa_c)
    if key not in mc._COARSE:
        mc._COARSE[key] = np.array([mc.tau(rc_n, 0.65, 1.53 + 0.003j, kappa_c, h, 0.5, rcoarse) for h in grid])
    tc = mc._COARSE[key]
    i = np.argmin(np.abs(grid - rh_med))
    fmf = tf / (tf + tf[i] * (1 - fmf_med) / fmf_med / tc[i] * tc)
    return af, fmf


def wd_slope(day, x, y):
    d = pd.DataFrame(dict(day=day, x=x, y=y))
    d = d[d.day.map(d.groupby("day").size()) >= 3]
    xd = d.x - d.groupby("day").x.transform("mean"); yd = d.y - d.groupby("day").y.transform("mean")
    return float((xd * yd).sum() / (xd * xd).sum() * 10)


def main():
    sites, hourly, obs, out = pd.read_csv(sys.argv[1]), Path(sys.argv[2]), pd.read_csv(sys.argv[3]), Path(sys.argv[4])
    grid = np.arange(5.0, 97.5, 1.0)
    obs = obs[(obs.spec == "primary") & (obs.cls == "pollution")].pivot(index="site", columns="var", values="beta")
    rows = []
    for s in sites[sites["class"] == "pollution"].site:
        d = pd.read_csv(hourly / f"{s}.csv")
        d = d[d.RH.between(1, 99)].dropna(subset=["alpha_f", "FMF"]).copy()
        lt = pd.to_datetime(d.hour) + pd.to_timedelta(float(sites.set_index("site").lon[s]) / 15, unit="h")
        d["day"] = lt.dt.date.astype(str)
        rh = d.RH.clip(upper=96).values
        for kappa, kc in itertools.product([0.1, 0.2, 0.3, 0.4, 0.6], [0.0, 0.1]):
            af, fmf = tables(d.alpha_f.median(), d.FMF.median(), float(np.median(rh)), 0.5, 1.53 + 0.01j, kappa, kc, grid)
            rows.append(dict(site=s, kappa=kappa, kappa_c=kc,
                             model_alpha_f=wd_slope(d.day.values, rh, np.interp(rh, grid, af)),
                             model_FMF=wd_slope(d.day.values, rh, np.interp(rh, grid, fmf)),
                             obs_alpha_f=obs.loc[s, "alpha_f"], obs_FMF=obs.loc[s, "FMF"]))
        print("done", s, flush=True)
    R = pd.DataFrame(rows); R.to_csv(out / "mie_withinday_sites.csv", index=False)
    S = R.groupby(["kappa", "kappa_c"]).agg(model_alpha_f=("model_alpha_f", "mean"), model_FMF=("model_FMF", "mean"),
                                            obs_alpha_f=("obs_alpha_f", "mean"), obs_FMF=("obs_FMF", "mean")).reset_index()
    S["ratio_model"] = S.model_FMF / S.model_alpha_f
    S["atten_alpha_f"] = S.obs_alpha_f / S.model_alpha_f
    S["atten_FMF"] = S.obs_FMF / S.model_FMF
    S["r_sites_alpha_f"] = [R[(R.kappa == k) & (R.kappa_c == c)][["model_alpha_f", "obs_alpha_f"]].corr(method="spearman").iloc[0, 1]
                            for k, c in zip(S.kappa, S.kappa_c)]
    S["r_sites_FMF"] = [R[(R.kappa == k) & (R.kappa_c == c)][["model_FMF", "obs_FMF"]].corr(method="spearman").iloc[0, 1]
                        for k, c in zip(S.kappa, S.kappa_c)]
    S.to_csv(out / "mie_withinday_summary.csv", index=False)
    print("observed FMF/alpha_f ratio (site means):", round(R.obs_FMF.mean() / R.obs_alpha_f.mean(), 3))
    print(S.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
