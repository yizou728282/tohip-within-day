"""PREREGISTERED EXPLORATORY ANALYSIS: external validation against measured f(RH).

Protocol: at sites with long-term humidified-nephelometer data co-located with or near AERONET (list from
Burgos et al., 2019), test whether the site-level alpha_f slope is negatively correlated with measured
f(RH=85%).

Specification fixed before computing any slope at the new sites:
  * f(RH): Burgos et al. (2019, Sci. Data 6, 157), Table 5, f(RH=85%/RHdry), PM10 or whole air, median (P50),
    total scattering at 550 nm (450/525/635 nm Ecotech sites per their Table 2).
  * AERONET match: nearest AERONET site within 50 km that has >= 300 days of Level 2.0 SDA data in 2015-2024
    (same threshold as the confirmation eligibility rule). Jungfraujoch, Chebogue Point, Black Forest,
    Finokalia, Shouxian, Kaashidhoo, Lin'an, Mace Head, Niamey, Nainital, Cape Cod, Point Reyes and Yosemite
    have no such site (or no f(RH) in the PM10 column).
  * Slope: frozen primary model of confirm_analysis.py (prereg-v2), 2015-2024, same thresholds.
  * Test: Spearman rank correlation between the alpha_f slope and f(RH=85%), one-sided (negative); Pearson
    and the FMF/AE slopes reported descriptively.

usage: python frh_validation.py FROZEN_SCRIPTS_DIR FRH_SITES_CSV OUT_CSV HOURLY_DIR [HOURLY_DIR ...]
"""
import sys
from pathlib import Path
import pandas as pd
from scipy import stats

sys.path.insert(0, sys.argv[1])
import confirm_analysis as ca


def main():
    sites, out, dirs = pd.read_csv(sys.argv[2]), sys.argv[3], [Path(p) for p in sys.argv[4:]]
    rows = []
    for s in sites.itertuples():
        f = next((d / f"{s.site}.csv" for d in dirs if (d / f"{s.site}.csv").exists()), None)
        if f is None:
            rows.append(dict(code=s.code, site=s.site, status="NO_HOURLY_FILE")); continue
        d = ca.load_site(f, f.parent, s.lon)       # RH_lag7 not needed; the hourly file serves as the RH source
        for var in ca.VARS:
            e = ca.keep_days(d.dropna(subset=[var]), ca.MIN_H)
            nd = e.day.nunique()
            if nd < ca.MIN_DAYS:
                rows.append(dict(code=s.code, site=s.site, var=var, n_days=nd, status="lt_min_days")); continue
            r = ca.fit(e, var, ["RH"], ["PW"], 10)
            rows.append(dict(code=s.code, site=s.site, var=var, f85=s.f85_pm10_p50, beta=r["beta"][0], se=r["se"][0],
                             lo=r["lo"][0], hi=r["hi"][0], n_days=r["n_days"], status="ok"))
        print("done", s.site, flush=True)
    R = pd.DataFrame(rows); R.to_csv(out, index=False)
    print(R.round(4).to_string(index=False))
    for var in ca.VARS:
        g = R[(R["var"] == var) & (R.status == "ok")]
        if len(g) >= 4:
            rs = stats.spearmanr(g.beta, g.f85); rp = stats.pearsonr(g.beta, g.f85)
            print(f"{var}: n={len(g)} Spearman rho={rs[0]:+.3f} (one-sided p for rho<0: {rs[1]/2 if rs[0] < 0 else 1 - rs[1]/2:.3f}); Pearson r={rp[0]:+.3f}")


if __name__ == "__main__":
    main()
