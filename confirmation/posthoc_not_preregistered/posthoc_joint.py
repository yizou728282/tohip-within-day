"""POST-HOC, NOT PREREGISTERED: current RH and lag-7-day RH in one within-day model (frozen fit/meta functions reused unchanged)."""
import sys
from pathlib import Path
import pandas as pd
sys.path.insert(0, sys.argv[1])
import confirm_analysis as ca
sites = ca.read_sites(sys.argv[2]); hourly, rhdir = Path(sys.argv[3]), sys.argv[4]
rows = []
for s in sites.itertuples():
    d = ca.load_site(hourly / f"{s.site_id}.csv", rhdir, s.lon)
    for var in ca.VARS:
        e = ca.keep_days(d.dropna(subset=[var, "RH_lag7"]), ca.MIN_H)
        if e.day.nunique() < ca.MIN_DAYS: continue
        r = ca.fit(e, var, ["RH", "RH_lag7"], ["PW"], 10)
        for i, nm in enumerate(["RH_now", "RH_lag7"]):
            rows.append(dict(site=s.site_id, cls=s.cls, var=var, term=nm, beta=r["beta"][i], se=r["se"][i], lo=r["lo"][i], hi=r["hi"][i]))
D = pd.DataFrame(rows); D.to_csv(sys.argv[5], index=False)
for (c, v, t), g in D.groupby(["cls", "var", "term"]):
    m = ca.meta(g); print(f"{c:10s} {v:8s} {t:8s} k={m['k']:2d} mu={m['mu']:+.5f} ci95=({m['ci95'][0]:+.5f},{m['ci95'][1]:+.5f}) nsig={m['n_sig']}")
