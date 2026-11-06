import sys, itertools
import numpy as np, pandas as pd
from scipy import stats
df = pd.read_csv(sys.argv[1], parse_dates=["date"])
order = ["marine","urban","smoke","dust"]
def resid(y, X):
    X = np.column_stack([np.ones(len(y))] + X)
    b, *_ = np.linalg.lstsq(X, y, rcond=None); return y - X @ b
rng = np.random.default_rng(0)
rows = []
for (site, reg), g in df.groupby(["site_id","regime"]):
    for v in ["AE","FMF"]:
        d = g.dropna(subset=[v,"RH","PW_cm","AOD"]); d = d[d.AOD>0]
        if len(d) < 30: continue
        if len(d) > 1500: d = d.sample(1500, random_state=0)
        y, rh, pw, la = d[v].values, d.RH.values.astype(float), d.PW_cm.values, np.log(d.AOD.values)
        out = dict(site=site, regime=reg, var=v, n=len(d), r_PW_RH=stats.spearmanr(pw, rh)[0])
        variants = {"raw": ([],[]), "PW": ([pw],[pw]), "PW+lnAOD": ([pw,la],[pw]),
                    "PW+lnAOD+season": None}
        for name, ctl in variants.items():
            if name == "PW+lnAOD+season":
                m = d.date.dt.month.values
                S = [np.sin(2*np.pi*m/12), np.cos(2*np.pi*m/12)]
                yr, xr = resid(y,[pw,la]+S), resid(rh,[pw]+S)
            elif not ctl[0]: yr, xr = y - y.mean(), rh - rh.mean()
            else: yr, xr = resid(y, ctl[0]), resid(rh, ctl[1])
            out[f"rho_{name}"] = stats.spearmanr(yr, xr)[0]
            out[f"g_{name}"] = stats.theilslopes(yr, xr)[0]
            out[f"sdRH_{name}"] = xr.std()
        out["dY_1sdRH_PW+lnAOD"] = out["g_PW+lnAOD"] * out["sdRH_PW+lnAOD"]
        out["sdY"] = y.std()
        rows.append(out)
R = pd.DataFrame(rows)
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
for v in ["AE","FMF"]:
    s = R[R["var"]==v]
    print(f"\n=== {v}: Spearman(Y*,RH*) by control set ===")
    print(s[["site","regime","n","r_PW_RH","rho_raw","rho_PW","rho_PW+lnAOD","rho_PW+lnAOD+season","dY_1sdRH_PW+lnAOD","sdY"]].round(3).to_string(index=False))
    for c in ["rho_raw","rho_PW","rho_PW+lnAOD","rho_PW+lnAOD+season"]:
        print(f"  median |{c}| = {s[c].abs().median():.3f}")
# order test: exact chance + label permutation
def frac(med, sign):
    S = {k: sign*med[k] for k in order}
    pairs = list(itertools.combinations(order,2))
    return sum(S[a] >= S[b] for a,b in pairs)/len(pairs)
perms = list(itertools.permutations(range(4)))
inv = [sum(p[i]>p[j] for i in range(4) for j in range(i+1,4)) for p in perms]
for k in range(7):
    print(f"P(random order fraction >= {k}/6) = {np.mean([(6-i)>=k for i in inv]):.3f}")
for v, sign, col in [("AE",-1,"g_PW+lnAOD"),("AE",-1,"g_PW"),("FMF",1,"g_PW+lnAOD"),("FMF",1,"g_PW")]:
    s = R[R["var"]==v].copy()
    obs = frac(s.groupby("regime")[col].median(), sign)
    null = []
    for _ in range(5000):
        s["perm"] = rng.permutation(s.regime.values)
        null.append(frac(s.groupby("perm")[col].median(), sign))
    print(f"{v} {col}: observed {obs:.2f}, label-permutation p(>=obs) = {np.mean(np.array(null)>=obs-1e-9):.3f}")
R.to_csv(sys.argv[2], index=False)
