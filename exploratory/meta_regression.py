"""PREREGISTERED EXPLORATORY ANALYSIS: sources of heterogeneity (protocol: "以气候态 FMF、平均 RH、到海岸距离做 meta 回归").

Mixed-effects meta-regression of the primary site slopes (per +10 % RH) on climatological median FMF, mean RH
of the analysed hours, and log10(1 + distance to coast in km), all standardised (z-scores). tau^2 by REML
(Fisher scoring), Knapp-Hartung standard errors with k - p degrees of freedom. Run per variable across all
72 confirmation sites, and within the pollution class for the H1 variables.
Choices not fixed by the protocol (made before running): moderator scaling, log of coast distance, the
all-site and pollution-only scopes.

usage: python meta_regression.py SITE_RESULTS_CSV SITE_COVARIATES_CSV OUT_CSV
"""
import sys
import numpy as np, pandas as pd
from scipy import stats


def reml_metareg(y, v, X, it=200):
    k, p = X.shape
    t2 = max(0.0, np.var(y, ddof=1) - v.mean())
    for _ in range(it):
        W = np.diag(1 / (v + t2))
        A = np.linalg.inv(X.T @ W @ X)
        P = W - W @ X @ A @ X.T @ W
        new = max(0.0, t2 + (y @ P @ P @ y - np.trace(P)) / np.trace(P @ P))
        if abs(new - t2) < 1e-14:
            t2 = new; break
        t2 = new
    W = np.diag(1 / (v + t2)); A = np.linalg.inv(X.T @ W @ X)
    b = A @ X.T @ W @ y
    s2 = float(((y - X @ b) ** 2 / (v + t2)).sum() / (k - p))
    se = np.sqrt(np.diag(A) * max(s2, 1e-300))
    return b, se, t2, k - p


def main():
    S, C, out = pd.read_csv(sys.argv[1]), pd.read_csv(sys.argv[2]), sys.argv[3]
    C["log_coast"] = np.log10(1 + C.dist_coast_km)
    mods = ["medFMF", "meanRH", "log_coast"]
    rows = []
    for scope in ["all", "pollution"]:
        for var in ["AE", "alpha_f", "FMF"]:
            if scope == "pollution" and var == "AE":
                continue
            d = S[(S.spec == "primary") & (S["var"] == var)].merge(C, on="site")
            if scope != "all":
                d = d[d["class"] == scope]
            Z = (d[mods] - d[mods].mean()) / d[mods].std(ddof=1)
            X = np.column_stack([np.ones(len(d)), Z.values])
            y, v = d.beta.values, d.se.values ** 2
            b, se, t2, df = reml_metareg(y, v, X)
            _, _, t2_0, _ = reml_metareg(y, v, X[:, :1])
            for name, bi, si in zip(["intercept"] + mods, b, se):
                t = bi / si
                rows.append(dict(scope=scope, var=var, k=len(d), term=name, coef_per_SD=bi, se=si, t=t, df=df,
                                 p=2 * stats.t.sf(abs(t), df), tau2=t2, tau2_null=t2_0,
                                 R2_tau=max(0.0, 1 - t2 / t2_0) if t2_0 > 0 else np.nan))
    R = pd.DataFrame(rows); R.to_csv(out, index=False)
    print(R.round(5).to_string(index=False))


if __name__ == "__main__":
    main()
