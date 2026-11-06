"""Derived from the frozen outputs (no code change): H1a-H3 recomputed under each preregistered sensitivity
specification, as the protocol requires ("敏感性分析（对 H1a–H3 逐一重算）"). The frozen script writes the pooled
estimates per specification (meta_results.csv) but evaluates the hypotheses only for the primary specification.
Same rules as the primary test: one-sided HKSJ t-tests, H3 with SE = sqrt(SE1^2 + SE2^2) and min(k)-1 df,
Holm over the four tests within each specification.

usage: python sensitivity_hypotheses.py META_RESULTS_CSV OUT_CSV
"""
import sys
import numpy as np, pandas as pd
from scipy import stats
M = pd.read_csv(sys.argv[1])
rows = []
for spec in ["primary", "sens_lnAOD", "sens_xk", "sens_sza70", "sens_5h", "sens_aodjump"]:
    g = lambda c, v: M[(M.cls == c) & (M["var"] == v) & (M.spec == spec)].iloc[0]
    H = {}
    a, b, f, m = g("pollution", "alpha_f"), g("pollution", "FMF"), g("marine", "AE"), g("pollution", "AE")
    H["H1a"] = (a.k, a.mu, a.p_less); H["H1b"] = (b.k, b.mu, b.p_greater); H["H2"] = (f.k, f.mu, f.p_less)
    diff, se, df = m.mu - f.mu, float(np.hypot(m.se, f.se)), min(m.k, f.k) - 1
    H["H3"] = (f"{m.k}/{f.k}", diff, stats.t.sf(diff / se, df))
    order = sorted(H, key=lambda h: H[h][2]); run = 0.0
    for i, h in enumerate(order):
        run = max(run, min(1.0, (len(order) - i) * H[h][2]))
        rows.append(dict(spec=spec, test=h, k=H[h][0], estimate=H[h][1], p=H[h][2], p_holm=run, supported=run < 0.05))
R = pd.DataFrame(rows).sort_values(["spec", "test"]); R.to_csv(sys.argv[2], index=False)
print(R.to_string(index=False, float_format=lambda x: f"{x:.4g}"))
